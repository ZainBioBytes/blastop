#!/usr/bin/env python3
"""
blastop.py
-----------------------
Filters BLAST output (custom fmt 6) and returns ONE best hit per contig.

Expected columns (0-indexed):
  0  qseqid        query contig name
  1  qlen          query length
  2  sseqid        subject accession
  3  slen          subject length
  4  pident        % identity
  5  length        alignment length (HSP)
  6  evalue
  7  qstart
  8  qend
  9  sstart
  10 send
  11 (gap/mismatch or N/A)
  12 stitle        subject title
  13 taxid
  14 last col

Filters applied:
  - pident >= MIN_IDENTITY
  - cumulative query coverage >= MIN_COVERAGE
    (coverage = sum of non-overlapping HSP lengths / qlen, per contig+subject)
  - (optional) stitle must contain at least one of the taxa in --taxa

Best hit selection:
  - Highest pident among passing hits
  - Tiebreak: longest cumulative alignment length

Usage:
  python blastop.py -i input.tsv -o besthits.tsv
  python blastop.py -i input.tsv -o besthits.tsv --min_identity 95 --min_coverage 70
  python blastop.py -i input.tsv -o besthits.tsv --taxa "Bacteroides fragilis"
  python blastop.py -i input.tsv -o besthits.tsv --taxa "Bacteroides,Phocaeicola,Campylobacter"
"""

import argparse
import sys
from collections import defaultdict


MIN_IDENTITY = 80.0   # %
MIN_COVERAGE = 50.0   # % of query length


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Filter BLAST hits by identity/coverage thresholds and select the single\n"
            "best hit per contig (highest %%identity, tiebreak: longest alignment).\n\n"
            "Optionally restrict output to hits matching specific taxa using --taxa."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Basic run with default thresholds (80%% identity, 50%% coverage)\n"
            "  python blastop.py -i blast.tsv -o besthits.tsv\n\n"
            "  # Stricter identity threshold\n"
            "  python blastop.py -i blast.tsv -o besthits.tsv --min_identity 95\n\n"
            "  # Filter to a single taxon (partial match, case-insensitive)\n"
            "  python blastop.py -i blast.tsv -o besthits.tsv --taxa 'Bacteroides fragilis'\n\n"
            "  # Filter to multiple taxa (comma-separated, no spaces around commas)\n"
            "  python blastop.py -i blast.tsv -o besthits.tsv \\\n"
            "      --taxa 'Bacteroides,Phocaeicola,Campylobacter' --min_identity 90\n"
        )
    )
    p.add_argument(
        "-i", "--input", required=True, metavar="FILE",
        help="Input BLAST TSV file (custom outfmt 6 with 15 columns)"
    )
    p.add_argument(
        "-o", "--output", required=True, metavar="FILE",
        help="Output TSV file for best hits"
    )
    p.add_argument(
        "--min_identity", type=float, default=MIN_IDENTITY, metavar="FLOAT",
        help=f"Minimum %% identity threshold (default: {MIN_IDENTITY})"
    )
    p.add_argument(
        "--min_coverage", type=float, default=MIN_COVERAGE, metavar="FLOAT",
        help=f"Minimum query coverage %% (default: {MIN_COVERAGE}). "
             "Computed as cumulative non-overlapping HSP length / query length."
    )
    p.add_argument(
        "--taxa", type=str, default=None, metavar="STR",
        help=(
            "Restrict hits to subjects whose title (stitle) OR accession (sseqid) "
            "contains any of the specified taxa/patterns. Provide a single name or a "
            "comma-separated list (no spaces around commas). Matching is "
            "case-insensitive substring search against both fields. "
            "Example: --taxa 'Bacteroides fragilis' or --taxa 'Bacteroides,Phocaeicola' "
            "or --taxa 'CP103158' to match by accession."
        )
    )
    return p.parse_args()


def merge_intervals(intervals):
    """Return total non-overlapping length of a list of (start, end) intervals."""
    if not intervals:
        return 0
    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return sum(e - s + 1 for s, e in merged)


def taxa_match(stitle, sseqid, taxa_list):
    """Return True if stitle OR sseqid contains any of the taxa (case-insensitive substring)."""
    stitle_lower = stitle.lower()
    sseqid_lower = sseqid.lower()
    return any(t in stitle_lower or t in sseqid_lower for t in taxa_list)


def main():
    args = parse_args()

    # Parse --taxa into a lowercase list for matching
    taxa_filter = None
    if args.taxa:
        taxa_filter = [t.strip().lower() for t in args.taxa.split(",") if t.strip()]
        print(f"Taxa filter active: {taxa_filter}")

    # ------------------------------------------------------------------ #
    # Pass 1: collect all HSPs per (contig, subject) pair                  #
    # ------------------------------------------------------------------ #
    # Structure:
    #   data[(qseqid, sseqid)] = {
    #       "qlen": int,
    #       "pident_max": float,
    #       "hsps": [(qstart, qend, pident, aln_len)],
    #       "stitle": str,
    #       "evalue": str,
    #       "slen": int,
    #       "taxid": str,
    #       "last": str,
    #   }
    data = defaultdict(lambda: {
        "qlen": 0, "pident_max": 0.0,
        "hsps": [], "stitle": "", "evalue": "",
        "slen": 0, "taxid": "", "last": ""
    })

    with open(args.input) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 15:
                continue
            try:
                qseqid  = cols[0]
                qlen    = int(cols[1])
                sseqid  = cols[2]
                slen    = int(cols[3])
                pident  = float(cols[4])
                aln_len = int(cols[5])
                evalue  = cols[6]
                qstart  = int(cols[7])
                qend    = int(cols[8])
                stitle  = cols[12]
                taxid   = cols[13]
                last    = cols[14]
            except (ValueError, IndexError):
                continue

            # Apply taxa filter early — skip HSPs that don't match
            if taxa_filter and not taxa_match(stitle, sseqid, taxa_filter):
                continue

            key = (qseqid, sseqid)
            d = data[key]
            d["qlen"]      = qlen
            d["slen"]      = slen
            d["stitle"]    = stitle
            d["evalue"]    = evalue
            d["taxid"]     = taxid
            d["last"]      = last
            d["hsps"].append((qstart, qend, pident, aln_len))
            if pident > d["pident_max"]:
                d["pident_max"] = pident

    # ------------------------------------------------------------------ #
    # Pass 2: compute coverage, apply filters, pick best hit per contig   #
    # ------------------------------------------------------------------ #
    # best_hit[qseqid] = (pident, cum_aln_len, key, coverage)
    best_hit = {}

    for (qseqid, sseqid), d in data.items():
        pident_max = d["pident_max"]
        if pident_max < args.min_identity:
            continue

        qlen        = d["qlen"]
        intervals   = [(qs, qe) for qs, qe, *_ in d["hsps"]]
        cum_aln_len = merge_intervals(intervals)
        coverage    = 100.0 * cum_aln_len / qlen if qlen > 0 else 0.0

        if coverage < args.min_coverage:
            continue

        if qseqid not in best_hit or \
           (pident_max, cum_aln_len) > (best_hit[qseqid][0], best_hit[qseqid][1]):
            best_hit[qseqid] = (pident_max, cum_aln_len, (qseqid, sseqid), coverage)

    # ------------------------------------------------------------------ #
    # Write output                                                         #
    # ------------------------------------------------------------------ #
    header = "\t".join([
        "qseqid", "qlen", "sseqid", "slen",
        "pident_max", "cum_aln_len", "query_coverage_pct",
        "evalue", "stitle", "taxid"
    ])

    passed  = 0
    skipped = 0

    with open(args.output, "w") as out:
        out.write(header + "\n")
        for qseqid in sorted(best_hit.keys(),
                             key=lambda x: int(x.replace("contig_", "")) if x.replace("contig_", "").isdigit() else x):
            pident_max, cum_aln_len, key, coverage = best_hit[qseqid]
            d = data[key]
            _, sseqid = key
            row = "\t".join([
                qseqid,
                str(d["qlen"]),
                sseqid,
                str(d["slen"]),
                f"{pident_max:.3f}",
                str(cum_aln_len),
                f"{coverage:.2f}",
                d["evalue"],
                d["stitle"],
                d["taxid"],
            ])
            out.write(row + "\n")
            passed += 1

    total_contigs = len(set(k[0] for k in data.keys()))
    skipped = total_contigs - passed

    print(f"Done.")
    print(f"  Total unique contigs in input : {total_contigs}")
    print(f"  Contigs passing filters       : {passed}")
    print(f"  Contigs filtered out          : {skipped}")
    print(f"  Filters: identity >= {args.min_identity}%, coverage >= {args.min_coverage}%")
    if taxa_filter:
        print(f"  Taxa filter: {args.taxa}")
    print(f"  Output written to: {args.output}")


if __name__ == "__main__":
    main()
