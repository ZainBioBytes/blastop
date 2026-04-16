# blastop

<p align="center">
  <img src="logo.svg" alt="blastop logo" width="600"/>
</p>

**BLAST Output Processor** — filter high-confidence hits and extract the single best hit per contig from BLAST tabular output.

Built for metagenomic/assembly workflows where BLAST produces multiple HSPs and hits per contig and you need a clean, filtered, one-hit-per-contig summary.

---

## Features

- Filters hits by **minimum % identity** and **minimum query coverage**
- Coverage is computed from **cumulative non-overlapping HSPs** (prevents double-counting)
- Selects **one best hit per contig** (highest % identity; tiebreak: longest alignment)
- Optional **taxa filter** — restrict output to specific organisms by name or accession (searches both `stitle` and `sseqid`, case-insensitive substring match)
- Supports **comma-separated multi-taxa filtering** in a single run

---

## Requirements

- Python 3.6+
- No external dependencies (standard library only)

---

## Usage

```bash
python blastop.py -i input.tsv -o besthits.tsv [Options]
```

### Options

```
-i / --input        Input BLAST TSV file (required)
-o / --output       Output TSV file (required)
--min_identity      Minimum % identity threshold (default: 80.0)
--min_coverage      Minimum query coverage % (default: 50.0)
--taxa              Restrict to taxa matching this string (substring, case-insensitive)
                    Comma-separated for multiple taxa. Searches stitle and sseqid.
```

### Examples

```bash
# Default thresholds (80% identity, 50% coverage)
python blastop.py -i blast.tsv -o besthits.tsv

# Stricter identity threshold
python blastop.py -i blast.tsv -o besthits.tsv --min_identity 95

# Filter to a single taxon
python blastop.py -i blast.tsv -o besthits.tsv --taxa "Campylobacter coli"

# Filter to multiple taxa (comma-separated, no spaces around commas)
python blastop.py -i blast.tsv -o besthits.tsv \
    --taxa "Bacteroides,Phocaeicola,Campylobacter" --min_identity 90

# Match by accession prefix
python blastop.py -i blast.tsv -o besthits.tsv --taxa "CP103"
```

---

## Input format

BLAST tabular output (`-outfmt 6`) with the following custom columns:

```
qseqid qlen sseqid slen pident length evalue qstart qend sstart send gaps stitle staxid score
```

Generate with:

```bash
blastn -query contigs.fasta -db nt \
    -outfmt "6 qseqid qlen sseqid slen pident length evalue qstart qend sstart send gaps stitle staxid score" \
    -out blast_results.tsv
```

---

## Output format

Tab-separated file with one row per contig (best hit only):

| Column | Description |
|---|---|
| qseqid | Query contig name |
| qlen | Query contig length (bp) |
| sseqid | Subject accession |
| slen | Subject sequence length |
| pident_max | Highest % identity across HSPs |
| cum_aln_len | Cumulative non-overlapping alignment length |
| query_coverage_pct | % of query covered by alignment |
| evalue | E-value of best HSP |
| stitle | Subject title |
| taxid | Taxonomy ID |

---

## License

MIT License. See [LICENSE](LICENSE).

---

## Citation

If you use blastop in your research, please cite this repository:

```
Zain Abedien. blastop: BLAST Output Processor. GitHub. https://github.com/ZainBioBytes/blastop
```
