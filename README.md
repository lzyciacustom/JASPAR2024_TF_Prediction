# TFBS-Mouse-mm39

Predict transcription factor binding sites (TFBS) in mouse gene promoters using a **local JASPAR2024 database** and the mm39 genome — no internet connection required after setup.

---

## Why this script?

The JASPAR web API is frequently unreliable. This pipeline runs entirely offline using a local copy of the JASPAR2024 SQLite database, making it suitable for HPC environments and reproducible analyses.

**What it does:**
- Fetches promoter sequence for any mouse gene (Ensembl ID or Entrez ID) from BSgenome mm39
- Scans the promoter against all JASPAR2024 CORE mouse PWMs
- Outputs a full site table, a per-TF summary, and three publication-ready plots

---

## Output files

| File | Description |
|------|-------------|
| `TFBS_all_sites.csv` | All predicted binding sites with genomic coordinates and PWM scores |
| `TF_summary.csv` | One row per TF — site count, max/mean score, coordinate range |
| `TFBS_filtered_sites.csv` | Sites matching your TF filter pattern (e.g. HIF-1, FOXO1) |
| `all_tfbs_plot.pdf` | Segment map of all TFs across the promoter |
| `filtered_tfbs_plot.pdf` | Segment map for filtered TFs only |
| `top10_tfbs_plot.pdf` | Position map + bar charts for top 10 TFs by PWM score |

---

## Requirements

### R packages

```r
install.packages(c("ggplot2", "dplyr", "gridExtra"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "biomaRt",
  "TFBSTools",
  "JASPAR2024",
  "GenomicRanges",
  "BSgenome.Mmusculus.UCSC.mm39",
  "Biostrings"
))

install.packages(c("DBI", "RSQLite"))
```

### JASPAR2024 local database

The script requires a local copy of the JASPAR2024 SQLite database (~500 MB).

Download from the JASPAR FTP:

```bash
wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024.sqlite3
```

Set the path in `config.R` (see below).

---

## Usage

### 1. Clone the repo

```bash
git clone https://github.com/YOUR_USERNAME/TFBS-Mouse-mm39.git
cd TFBS-Mouse-mm39
```

### 2. Edit `config.R`

All user-facing parameters are in one place:

```r
ENTREZ_ID      <- "208188"           # Entrez gene ID
UPSTREAM       <- 2000               # bp upstream of TSS to define promoter
JASPAR_DB_PATH <- "/path/to/JASPAR2024.sqlite3"
OUTDIR         <- "."                # where to save outputs
MIN_SCORE      <- "85%"             # PWM match threshold (85-95% recommended)
TF_FILTER      <- "HIF-1|FOXO1|Nrf1"  # regex for filtered plot
```

### 3. Run

```r
source("predict_tfbs.R")
```

Or from the command line:

```bash
Rscript predict_tfbs.R
```

---

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `ENTREZ_ID` | `"208188"` | Mouse Entrez gene ID |
| `UPSTREAM` | `2000` | Promoter window size (bp) |
| `MIN_SCORE` | `"85%"` | Lower = more hits, more false positives. Raise to 90–95% for stricter results |
| `TF_FILTER` | `"HIF-1\|FOXO1\|Nrf1"` | Regex pattern; case-insensitive |
| `JASPAR_DB_PATH` | — | Full path to local `.sqlite3` file |
| `OUTDIR` | `"."` | Output directory |

---

## Example output

Running on gene **Entrez 208188** with a 2000 bp upstream window and 85% score threshold:

- ~1000+ predicted binding sites across ~100 TFs
- Top TFs ranked by max PWM score in `TF_summary.csv`
- Three PDF plots generated automatically

---

## Limitations

- Biomart query requires an internet connection at runtime (for coordinate lookup). To run fully offline, replace with a pre-saved gene coordinate table.
- PWM scanning at 85% produces false positives; results should be validated against published ChIP-seq data where available.
- Currently configured for **mouse (mm39)** only. Human adaptation requires swapping `BSgenome.Mmusculus.UCSC.mm39` and the JASPAR species filter.

---

## Repo structure

```
TFBS-Mouse-mm39/
├── README.md
├── config.R            # all user parameters — edit this
├── predict_tfbs.R      # main pipeline
└── example/
    └── TF_summary_example.csv   # example output for gene 208188
```

---

## Citation

If you use this pipeline, please cite:

- **JASPAR 2024**: Castro-Mondragon et al., *Nucleic Acids Research*, 2024
- **TFBSTools**: Tan & Lenhard, *Bioinformatics*, 2016
- **BSgenome / Bioconductor**: Huber et al., *Nature Methods*, 2015

---

## License

MIT
