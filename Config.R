# ==============================================================================
# config.R — Edit this file before running predict_tfbs.R
# ==============================================================================
 
# Gene to analyze (Entrez ID)
ENTREZ_ID      <- "208188"
 
# Promoter window: bp upstream of TSS
UPSTREAM       <- 2000
 
# Path to local JASPAR2024 SQLite database
# Download: https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024.sqlite3
JASPAR_DB_PATH <- "E:/JASPAR2024.sqlite3"
 
# Output directory (will be created if it doesn't exist)
OUTDIR         <- "."
 
# PWM match threshold: "85%" to "95%"
# Lower = more hits but more false positives
MIN_SCORE      <- "85%"
 
# Regex pattern for filtered plot (case-insensitive, use | for multiple)
TF_FILTER      <- "HIF-1|FOXO1|Nrf1"
 
