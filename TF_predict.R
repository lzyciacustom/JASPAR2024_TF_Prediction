# ==============================================================================
# TFBS Prediction Pipeline
# Mouse genome mm39 | JASPAR2024 local SQLite
# ==============================================================================

library(biomaRt)
library(TFBSTools)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm39)
library(ggplot2)
library(DBI)
library(RSQLite)
library(dplyr)

# ==============================================================================
# SECTION 1: CONSTANTS - only edit this section
# ==============================================================================
GENE_ID      <- "ENSMUSG00000051136"
ENTREZ_ID    <- "208188"
UPSTREAM     <- 2000          
JASPAR_DB_PATH <- "E:/JASPAR2024.sqlite3"
OUTDIR       <- "."           
MIN_SCORE    <- "85%"         
# TFs to highlight in filtered plot (regex, case-insensitive)
TF_FILTER    <- "HIF-1|FOXO1|Nrf1"

# ==============================================================================
# SECTION 2: RETRIEVE GENE INFO AND PROMOTER SEQUENCE
# ==============================================================================
message("=== Retrieving gene information from Ensembl ===")

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_info <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand"),
  filters    = "entrezgene_id",
  values     = ENTREZ_ID,
  mart       = ensembl
)

if (nrow(gene_info) == 0) stop("Gene not found in Ensembl database.")

chr <- paste0("chr", gene_info$chromosome_name)

# Calculate promoter coordinates based on strand
if (gene_info$strand == 1) {
  promoter_start <- gene_info$start_position - 1 - UPSTREAM
  promoter_end   <- gene_info$start_position - 1
} else {
  promoter_start <- gene_info$end_position
  promoter_end   <- gene_info$end_position + UPSTREAM
}

message(sprintf("Gene: %s | Chr: %s | Promoter: %d-%d | Strand: %s",
                ENTREZ_ID, chr, promoter_start, promoter_end,
                ifelse(gene_info$strand == 1, "+", "-")))

# Create GRanges and fetch sequence
promoter_range <- GRanges(
  seqnames = chr,
  ranges   = IRanges(start = promoter_start, end = promoter_end),
  strand   = ifelse(gene_info$strand == 1, "+", "-")
)

promoter_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm39, promoter_range)
message("Promoter sequence length: ", width(promoter_seq), " bp")

# ==============================================================================
# SECTION 3: LOAD JASPAR MATRICES
# ==============================================================================
message("=== Loading JASPAR2024 matrices ===")

if (!file.exists(JASPAR_DB_PATH)) {
  stop("JASPAR2024 SQLite database not found: ", JASPAR_DB_PATH)
}

jaspar_db <- dbConnect(SQLite(), dbname = JASPAR_DB_PATH)

pfm_mouse <- tryCatch({
  getMatrixSet(jaspar_db, opts = list(species = "Mus musculus", collection = "CORE"))
}, error = function(e) {
  dbDisconnect(jaspar_db)
  stop("Failed to retrieve JASPAR matrices: ", e$message)
})

message("Loaded ", length(pfm_mouse), " mouse TF matrices from JASPAR2024")

# ==============================================================================
# SECTION 4: SCAN FOR TFBS
# ==============================================================================
message("=== Scanning promoter for TFBS (min.score = ", MIN_SCORE, ") ===")

site_results <- lapply(pfm_mouse, function(pfm) {
  pwm <- toPWM(pfm)
  searchSeq(pwm, promoter_seq, min.score = MIN_SCORE)
})

# Filter out empty results
site_results_filtered <- Filter(function(x) length(x) > 0, site_results)
message("TFs with at least one hit: ", length(site_results_filtered))

# ==============================================================================
# SECTION 5: EXTRACT RESULTS INTO DATA FRAME
# ==============================================================================
message("=== Extracting TFBS positions ===")

tfbs_positions <- do.call(rbind, lapply(site_results_filtered, function(site_set_list) {
  site_set <- site_set_list[[1]]
  if (length(site_set) > 0) {
    data.frame(
      start   = start(site_set@views),
      end     = end(site_set@views),
      score   = site_set@score,
      TF      = site_set@pattern@name,
      TF_ID   = site_set@pattern@ID,
      strand  = site_set@strand,
      seqname = site_set@seqname,
      stringsAsFactors = FALSE
    )
  }
}))

if (is.null(tfbs_positions) || nrow(tfbs_positions) == 0) {
  dbDisconnect(jaspar_db)
  stop("No TFBS found. Try lowering MIN_SCORE threshold.")
}

# Add absolute genomic coordinates
tfbs_positions$genomic_start <- promoter_start + tfbs_positions$start
tfbs_positions$genomic_end   <- promoter_start + tfbs_positions$end

message("Total TFBS hits: ", nrow(tfbs_positions))
message("Unique TFs: ", length(unique(tfbs_positions$TF)))

# ==============================================================================
# SECTION 6: SAVE CSV OUTPUTS
# ==============================================================================
message("=== Saving CSV outputs ===")

# Full TFBS table
write.csv(tfbs_positions,
          file = file.path(OUTDIR, "TFBS_all_sites.csv"),
          row.names = FALSE)

# TF summary (one row per TF)
tf_summary <- tfbs_positions %>%
  group_by(TF, TF_ID) %>%
  summarise(
    n_sites           = n(),
    max_score         = max(score),
    mean_score        = round(mean(score), 4),
    min_genomic_start = min(genomic_start),
    max_genomic_end   = max(genomic_end),
    strands_found     = paste(sort(unique(strand)), collapse = "/"),
    .groups = "drop"
  ) %>%
  arrange(desc(max_score))

write.csv(tf_summary,
          file = file.path(OUTDIR, "TF_summary.csv"),
          row.names = FALSE)

message("Saved: TFBS_all_sites.csv (", nrow(tfbs_positions), " rows)")
message("Saved: TF_summary.csv (", nrow(tf_summary), " unique TFs)")

# ==============================================================================
# SECTION 7: PLOT - ALL TFBS
# ==============================================================================
message("=== Generating plots ===")

n_tfs       <- length(unique(tfbs_positions$TF))
plot_height <- max(8, n_tfs * 0.25)  # dynamic height: ~0.25 inch per TF

pdf(file.path(OUTDIR, "all_tfbs_plot.pdf"), width = 14, height = plot_height)
print(
  ggplot(tfbs_positions,
         aes(x = genomic_start, xend = genomic_end,
             y = TF, yend = TF, color = score)) +
    geom_segment(linewidth = 2.5) +
    scale_color_gradient(low = "#91bfdb", high = "#d73027") +
    scale_x_continuous(
      limits = c(promoter_start, promoter_end),
      expand = c(0, 0)
    ) +
    theme_minimal(base_size = 10) +
    labs(
      title  = paste0("Predicted TFBS — All TFs (", n_tfs, " factors)"),
      x      = "Genomic Position (bp)",
      y      = "Transcription Factor",
      color  = "PWM Score"
    ) +
    theme(
      axis.text.y  = element_text(size = max(4, 8 - n_tfs %/% 50)),
      axis.title   = element_text(size = 11),
      plot.title   = element_text(size = 13, face = "bold"),
      panel.grid.major.y = element_line(color = "grey90")
    )
)
dev.off()

message("Saved: all_tfbs_plot.pdf")

# ==============================================================================
# SECTION 8: PLOT - FILTERED TFBS
# ==============================================================================
tfbs_filtered <- tfbs_positions %>%
  filter(grepl(TF_FILTER, TF, ignore.case = TRUE)) %>%
  filter(genomic_start >= promoter_start & genomic_end <= promoter_end)

message("Filtered TFs matching '", TF_FILTER, "': ",
        length(unique(tfbs_filtered$TF)), " factors, ",
        nrow(tfbs_filtered), " sites")
print(unique(tfbs_filtered$TF))

write.csv(tfbs_filtered,
          file = file.path(OUTDIR, "TFBS_filtered_sites.csv"),
          row.names = FALSE)

if (nrow(tfbs_filtered) == 0) {
  message("No binding sites found for TF_FILTER pattern. Skipping filtered plot.")
} else {
  n_tfs_filt <- length(unique(tfbs_filtered$TF))
  
  pdf(file.path(OUTDIR, "filtered_tfbs_plot.pdf"),
      width = 12, height = max(5, n_tfs_filt * 0.5 + 2))
  print(
    ggplot(tfbs_filtered,
           aes(x = genomic_start, xend = genomic_end,
               y = TF, yend = TF, color = score)) +
      geom_segment(linewidth = 3) +
      scale_color_gradient(low = "#91bfdb", high = "#d73027") +
      scale_x_continuous(
        limits = c(promoter_start, promoter_end),
        expand = c(0.02, 0)
      ) +
      theme_minimal(base_size = 11) +
      labs(
        title  = paste0("Predicted TFBS — ", TF_FILTER),
        x      = "Genomic Position (bp)",
        y      = "Transcription Factor",
        color  = "PWM Score"
      ) +
      theme(
        axis.text.y = element_text(size = 9),
        axis.title  = element_text(size = 11),
        plot.title  = element_text(size = 13, face = "bold")
      )
  )
  dev.off()
  
  message("Saved: filtered_tfbs_plot.pdf")
}

# ==============================================================================
# SECTION 9: CLEANUP
# ==============================================================================
dbDisconnect(jaspar_db)
