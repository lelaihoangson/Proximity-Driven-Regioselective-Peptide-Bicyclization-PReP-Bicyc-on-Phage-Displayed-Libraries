#!/usr/bin/env Rscript

# ==============================================================================
# Supplementary Scripts: Peptide Library Amino-Acid and Venn-Diagram Analyses
# ==============================================================================
# Suggested filename:
#   supplementary_peptide_library_analysis.R
#
# Contents:
#   1. Amino-acid analysis of paired-end FASTQ reads
#   2. Venn diagram comparing peptide sets from selection rounds R1-R4
#
# Before running:
#   - Edit the USER CONFIGURATION section below.
#   - Install required packages if needed.
# ============================================================================== 

# ---- Required packages --------------------------------------------------------
required_packages <- c(
  "microseq",
  "RColorBrewer",
  "dplyr",
  "stringr",
  "gplots",
  "ggVennDiagram",
  "ggplot2",
  "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(microseq)
  library(RColorBrewer)
  library(dplyr)
  library(stringr)
  library(gplots)
  library(ggVennDiagram)
  library(ggplot2)
  library(grid)
})

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================

# Input FASTQ files for amino-acid analysis
forward_fastq <- "/Users/gqjin_mac/Documents/Sequencing/20240724/20240808/24355Wns_JinPD1R1_S1_L001_R1_001.fastq"
reverse_fastq <- "/Users/gqjin_mac/Documents/Sequencing/20240724/20240808/24355Wns_JinPD1R1_S1_L001_R2_001.fastq"

# Output directory
output_dir <- "/Users/gqjin_mac/Documents/Sequencing/20240724/PD1 analysis"

# Prefix used for output files
sample_name <- "PD1R1"

# Library definition
library_regex <- "GCCCAG.{45}GCGGCG.{6}"
library_start_nt <- 19L
library_codons <- 11L

# Required peptide pattern: CX4CX4C
required_cysteine_positions <- c(1L, 6L, 11L)

# Input peptide files for Venn-diagram analysis
venn_files <- c(
  R1 = file.path(output_dir, "AAsPD1R1_weblogo.txt"),
  R2 = file.path(output_dir, "AAsPD1R2_weblogo.txt"),
  R3 = file.path(output_dir, "AAsPD1R3_weblogo.txt"),
  R4 = file.path(output_dir, "AAsPD1R4_weblogo.txt")
)

# Set to FALSE to skip either analysis
run_amino_acid_analysis <- TRUE
run_venn_analysis <- TRUE

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

check_input_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "The following input files were not found:\n",
      paste(" -", missing, collapse = "\n")
    )
  }
}

translate_codon <- function(codon) {
  genetic_code <- c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    TAT = "Y", TAC = "Y", TAA = NA,  TAG = "TAG",
    TGT = "C", TGC = "C", TGA = NA,  TGG = "W",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )

  unname(genetic_code[toupper(codon)])
}

translate_library_region <- function(sequence, start_nt, n_codons) {
  end_nt <- start_nt + n_codons * 3L - 1L

  if (nchar(sequence) < end_nt) {
    return(rep(NA_character_, n_codons))
  }

  region <- substr(sequence, start_nt, end_nt)
  starts <- seq.int(1L, nchar(region), by = 3L)
  codons <- substring(region, starts, starts + 2L)
  vapply(codons, translate_codon, FUN.VALUE = character(1), USE.NAMES = FALSE)
}

has_required_cysteines <- function(row, positions) {
  peptide <- as.character(row)
  !any(is.na(peptide)) && all(peptide[positions] == "C")
}

save_heatmap_png <- function(
    matrix_data,
    filename,
    breaks,
    colors,
    xlab = "Position in Library",
    ylab = "Amino Acid",
    width = 1800,
    height = 1400,
    res = 200) {

  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)

  heatmap.2(
    matrix_data,
    Rowv = NA,
    Colv = NA,
    col = colors,
    density.info = "none",
    scale = "none",
    trace = "none",
    breaks = breaks,
    xlab = xlab,
    ylab = ylab,
    margins = c(6, 6),
    dendrogram = "none"
  )
}

# ==============================================================================
# SUPPLEMENTARY SCRIPT 1: AMINO-ACID ANALYSIS
# ==============================================================================

run_amino_acid_workflow <- function() {
  check_input_files(c(forward_fastq, reverse_fastq))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Reading paired-end FASTQ files...")
  forward_reads <- readFastq(forward_fastq)
  reverse_reads <- readFastq(reverse_fastq)

  forward_sequences <- forward_reads[[2]]
  reverse_sequences <- reverseComplement(reverse_reads[[2]], reverse = TRUE)

  message("Extracting sequences matching the library pattern...")
  forward_matches <- gregexpr(library_regex, forward_sequences, extract = TRUE)
  reverse_matches <- gregexpr(library_regex, reverse_sequences, extract = TRUE)

  n_pairs <- min(length(forward_matches), length(reverse_matches))
  accepted_sequences <- character(n_pairs)

  library_end_nt <- library_start_nt + library_codons * 3L - 1L

  for (i in seq_len(n_pairs)) {
    fwd <- forward_matches[[i]]
    rev <- reverse_matches[[i]]

    if (length(fwd) == 0L || length(rev) == 0L || is.na(fwd[1]) || is.na(rev[1])) {
      next
    }

    fwd <- fwd[1]
    rev <- rev[1]

    if (identical(fwd, rev)) {
      accepted_sequences[i] <- fwd
      next
    }

    fwd_chars <- strsplit(fwd, "", fixed = TRUE)[[1]]
    rev_chars <- strsplit(rev, "", fixed = TRUE)[[1]]

    if (length(fwd_chars) != length(rev_chars)) {
      next
    }

    mismatch_positions <- which(fwd_chars != rev_chars)

    # Accept exactly one mismatch only when it lies outside the library region.
    if (
      length(mismatch_positions) == 1L &&
      (mismatch_positions < library_start_nt || mismatch_positions > library_end_nt)
    ) {
      accepted_sequences[i] <- fwd
    }
  }

  accepted_sequences <- accepted_sequences[nzchar(accepted_sequences)]

  if (length(accepted_sequences) == 0L) {
    stop("No concordant paired-end library sequences passed the filtering criteria.")
  }

  message("Accepted paired reads: ", length(accepted_sequences))

  # Codon-spaced sequences, retained for compatibility with the original analysis.
  codon_strings <- gsub("(...)", "\\1 ", accepted_sequences)
  codon_strings <- trimws(codon_strings)

  sequence_counts <- as.data.frame(
    sort(table(codon_strings), decreasing = TRUE),
    stringsAsFactors = FALSE
  )
  colnames(sequence_counts) <- c("CodonSequence", "Count")

  amino_acid_matrix <- t(vapply(
    accepted_sequences,
    translate_library_region,
    start_nt = library_start_nt,
    n_codons = library_codons,
    FUN.VALUE = character(library_codons)
  ))

  amino_acid_df <- as.data.frame(amino_acid_matrix, stringsAsFactors = FALSE)
  colnames(amino_acid_df) <- paste0("Position_", seq_len(library_codons))

  # Count unique amino-acid sequences.
  unique_amino_acids <- amino_acid_df %>%
    count(across(everything()), name = "n", sort = TRUE) %>%
    filter(if_all(-n, ~ !is.na(.) & . != "0"))

  # Retain only CX4CX4C peptides: cysteine at positions 1, 6, and 11.
  cx4cx4c_mask <- apply(
    unique_amino_acids[, seq_len(library_codons), drop = FALSE],
    1,
    has_required_cysteines,
    positions = required_cysteine_positions
  )
  unique_cx4cx4c <- unique_amino_acids[cx4cx4c_mask, , drop = FALSE]

  # Percentage of reads containing at least one TAG codon.
  contains_tag <- apply(amino_acid_df, 1, function(x) any(x == "TAG", na.rm = TRUE))
  percent_tag <- mean(contains_tag) * 100

  # Remove all reads containing TAG.
  amino_acid_no_tag <- amino_acid_df[!contains_tag, , drop = FALSE]

  aa_levels <- c(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "TAG"
  )

  aa_table <- apply(
    amino_acid_df,
    2,
    function(x) table(factor(x, levels = aa_levels))
  )
  aa_table <- as.matrix(aa_table / nrow(amino_acid_df))
  colnames(aa_table) <- seq_len(library_codons)

  random_aa <- matrix(
    0,
    nrow = length(aa_levels),
    ncol = library_codons,
    dimnames = list(aa_levels, seq_len(library_codons))
  )
  random_aa[c("A", "G", "P", "T", "V"), ] <- 2 / 32
  random_aa[c("C", "H", "Q", "N", "K", "Y", "D", "E", "W", "I", "M", "TAG", "F"), ] <- 1 / 32
  random_aa[c("L", "S", "R"), ] <- 3 / 32

  library_bias <- (aa_table - random_aa) / random_aa

  aa_no_tag_table <- apply(
    amino_acid_no_tag,
    2,
    function(x) table(factor(x, levels = aa_levels))
  )
  aa_no_tag_table <- as.matrix(aa_no_tag_table / nrow(amino_acid_no_tag))
  colnames(aa_no_tag_table) <- seq_len(library_codons)

  # Write tables.
  write.csv(
    sequence_counts,
    file.path(output_dir, paste0(sample_name, "_codon_sequence_counts.csv")),
    row.names = FALSE
  )
  write.csv(
    unique_amino_acids,
    file.path(output_dir, paste0(sample_name, "_unique_amino_acid_sequences.csv")),
    row.names = FALSE
  )
  write.csv(
    unique_cx4cx4c,
    file.path(output_dir, paste0(sample_name, "_CX4CX4C_sequences.csv")),
    row.names = FALSE
  )
  write.csv(
    amino_acid_df,
    file.path(output_dir, paste0("AAs", sample_name, ".csv")),
    row.names = TRUE
  )
  write.csv(
    aa_table,
    file.path(output_dir, paste0("AAtable", sample_name, ".csv")),
    row.names = TRUE
  )
  write.csv(
    as.data.frame(library_bias),
    file.path(output_dir, paste0("librarybias", sample_name, ".csv")),
    row.names = TRUE
  )

  # Save heatmaps.
  heatmap_colors <- colorRampPalette(brewer.pal(9, "Blues"))(100)

  save_heatmap_png(
    aa_table,
    file.path(output_dir, paste0(sample_name, "_amino_acid_frequency_heatmap.png")),
    breaks = seq(0.0, 0.6, length.out = 101),
    colors = heatmap_colors
  )

  save_heatmap_png(
    random_aa,
    file.path(output_dir, "NNK_expected_amino_acid_frequency_heatmap.png"),
    breaks = seq(0.0, 0.6, length.out = 101),
    colors = heatmap_colors
  )

  save_heatmap_png(
    library_bias,
    file.path(output_dir, paste0(sample_name, "_library_bias_heatmap.png")),
    breaks = seq(-1, 4, length.out = 101),
    colors = heatmap_colors,
    ylab = "Amino Acid"
  )

  save_heatmap_png(
    aa_no_tag_table,
    file.path(output_dir, paste0(sample_name, "_amino_acid_frequency_no_TAG_heatmap.png")),
    breaks = seq(0.0, 0.3, length.out = 101),
    colors = heatmap_colors
  )

  summary_file <- file.path(output_dir, paste0(sample_name, "_analysis_summary.txt"))
  writeLines(
    c(
      paste0("Accepted paired reads: ", length(accepted_sequences)),
      paste0("Unique amino-acid sequences: ", nrow(unique_amino_acids)),
      paste0("Unique CX4CX4C sequences: ", nrow(unique_cx4cx4c)),
      paste0("Reads containing TAG (%): ", round(percent_tag, 4))
    ),
    summary_file
  )

  message("Amino-acid analysis completed.")
  message("TAG-containing reads: ", round(percent_tag, 4), "%")
}

# ==============================================================================
# SUPPLEMENTARY SCRIPT 2: VENN DIAGRAM FOR R1-R4
# ==============================================================================

run_venn_workflow <- function() {
  check_input_files(venn_files)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  peptide_sets <- lapply(venn_files, function(filename) {
    data <- read.table(
      filename,
      header = FALSE,
      stringsAsFactors = FALSE,
      col.names = "Peptide",
      strip.white = TRUE,
      blank.lines.skip = TRUE
    )

    unique(data$Peptide[nzchar(data$Peptide) & !is.na(data$Peptide)])
  })

  custom_colors <- c(
    "#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#FB9A99",
    "#FFFF99", "#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A"
  )

  venn_plot <- ggVennDiagram(
    peptide_sets,
    label_alpha = 1,
    edge_size = 0.8,
    label = "count",
    category.names = names(peptide_sets)
  ) +
    scale_fill_gradientn(colors = custom_colors, name = "Peptide Count") +
    theme_void() +
    labs(title = "Venn Diagram of Peptides in R1, R2, R3, and R4") +
    theme(
      plot.title = element_text(
        size = 25,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 12)
      ),
      legend.position = "right",
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold"),
      legend.background = element_blank(),
      legend.key.size = unit(1.0, "cm")
    )

  print(venn_plot)

  venn_output <- file.path(output_dir, "AAsPD1R1-4_VennDiagram_Styled.png")
  ggsave(
    venn_output,
    plot = venn_plot,
    width = 10,
    height = 10,
    dpi = 600
  )

  message("Venn diagram saved to: ", venn_output)
}

# ==============================================================================
# RUN ANALYSES
# ==============================================================================

if (run_amino_acid_analysis) {
  run_amino_acid_workflow()
}

if (run_venn_analysis) {
  run_venn_workflow()
}

message("All requested supplementary analyses completed successfully.")
