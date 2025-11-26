# Set working directory (adjust to your actual path)
setwd("/path_to_file")

# Install required packages (uncomment on first run)
# install.packages(c("adegenet", "poppr", "ape", "ggplot2", "dplyr", "pheatmap"))

# Load required libraries
library(adegenet)
library(poppr)
library(ape)
library(ggplot2)
library(dplyr)
library(pheatmap)

# ----------------------------
# 1. Read and transpose genotype data
# ----------------------------
raw_tall <- read.csv("genotype of 38 InDel markers.csv", header = TRUE, row.names = 1, check.names = FALSE)
raw_data <- t(raw_tall)  # Dimensions: 117 accessions × 38 markers
cat("✅ Data loaded successfully. '0' is treated as a valid allele (no amplification state).\n")

# ----------------------------
# 2. Compute Simple Matching Coefficient (SMC) matrix
# ----------------------------
geno_mat <- as.matrix(raw_data)
n_var <- nrow(geno_mat)
n_loci <- ncol(geno_mat)

smc_matrix <- matrix(0, nrow = n_var, ncol = n_var)
rownames(smc_matrix) <- rownames(geno_mat)
colnames(smc_matrix) <- rownames(geno_mat)

for (i in 1:n_var) {
  for (j in i:n_var) {
    same <- sum(geno_mat[i, ] == geno_mat[j, ])
    sim <- same / n_loci
    smc_matrix[i, j] <- sim
    smc_matrix[j, i] <- sim
  }
}
cat("✅ Genetic similarity matrix (SMC) computed.\n")

# ----------------------------
# 3. Export lower triangle of SMC matrix (rounded to 2 decimal places)
# ----------------------------
output_matrix <- matrix("", nrow = nrow(smc_matrix), ncol = ncol(smc_matrix),
                        dimnames = dimnames(smc_matrix))
lower_idx <- lower.tri(smc_matrix, diag = TRUE)
output_matrix[lower_idx] <- sprintf("%.2f", smc_matrix[lower_idx])

write.csv(output_matrix, "rice_accessions_genetic_similarity_lower_triangle.csv", quote = FALSE)
cat("✅ Lower-triangle genetic similarity matrix (2 decimal places) saved.\n")
cat("\n--- Genetic similarity matrix (first 10 accessions) ---\n")
print(output_matrix[1:10, 1:10], quote = FALSE)

# ----------------------------
# 4. Construct genetic distance matrix
# ----------------------------
dist_matrix <- as.dist(1 - smc_matrix)

# ----------------------------
# 7. Generate digital fingerprints and detect duplicates
# ----------------------------
raw_data_char <- apply(geno_mat, c(1,2), as.character)
digital_fingerprint <- apply(raw_data_char, 1, function(x) paste(x, collapse = ""))

fingerprint_df <- data.frame(
  Variety_ID = names(digital_fingerprint),
  Digital_Fingerprint = digital_fingerprint,
  stringsAsFactors = FALSE
)

# Group by fingerprint to identify duplicates (frequency >= 2)
dup_groups <- fingerprint_df %>%
  group_by(Digital_Fingerprint) %>%
  filter(n() >= 2) %>%
  group_split()

if (length(dup_groups) == 0) {
  cat("✅ All 117 accessions have unique DNA fingerprints ('0' treated as valid state).\n")
} else {
  cat("⚠️ Duplicate fingerprints detected: ", length(dup_groups), " duplicate groups found:\n\n")
  
  # Prepare output lines
  output_lines <- c("⚠️ Duplicate fingerprint groups:\n")
  
  for (i in seq_along(dup_groups)) {
    group_varieties <- dup_groups[[i]]$Variety_ID
    group_str <- paste(group_varieties, collapse = ", ")
    msg <- paste0("Group ", i, ": ", group_str)
    cat(msg, "\n")
    output_lines <- c(output_lines, msg)
  }
  
  # Save to file
  writeLines(output_lines, "duplicate_fingerprint_groups.txt")
  cat("\n✅ Duplicate group details saved to 'duplicate_fingerprint_groups.txt'\n")
}

# Save full fingerprint table
write.csv(fingerprint_df, "rice_accessions_digital_fingerprints.csv", row.names = FALSE, quote = FALSE)

# Save as TXT to avoid Excel scientific notation issues
write.table(
  fingerprint_df,
  file = "rice_accessions_digital_fingerprints.txt",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  eol = "\n"
)
cat("✅ Digital fingerprints saved as TXT (to prevent Excel formatting issues).\n")

# Safely save CSV with quotes to preserve long numeric-like strings
write.table(
  fingerprint_df,
  file = "rice_accessions_digital_fingerprints_quoted.csv",
  sep = ",",
  row.names = FALSE,
  quote = TRUE,
  col.names = TRUE,
  eol = "\n"
)

# ----------------------------
# 8. Plot heatmap (export as TIFF + A4 PDF landscape)
# ----------------------------
heat_plot <- function() {
  if (!require("pheatmap")) stop("pheatmap package is required")
  if (!require("grid")) stop("grid package is required")
  
  pheatmap(
    mat = geno_mat,
    color = colorRampPalette(c("gray80", "green", "white", "red", "orange"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    main = "38 InDel Markers",
    fontsize = 8
  )
  
  # Add y-axis label on the left
  grid.text(
    label = "114 Shanlan upland rice landraces",
    x = unit(-1.5, "cm"),
    y = unit(0.5, "npc"),
    rot = 90,
    gp = gpar(fontsize = 10, fontfamily = "Times New Roman")
  )
}

# --- Export as high-resolution TIFF ---
tiff("Heatmap illustrating the genotypic variation of 114 Shanlan upland rice landraces.tiff", width = 3508, height = 2480, res = 300, compression = "lzw")
heat_plot()
dev.off()

# --- Export as A4 PDF (landscape orientation) ---
pdf("Heatmap illustrating the genotypic variation of 114 Shanlan upland rice landraces.pdf", width = 11.69, height = 8.27)
heat_plot()
dev.off()
cat("✅ Heatmap saved as TIFF and A4 PDF (landscape).\n")

# ----------------------------
# 9. Identify minimal marker set for unique fingerprinting
# ----------------------------
geno_mat_char <- apply(geno_mat, c(1,2), as.character)
n_var <- nrow(geno_mat_char)
selected_markers <- character(0)
current_fingerprints <- rep("", n_var)
remaining_markers <- colnames(geno_mat_char)

cat("\nStarting minimal marker set selection ('0' treated as valid allele)...\n")
while (length(selected_markers) < length(remaining_markers)) {
  if (length(unique(current_fingerprints)) == n_var) {
    cat("✅ Unique fingerprints achieved! Total markers used:", length(selected_markers), "\n")
    break
  }
  best_marker <- NULL
  max_unique_count <- -1
  for (marker in remaining_markers) {
    if (marker %in% selected_markers) next
    new_char <- geno_mat_char[, marker]
    temp_fingerprints <- paste0(current_fingerprints, new_char)
    unique_count <- length(unique(temp_fingerprints))
    if (unique_count > max_unique_count) {
      max_unique_count <- unique_count
      best_marker <- marker
    }
  }
  if (is.null(best_marker)) break
  selected_markers <- c(selected_markers, best_marker)
  current_fingerprints <- paste0(current_fingerprints, geno_mat_char[, best_marker])
  cat("Selected", length(selected_markers), "markers; can distinguish", max_unique_count, "accessions\n")
  remaining_markers <- setdiff(remaining_markers, best_marker)
}

fingerprint_df_min <- data.frame(
  Variety = rownames(geno_mat),
  Fingerprint = current_fingerprints,
  stringsAsFactors = FALSE
)
dups_min <- fingerprint_df_min[duplicated(fingerprint_df_min$Fingerprint), "Variety"]
if (length(dups_min) > 0) {
  cat("⚠️ Duplicates remain in minimal set:\n")
  print(dups_min)
} else {
  cat("✅ Minimal marker set uniquely distinguishes all accessions!\n")
}

# Save results of minimal set
write.csv(fingerprint_df_min, "rice_minimal_fingerprint_set.csv", row.names = FALSE)
writeLines(
  paste(fingerprint_df_min$Variety, fingerprint_df_min$Fingerprint, sep = ": "),
  "rice_minimal_fingerprint_set.txt"
)
writeLines(selected_markers, "rice_minimal_marker_list.txt")

# ----------------------------
# End of script
# ----------------------------
cat("\n🎉 All analyses completed! Results (including TIFF and A4 PDF) saved.\n")