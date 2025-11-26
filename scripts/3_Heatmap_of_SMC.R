# Set working directory (adjust to your actual path)
setwd("/path_to_file")

# Load required package
library(pheatmap)

# Read the lower-triangle genetic similarity matrix
data <- read.csv("./rice_accessions_genetic_similarity_lower_triangle.csv", row.names = 1, check.names = FALSE)

# Convert to matrix and complete it into a symmetric matrix
mat_lower <- as.matrix(data)
mat_lower[is.na(mat_lower)] <- 0  # Temporarily replace missing values with 0

# Construct full symmetric similarity matrix
mat_full <- pmax(mat_lower, t(mat_lower), na.rm = TRUE)

# Ensure diagonal elements are 1 (similarity of an accession with itself)
diag(mat_full) <- 1

# --- Export as high-resolution TIFF ---
tiff("Simple matching coefficient heatmap of Shanlan upland rice landraces.tiff", width = 3508, height = 2480, res = 300, compression = "lzw")

# Plot heatmap and save as TIFF
pheatmap(
  mat_full,
  color = colorRampPalette(c("white", "orange", "red"))(300),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,   # Set to FALSE to avoid label overlap when many accessions
  show_colnames = FALSE,
  border_color = NA,
  main = "Simple matching coefficient heatmap of Shanlan upland rice landraces"
)

# Close TIFF graphics device
dev.off()

# --- Export as A4 PDF (landscape orientation) ---
pdf("Simple matching coefficient heatmap of Shanlan upland rice landraces.pdf", width = 11.69, height = 8.27)  # A4 landscape size in inches

# Plot heatmap and save as PDF
pheatmap(
  mat_full,
  color = colorRampPalette(c("white", "orange", "red"))(300),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,   # Set to FALSE to avoid label overlap when many accessions
  show_colnames = FALSE,
  border_color = NA,
  main = "Simple matching coefficient heatmap of Shanlan upland rice landraces"
)

# Close PDF graphics device
dev.off()

cat("✅ Heatmap saved as TIFF and A4 PDF (landscape).\n")