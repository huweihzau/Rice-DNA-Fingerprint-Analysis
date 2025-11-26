# If ggplot2 and ggcorrplot packages are not installed, install them first
# install.packages("ggplot2")
# install.packages("ggcorrplot")

# Clear all objects from the workspace
rm(list = ls())

# Load required libraries
library(corrplot)
library(ggplot2)
library(ggcorrplot)
library(vcd)
library(psych)
library(ggrepel)
library(readxl)
library(dplyr)
library(reshape2)
library(gridExtra)
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(showtext)

# Enable showtext for custom font support
showtext_auto()

# Add Times New Roman font (adjust path if needed for your system)
font_add("Times New Roman", regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")

# Re-enable showtext after adding the font
showtext_auto()

# Set default plotting font to Times New Roman
par(family = "Times New Roman")

# Set working directory
setwd("/path_to_file")

# Read phenotype data from CSV file
data <- read.table("phenotype.csv", header = TRUE, sep = ",")  # Remove row.names=1 if it causes errors

# Function to compute correlation matrix and corresponding p-values matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate Pearson correlation matrix (using complete observations only)
cor_matrix <- cor(data, use = "complete.obs")

# Calculate p-values matrix for correlations
p_values <- cor.mtest(data)

# Set up PDF output device
pdf("./Correlation analysis of major traits in Shanlan upland rice landraces.pdf", width = 8.27, height = 11.69)  # A4 size in inches

# Define a smooth diverging color palette: dark blue → white → orange
color_palette <- colorRampPalette(c("#0E3461", "#FFFFFF", "#D87C4C"))(300)

# Plot the full correlation matrix using circles
corrplot(cor_matrix, 
         method = "circle", 
         type = "full", 
         outline = FALSE, 
         diag = TRUE, 
         mar = c(0, 0, 1, 0), 
         bg = "white", 
         add = FALSE, 
         is.corr = TRUE, 
         addgrid.col = "darkgray", 
         order = "original", 
         tl.pos = NULL, 
         tl.cex = 1, 
         tl.col = "black", 
         tl.srt = 45,        # Rotate axis labels by 45 degrees
         tl.font = 2,        # Bold axis labels
         cl.pos = "r", 
         col = color_palette, 
         cl.lim = c(-1, 1))  # Correlation range from -1 to 1

# Overlay numerical correlation coefficients in the lower triangle
num_color_palette <- colorRampPalette(c("#0E3461", "#D87C4C"))(300)
corrplot(cor_matrix, 
         method = "number", 
         type = "lower", 
         order = "original", 
         col = num_color_palette, 
         diag = FALSE, 
         tl.pos = "n", 
         cl.pos = NULL, 
         add = TRUE)

# Add significance stars: * for p < 0.05, ** for p < 0.01
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i < j) {
      if (p_values[i, j] < 0.01) {
        text(j, nrow(cor_matrix) - i + 1, "**", col = "black", cex = 1.5)
      } else if (p_values[i, j] < 0.05 && p_values[i, j] >= 0.01) {
        text(j, nrow(cor_matrix) - i + 1, "*", col = "black", cex = 1.5)
      }
    }
  }
}

# Close the PDF device to save the plot
dev.off()

# Save correlation matrix to CSV
write.csv(cor_matrix, file = "correlation_matrix.csv")

# Save p-values matrix to CSV
write.csv(p_values, file = "p_values_matrix.csv")