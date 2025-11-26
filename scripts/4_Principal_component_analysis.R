# Load required packages
library(ade4)
library(ggplot2)
library(cluster)

# Set working directory (adjust to your actual path)
setwd("/path_to_file")

# Read lower-triangle genetic similarity matrix and reconstruct full symmetric matrix
data <- read.csv("./rice_accessions_genetic_similarity_lower_triangle.csv", row.names = 1, check.names = FALSE)
mat_lower <- as.matrix(data)
mat_lower[is.na(mat_lower)] <- 0  # Replace missing values with 0 temporarily
mat_full <- pmax(mat_lower, t(mat_lower), na.rm = TRUE)  # Make symmetric
diag(mat_full) <- 1  # Set diagonal to 1 (self-similarity)

# Convert to distance matrix and perform Principal Coordinates Analysis (PCoA)
dist_matrix <- as.dist(1 - mat_full)
pcoa_results <- dudi.pco(dist_matrix, scannf = FALSE, nf = 2)
pcoa_df <- data.frame(
  Axis1 = pcoa_results$li[, 1],
  Axis2 = pcoa_results$li[, 2],
  Variety = rownames(pcoa_results$li)
)

# Perform K-means clustering (3 clusters)
set.seed(123)
kmeans_result <- kmeans(pcoa_df[, c("Axis1", "Axis2")], centers = 3)
pcoa_df$Cluster <- factor(kmeans_result$cluster)

# Plot: all points + 75% confidence ellipses around core points of each cluster
pcoa_plot <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(
    geom = "polygon",
    level = 0.75,
    alpha = 0.15,
    aes(fill = Cluster),
    show.legend = FALSE
  ) +
  labs(
    title = "Principal Coordinates Analysis (PCoA)",
    x = "First principal coordinate",
    y = "Second principal coordinate",
    color = "Subpopulation"
  ) +
  scale_color_manual(values = c("red", "forestgreen", "blue")) +
  scale_fill_manual(values = c("red", "forestgreen", "blue")) +
  theme_minimal() +
  coord_fixed() +
  theme(
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),      # Plot title
    axis.title = element_text(size = 30, face = "bold"),                   # Axis labels (x/y)
    axis.text = element_text(size = 30),                                   # Axis tick labels
    legend.title = element_text(size = 30, face = "bold"),                 # Legend title ("Subpopulation")
    legend.text = element_text(size = 30)                                  # Legend labels ("1", "2", "3")
  )

# Save plots
ggsave("Principal component analysis (PCA) of 114 Shanlan upland rice landraces.tiff", plot = pcoa_plot, width = 6.27, height = 8.27, dpi = 300, device = "tiff")
ggsave("Principal component analysis (PCA) of 114 Shanlan upland rice landraces.pdf", plot = pcoa_plot, width = 11.69, height = 8.27, device = "pdf")

cat("✅ All accessions plotted; 75% confidence ellipses drawn around core regions of each cluster.\n")