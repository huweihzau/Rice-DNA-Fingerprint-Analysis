# 1. Load required libraries
# ------------------------------------------------------------------------------
library(ggplot2)
library(ape)
library(ggtree)
library(dplyr)
# If you need the encircle effect, ggtree may depend on ggalt or ggforce.
# Usually, ggtree handles this automatically.
# If you encounter errors, try: install.packages("ggalt") or install.packages("ggforce")

# 2. Load data (using your full similarity matrix)
# ------------------------------------------------------------------------------
sim_file <- "SMC_matrix_full.csv"
if (!file.exists(sim_file)) stop("Error: File 'SMC_matrix_full.csv' not found!")

# Read the data
sim_data <- read.csv(sim_file, row.names = 1, check.names = FALSE)
sim_matrix <- as.matrix(sim_data)

# 3. Transform similarity to genetic distance
# ------------------------------------------------------------------------------
# Distance = 1 - Similarity
dist_matrix <- as.dist(1 - sim_matrix)

# 4. Build UPGMA clustering tree
# ------------------------------------------------------------------------------
upgma_tree <- hclust(dist_matrix, method = "average")
phylo_tree <- as.phylo(upgma_tree)

# 5. Group samples (automatically split into 3 clusters)
# ------------------------------------------------------------------------------
groups <- cutree(upgma_tree, k = 3)

# Helper function: find the MRCA node for each group
get_mrca_node <- function(tree, group_vector, target_group) {
  tips <- names(group_vector[group_vector == target_group])
  if (length(tips) == 1) {
    return(which(tree$tip.label == tips)) 
  } else {
    return(MRCA(tree, tips)) 
  }
}

node_g1 <- get_mrca_node(phylo_tree, groups, 1)
node_g2 <- get_mrca_node(phylo_tree, groups, 2)
node_g3 <- get_mrca_node(phylo_tree, groups, 3)

# 6. Visualization (circular layout with encircling highlights)
# ------------------------------------------------------------------------------
# Key: type = "encircle" draws a circular outline around clades,
# instead of filling a wedge from the center.

# Create base tree plot
p <- ggtree(phylo_tree, layout = "circular", open.angle = 10, size = 0.5)

# Add encircling highlights
# alpha controls transparency; expand controls distance from branches
p <- p + 
  geom_hilight(node = node_g1, fill = "#F8766D", alpha = 0.4, type = "encircle", expand = 0.02) + 
  geom_hilight(node = node_g2, fill = "#00BA38", alpha = 0.4, type = "encircle", expand = 0.02) + 
  geom_hilight(node = node_g3, fill = "#619CFF", alpha = 0.4, type = "encircle", expand = 0.02)

# Add tip labels
p <- p + geom_tiplab(size = 2, offset = 0.01, align = TRUE, linesize = 0.2) 

# Add title and theme settings
tree_plot <- p +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none" 
  ) +
  labs(title = "UPGMA Clustering (Encircle Style)")

# Display the plot
print(tree_plot)

# 7. Save output
# ------------------------------------------------------------------------------
if (!dir.exists("output")) dir.create("output")
ggsave("output/Phylogenetic_tree_Encircle.pdf", plot = tree_plot, width = 10, height = 10)
ggsave("output/Phylogenetic_tree_Encircle.png", plot = tree_plot, width = 10, height = 10, dpi = 300)

message("Plotting complete! Saved as encircle-style phylogenetic tree.")