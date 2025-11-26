setwd("/path_to_file")

# Load required package (install if missing)
if (!require(igraph)) install.packages("igraph")
library(igraph)

# ========== 1. Read lower-triangle SMC matrix ==========
cat("📥 Reading SMC matrix...\n")
s <- read.csv("rice_accessions_genetic_similarity_lower_triangle.csv", row.names = 1, check.names = FALSE, na.strings = "")
mat_lower <- as.matrix(s)

# ========== 2. Reconstruct full symmetric matrix ==========
cat("🔄 Reconstructing symmetric matrix...\n")

# Initialize a full NA matrix
n <- nrow(mat_lower)
mat_full <- matrix(NA, n, n)
rownames(mat_full) <- rownames(mat_lower)
colnames(mat_full) <- colnames(mat_lower)

# Fill lower triangle (including diagonal)
mat_full[lower.tri(mat_full, diag = TRUE)] <- mat_lower[lower.tri(mat_lower, diag = TRUE)]

# Mirror to upper triangle for symmetry
mat_full[upper.tri(mat_full)] <- t(mat_full)[upper.tri(mat_full)]

# Enforce diagonal = 1 (self-similarity in SMC)
diag(mat_full) <- 1

# Check for remaining NAs
if (any(is.na(mat_full))) {
  stop("❌ Matrix reconstruction failed: missing values still present!")
}

cat("✅ SMC matrix successfully reconstructed as a symmetric matrix.\n")

# ========== 3. Set similarity threshold ==========
threshold <- 0.85  # Can be adjusted (e.g., 0.75, 0.85)
cat("🎯 Using similarity threshold:", threshold, "\n")

# ========== 4. Build undirected graph (ensure adjacency matrix is symmetric) ==========
cat("🕸️ Building similarity network...\n")

# Create logical adjacency matrix
adj_raw <- mat_full >= threshold

# Enforce symmetry (guard against floating-point inconsistencies)
adj <- adj_raw | t(adj_raw)

# Remove self-loops
diag(adj) <- FALSE

# Safety check: confirm symmetry
if (!isSymmetric(adj)) {
  stop("❌ Adjacency matrix is not symmetric; cannot build undirected graph!")
}

# Construct graph
g <- graph_from_adjacency_matrix(
  adj,
  mode = "undirected",
  weighted = NULL,
  diag = FALSE
)

cat("✅ Graph built successfully. Nodes:", vcount(g), ", Edges:", ecount(g), "\n")

# ========== 5. Identify connected components (similarity groups) ==========
cat("🧩 Identifying similarity groups...\n")
clusters_list <- clusters(g)
groups <- split(names(clusters_list$membership), clusters_list$membership)

# ========== 6. Select core germplasm (first accession alphabetically per group) ==========
core_accessions <- sapply(groups, function(x) sort(x)[1])

# ========== 7. Identify unique accessions (singletons) ==========
singleton_groups <- sapply(groups, length) == 1
unique_accessions <- core_accessions[singleton_groups]
num_unique <- length(unique_accessions)

# ========== 8. Print summary ==========
cat("\n✅ Number of core accessions:", length(core_accessions), "\n")
cat("✅ Number of unique accessions (no similar partners):", num_unique, "\n")
cat("✅ Core accession list:\n")
print(core_accessions)

# ========== 9. Save results ==========
writeLines(core_accessions, "core_accessions.txt")
writeLines(unique_accessions, "unique_accessions.txt")
write.csv(mat_full, "SMC_matrix_full.csv", quote = FALSE)

# Save detailed grouping info
group_df <- data.frame(
  GroupID = rep(names(groups), lengths(groups)),
  Accession = unlist(groups),
  IsCore = unlist(lapply(groups, function(x) sort(x)[1] == x))
)
write.csv(group_df, "similarity_groups.csv", row.names = FALSE)

cat("\n📁 Saved the following files:\n")
cat("- core_accessions.txt: List of core accessions\n")
cat("- unique_accessions.txt: List of unique (singleton) accessions\n")
cat("- SMC_matrix_full.csv: Full reconstructed SMC matrix\n")
cat("- similarity_groups.csv: Detailed similarity group membership (with core flag)\n")

cat("\n🎉 Analysis completed!\n")

# ========== 10. Visualize network ==========
# Color by connected component
comp <- clusters(g)
colors <- rainbow(max(comp$membership))
V(g)$color <- colors[comp$membership]

# Label only core accessions
core_set <- core_accessions
V(g)$label <- ifelse(V(g)$name %in% core_set, V(g)$name, "")

# Generate compact layout
l <- layout_with_fr(g, start.temp = 0.4, niter = 800)
l <- l * 0.7  # Compress layout

# Save as PDF
pdf("Network visualization of genetic similarity among Shanlan upland rice landraces.pdf", width = 10, height = 8)
plot(g,
     layout = l,
     vertex.size = 10,
     vertex.label.cex = 0.5,
     edge.width = 0.85,
     vertex.color = ifelse(V(g)$name %in% unique_accessions, "gold", "lightblue"),
     main = "Core germplasm of Shanlan upland rice landraces (SMC ≥ 0.8)"
)
dev.off()

# Save as high-resolution TIFF
tiff("Network visualization of genetic similarity among Shanlan upland rice landraces.tiff", width = 10, height = 8, units = "in", res = 300)
plot(g,
     layout = l,
     vertex.size = 10,
     vertex.label.cex = 0.5,
     edge.width = 0.85,
     vertex.color = ifelse(V(g)$name %in% unique_accessions, "gold", "lightblue"),
     main = "Core germplasm of Shanlan upland rice landraces (SMC ≥ 0.8)"
)
dev.off()

# ========== 11. Generate "Core Accession – Similar Accessions" table ==========
cat("📊 Generating core accession to similar accessions mapping...\n")

# Sort each group alphabetically
core_to_group <- lapply(groups, sort)
core_names <- sapply(core_to_group, function(x) x[1])

# Build wide-format data frame
max_group_size <- max(lengths(core_to_group))
group_matrix <- sapply(core_to_group, function(x) {
  c(x, rep("", max_group_size - length(x)))
})

group_df_wide <- as.data.frame(t(group_matrix), stringsAsFactors = FALSE)
colnames(group_df_wide)[1] <- "Core_Accession"
if (ncol(group_df_wide) > 1) {
  colnames(group_df_wide)[-1] <- paste0("Member_", 2:ncol(group_df_wide))
}

# Sort by core accession name
group_df_wide <- group_df_wide[order(group_df_wide$Core_Accession), , drop = FALSE]
row.names(group_df_wide) <- NULL

# Save
write.csv(group_df_wide, "core_accession_similarity_groups.csv", row.names = FALSE, quote = FALSE)
cat("✅ Saved core-to-members mapping: core_accession_similarity_groups.csv\n")

# ========== 12. Count total similar pairs above threshold ==========
upper_mask <- upper.tri(mat_full, diag = FALSE)
num_pairs_above_threshold <- sum(mat_full[upper_mask] >= threshold, na.rm = TRUE)
cat("🔍 Number of non-diagonal accession pairs with SMC ≥", threshold, ":", num_pairs_above_threshold, "\n")