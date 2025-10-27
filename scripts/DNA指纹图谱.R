# 设置工作目录（请根据你的实际路径调整）
setwd("/your_file_path")

# 安装必要包（首次运行时取消注释）
# install.packages(c("adegenet", "poppr", "ape", "ggplot2", "dplyr", "pheatmap"))

# 加载所需包
library(adegenet)
library(poppr)
library(ape)
library(ggplot2)
library(dplyr)
library(pheatmap)

# ----------------------------
# 1. 读取并转置数据
# ----------------------------
raw_tall <- read.csv("DNA指纹图谱.csv", header = TRUE, row.names = 1, check.names = FALSE)
raw_data <- t(raw_tall)  # 117 × 38
cat("✅ 数据已加载，0 视为有效等位基因（无扩增状态）。\n")

# ----------------------------
# 2. 计算遗传相似性系数（SMC）
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
cat("✅ 遗传相似性系数（SMC）计算完成。\n")

# ----------------------------
# 3. 输出下三角遗传相似性系数矩阵（保留2位小数）
# ----------------------------
output_matrix <- matrix("", nrow = nrow(smc_matrix), ncol = ncol(smc_matrix),
                        dimnames = dimnames(smc_matrix))
lower_idx <- lower.tri(smc_matrix, diag = TRUE)
output_matrix[lower_idx] <- sprintf("%.2f", smc_matrix[lower_idx])

write.csv(output_matrix, "水稻品种_遗传相似性系数_下三角.csv", quote = FALSE)
cat("✅ 下三角遗传相似性系数矩阵（保留2位小数）已保存。\n")
cat("\n--- 遗传相似性系数（下三角，前10个品种）---\n")
print(output_matrix[1:10, 1:10], quote = FALSE)

# ----------------------------
# 4. 构建遗传距离矩阵
# ----------------------------
dist_matrix <- as.dist(1 - smc_matrix)

# ----------------------------
# 5. 构建 NJ 树（输出 TIFF + A4 PDF）
# ----------------------------
nj_tree <- nj(dist_matrix)

# --- TIFF ---
# --- TIFF ---
tiff("rice_NJ_dendrogram_phylogram.tiff", width = 2480, height = 3508, res = 300, compression = "lzw")
plot(nj_tree, type = "phylogram", cex = 0.55,
     main = "水稻品种 NJ 聚类树 (38 InDel 标记)",
     font = 1)  # ← 添加这一行
dev.off()

# --- A4 PDF ---
pdf("rice_NJ_dendrogram_A4_phylogram.pdf", width = 8.27, height = 11.69)
plot(nj_tree, type = "phylogram", cex = 0.6,
     main = "水稻品种 NJ 聚类树 (38 InDel 标记)",
     font = 1)  # ← 添加这一行
dev.off()

# ----------------------------
# 6. PCoA 分析（输出 TIFF + A4 PDF）
# ----------------------------
pcoa_results <- dudi.pco(dist_matrix, scannf = FALSE, nf = 2)
pcoa_df <- data.frame(
  Axis1 = pcoa_results$li[, 1],
  Axis2 = pcoa_results$li[, 2],
  Variety = rownames(pcoa_results$li)
)

pcoa_plot <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2)) +
  geom_point(size = 2.5, color = "darkred", alpha = 0.8) +
  labs(
    title = "水稻品种 PCoA 分析 (基于 38 个 InDel 标记)",
    x = paste0("PCoA 轴 1 (", round(pcoa_results$eig[1] / sum(pcoa_results$eig) * 100, 1), "%)"),
    y = paste0("PCoA 轴 2 (", round(pcoa_results$eig[2] / sum(pcoa_results$eig) * 100, 1), "%)")
  ) +
  theme_minimal() +
  coord_fixed()

# --- TIFF ---
ggsave("rice_PCoA.tiff", plot = pcoa_plot, width = 6.27, height = 11.69, dpi = 300, device = "tiff")

# --- A4 PDF ---
ggsave("rice_PCoA_A4.pdf", plot = pcoa_plot, width = 6.27, height = 11.69, device = "pdf")
cat("✅ PCoA 图已保存为 TIFF 和 A4 PDF。\n")

# ----------------------------
# 7. 生成数字指纹并检测重复组
# ----------------------------
raw_data_char <- apply(geno_mat, c(1,2), as.character)
digital_fingerprint <- apply(raw_data_char, 1, function(x) paste(x, collapse = ""))

fingerprint_df <- data.frame(
  Variety_ID = names(digital_fingerprint),
  Digital_Fingerprint = digital_fingerprint,
  stringsAsFactors = FALSE
)

# 按指纹分组，找出重复组（出现次数 >=2）
dup_groups <- fingerprint_df %>%
  group_by(Digital_Fingerprint) %>%
  filter(n() >= 2) %>%
  group_split()

if (length(dup_groups) == 0) {
  cat("✅ 所有 117 个品种均有唯一 DNA 指纹（0 视为有效状态）。\n")
} else {
  cat("⚠️ 发现重复指纹，共", length(dup_groups), "组重复：\n\n")
  
  # 准备输出文本
  output_lines <- c("⚠️ 重复指纹分组结果：\n")
  
  for (i in seq_along(dup_groups)) {
    group_varieties <- dup_groups[[i]]$Variety_ID
    group_str <- paste(group_varieties, collapse = ", ")
    msg <- paste0("组 ", i, ": ", group_str)
    cat(msg, "\n")
    output_lines <- c(output_lines, msg)
  }
  
  # 保存到文件
  writeLines(output_lines, "重复指纹分组.txt")
  cat("\n✅ 重复分组详情已保存为 '重复指纹分组.txt'\n")
}

write.csv(fingerprint_df, "水稻品种_数字指纹.csv", row.names = FALSE, quote = FALSE)
# 保存为 TXT：每行 "品种ID,指纹"
write.table(
  fingerprint_df,
  file = "水稻品种_数字指纹.txt",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  eol = "\n"
)
cat("✅ 数字指纹已保存为 TXT（避免 Excel 科学计数法问题）。\n")
# 安全写入 CSV：确保指纹作为文本
write.table(
  fingerprint_df,
  file = "水稻品种_数字指纹111.csv",
  sep = ",",
  row.names = FALSE,
  quote = TRUE,          # 用双引号包围字段（包括指纹）
  col.names = TRUE,
  eol = "\n"
)

# ----------------------------
# 8. 绘制热图（输出 TIFF + A4 PDF）
# ----------------------------
# 由于热图较宽，A4 横向更合适（11.69 × 8.27）
heat_plot <- function() {
  pheatmap(
    mat = geno_mat,
    color = colorRampPalette(c("gray80", "green", "white", "red", "orange"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    main = "水稻 DNA 指纹热图 (117 品种 × 38 InDel)",
    fontsize = 8
  )
}

# --- TIFF ---
tiff("rice_DNA_fingerprint_heatmap.tiff", width = 3508, height = 2480, res = 300, compression = "lzw")
heat_plot()
dev.off()

# --- A4 PDF（横向）---
pdf("rice_DNA_fingerprint_heatmap_A4.pdf", width = 11.69, height = 8.27)
heat_plot()
dev.off()
cat("✅ 热图已保存为 TIFF 和 A4 PDF（横向）。\n")

# ----------------------------
# 9. 筛选最小标记集
# ----------------------------
geno_mat_char <- apply(geno_mat, c(1,2), as.character)
n_var <- nrow(geno_mat_char)
selected_markers <- character(0)
current_fingerprints <- rep("", n_var)
remaining_markers <- colnames(geno_mat_char)

cat("\n开始筛选最小标记集（0 视为有效等位基因）...\n")
while (length(selected_markers) < length(remaining_markers)) {
  if (length(unique(current_fingerprints)) == n_var) {
    cat("✅ 已找到唯一指纹！共使用", length(selected_markers), "个标记。\n")
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
  cat("已选", length(selected_markers), "个标记，可区分", max_unique_count, "个品种\n")
  remaining_markers <- setdiff(remaining_markers, best_marker)
}

fingerprint_df_min <- data.frame(
  Variety = rownames(geno_mat),
  Fingerprint = current_fingerprints,
  stringsAsFactors = FALSE
)
dups_min <- fingerprint_df_min[duplicated(fingerprint_df_min$Fingerprint), "Variety"]
if (length(dups_min) > 0) {
  cat("⚠️ 最小集仍有重复品种:\n")
  print(dups_min)
} else {
  cat("✅ 最小标记集可唯一区分所有品种！\n")
}
write.csv(fingerprint_df_min, "水稻_最小指纹集_指纹.csv", row.names = FALSE)
writeLines(
  paste(fingerprint_df_min$Variety, fingerprint_df_min$Fingerprint, sep = ": "),
  "水稻_最小指纹集_指纹.txt"
)
writeLines(selected_markers, "水稻_最小指纹集_标记列表.txt")

# ----------------------------
# 结束
# ----------------------------
cat("\n🎉 全部分析完成！所有结果（含 TIFF + A4 PDF 图）已保存。\n")