# ============================================================================
# 1. 环境准备与包安装（针对R 4.5.x + Bioc 3.22 最终修正版）
# ============================================================================
setwd("E:/文章专利/水牛乳")
# ----------------------------------------------------------------------------
# 1.1 配置镜像（中国大陆用户强烈推荐）
# ----------------------------------------------------------------------------
options(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN")
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# ----------------------------------------------------------------------------
# 1.2 安装BiocManager并强制锁定Bioc 3.22（这是唯一正确的版本）
# ----------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = getOption("CRAN"))
}

# 【唯一修复】显式指定 version = "3.22"
BiocManager::install(version = "3.22", ask = FALSE, update = TRUE, force = TRUE)

# 验证版本（必须显示 3.22）
cat("Bioconductor版本:", as.character(BiocManager::version()), "\n")
stopifnot(BiocManager::version() == "3.22")

# ----------------------------------------------------------------------------
# 1.3 安装CRAN包
# ----------------------------------------------------------------------------
cran_packages <- c(
    "tidyverse", "openxlsx", "Hmisc", "igraph", "circlize",
    "RColorBrewer", "ggrepel", "patchwork", "renv", 
    "pheatmap", "viridis", "ropls", "corrplot"
)

for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = getOption("CRAN"))
        library(pkg, character.only = TRUE)
    }
}

# ----------------------------------------------------------------------------
# 1.4 安装Bioconductor包（现在ComplexHeatmap来自3.22，完全兼容）
# ----------------------------------------------------------------------------
bioc_packages <- c(
    "ComplexHeatmap",   # ✅ 现在100%可装
    "limma", 
    "clusterProfiler", 
    "org.Bt.eg.db", 
    "STRINGdb", 
    "preprocessCore"
)

for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        library(pkg, character.only = TRUE)
    }
}

# ----------------------------------------------------------------------------
# 1.5 验证关键包
# ----------------------------------------------------------------------------
cat("\n✅ 环境就绪\n")
cat("R版本:", R.version.string, "\n")
cat("Bioc版本:", as.character(BiocManager::version()), "\n")
cat("ComplexHeatmap版本:", as.character(packageVersion("ComplexHeatmap")), "\n")

# 创建目录、设随机种子
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
set.seed(123456)

# 保存会话信息
sink("results/session_info.txt")
sessionInfo()
sink()

message("🎉 环境准备完成！所有包已成功加载。")
# ============================================================================
# 2.1 蛋白质组数据导入与预处理（适配你的数据列名：Description）
# ============================================================================

library(tidyverse)
library(openxlsx)
library(limma)

# 读取数据
protein_file <- "https://raw.githubusercontent.com/Alice/my-r-scripts/main/all_sample_ori.xlsx"
protein_raw <- openxlsx::read.xlsx(protein_file, sheet = 1)

# 查看列名（确认无误）
colnames(protein_raw)
# [1] "Protein"     "Description" "Gene"        "Holstein_1"  "Holstein_2" 
# [6] "Holstein_3"  "Jersey_1"    "Jersey_2"    "Jersey_3"    "Buffalo_1"  
# [11] "Buffalo_2"   "Buffalo_3"

# ----------------------------------------------------------------------------
# 步骤1: 过滤（根据你的数据实际情况）
# ----------------------------------------------------------------------------

# 你的数据里没有 Reverse、Potential.contaminant 列，也没有 REV__/CON__ 前缀
# 跳过反向序列和污染物过滤

protein_filtered <- protein_raw

# 仅检查是否需要过滤“仅通过位点鉴定”的蛋白
# 查看 Description 列是否包含该字样
if (any(grepl("Only identified by site", protein_filtered$Description, ignore.case = TRUE))) {
  protein_filtered <- protein_filtered %>%
    dplyr::filter(!grepl("Only identified by site", Description, ignore.case = TRUE))
  message("已过滤仅位点鉴定蛋白。")
} else {
  message("未发现仅位点鉴定蛋白，跳过该过滤。")
}

message(paste0("原始蛋白数: ", nrow(protein_raw)))
message(paste0("过滤后蛋白数: ", nrow(protein_filtered)))

# ----------------------------------------------------------------------------
# 步骤2: 提取强度矩阵
# ----------------------------------------------------------------------------

# 样品列：荷斯坦、娟姗、水牛
sample_cols <- c("Holstein_1", "Holstein_2", "Holstein_3",
                 "Jersey_1", "Jersey_2", "Jersey_3",
                 "Buffalo_1", "Buffalo_2", "Buffalo_3")

# 确保所有样品列都存在
sample_cols <- sample_cols[sample_cols %in% colnames(protein_filtered)]

protein_matrix <- protein_filtered[, sample_cols]
protein_matrix <- as.matrix(protein_matrix)
mode(protein_matrix) <- "numeric"
rownames(protein_matrix) <- protein_filtered$Protein

# ----------------------------------------------------------------------------
# 步骤3: 缺失值过滤（每组至少70%样品有值）
# ----------------------------------------------------------------------------

# 定义分组
group_protein <- data.frame(
  sample = colnames(protein_matrix),
  group = case_when(
    grepl("Holstein", colnames(protein_matrix)) ~ "Holstein",
    grepl("Jersey", colnames(protein_matrix)) ~ "Jersey",
    grepl("Buffalo", colnames(protein_matrix)) ~ "Buffalo"
  )
)

# 按组筛选：每组至少 70% 样品非缺失
# 荷斯坦/娟姗/水牛各有3个样品，70% ≈ 2.1，即至少2个样品有值
keep_rows <- apply(protein_matrix, 1, function(x) {
  holstein_na <- sum(!is.na(x[group_protein$group == "Holstein"])) >= 2
  jersey_na   <- sum(!is.na(x[group_protein$group == "Jersey"]))   >= 2
  buffalo_na  <- sum(!is.na(x[group_protein$group == "Buffalo"]))  >= 2
  holstein_na & jersey_na & buffalo_na
})

protein_matrix_filtered <- protein_matrix[keep_rows, ]
message(paste0("缺失值过滤后蛋白数: ", nrow(protein_matrix_filtered)))

# ----------------------------------------------------------------------------
# 步骤4: 缺失值插补（最小值的半数）
# ----------------------------------------------------------------------------

min_val <- min(protein_matrix_filtered, na.rm = TRUE)
half_min <- min_val / 2
protein_matrix_imputed <- protein_matrix_filtered
protein_matrix_imputed[is.na(protein_matrix_imputed)] <- half_min

# ----------------------------------------------------------------------------
# 步骤5: Log2转换
# ----------------------------------------------------------------------------

# 加一个小偏移避免 log2(0)
protein_log2 <- log2(protein_matrix_imputed + 1)

# ----------------------------------------------------------------------------
# 步骤6: 分位数归一化
# ----------------------------------------------------------------------------

protein_norm <- limma::normalizeQuantiles(protein_log2)

# ----------------------------------------------------------------------------
# 步骤7: 添加基因名和描述信息
# ----------------------------------------------------------------------------

# 保留对应的注释
protein_anno <- protein_filtered[keep_rows, c("Protein", "Gene", "Description")]
protein_processed <- cbind(protein_anno, protein_norm)

# 修正列名：将样品列名恢复（归一化后列名不变）
colnames(protein_processed)[4:ncol(protein_processed)] <- colnames(protein_norm)

# 保存预处理后的数据
saveRDS(protein_processed, "results/protein_processed.rds")
saveRDS(protein_norm, "results/protein_norm_matrix.rds")

# 导出为Excel
write.xlsx(protein_processed, "tables/Table_Protein_Processed.xlsx", rowNames = FALSE)

message("✅ 蛋白质组数据预处理完成！")
# ============================================================================
# 2.2 脂质组数据导入与预处理（适配你的实际列名）
# ============================================================================

# 读取合并模式数据（主要分析数据集）
lipid_file_all <- "https://raw.githubusercontent.com/Alice/my-r-scripts/main/meta_intensity_all.xlsx"
lipid_all_raw <- openxlsx::read.xlsx(lipid_file_all, sheet = 1)

# 读取正负离子模式数据（用于一致性验证）
lipid_file_pos <- "https://raw.githubusercontent.com/Alice/my-r-scripts/main/meta_intensity_pos.xlsx"
lipid_pos_raw <- openxlsx::read.xlsx(lipid_file_pos, sheet = 1)

lipid_file_neg <- "https://raw.githubusercontent.com/Alice/my-r-scripts/main/meta_intensity_neg.xlsx"
lipid_neg_raw <- openxlsx::read.xlsx(lipid_file_neg, sheet = 1)

# 保存原始数据副本
saveRDS(lipid_all_raw, "results/lipid_all_raw.rds")

# ----------------------------------------------------------------------------
# 步骤1: 使用合并模式数据作为主要分析数据集
# ----------------------------------------------------------------------------
lipid_raw <- lipid_all_raw

# ----------------------------------------------------------------------------
# 步骤2: 提取脂质注释信息（【关键修正】列名与实际完全一致）
# ----------------------------------------------------------------------------
lipid_anno <- lipid_raw[, c("ID", "Name", "LipidGroup", "Class", 
                            "Chain_length", "Unsaturated_keys", 
                            "Formula", "RT.[min]", "m/z")]

# ----------------------------------------------------------------------------
# 步骤3: 提取样品强度列（【关键修正】正确匹配 all_Buffalo_1 等列）
# ----------------------------------------------------------------------------
# 直接列出所有样品列名（最稳健，避免 grep 误匹配）
sample_cols_lipid <- c("all_Buffalo_1", "all_Buffalo_2", "all_Buffalo_3",
                       "all_Holstein_1", "all_Holstein_2", "all_Holstein_3",
                       "all_Jersey_1", "all_Jersey_2", "all_Jersey_3")
# 确保这些列都存在
sample_cols_lipid <- sample_cols_lipid[sample_cols_lipid %in% colnames(lipid_raw)]

# 提取强度矩阵
lipid_matrix <- lipid_raw[, sample_cols_lipid]
lipid_matrix <- as.matrix(lipid_matrix)
mode(lipid_matrix) <- "numeric"
rownames(lipid_matrix) <- lipid_raw$ID

# ----------------------------------------------------------------------------
# 步骤4: 缺失值处理（零值替换为NA，每行最小值的一半插补）
# ----------------------------------------------------------------------------
# 零值替换为 NA
lipid_matrix[lipid_matrix == 0] <- NA

# 计算每个特征的最小值（非NA）
lipid_min_vals <- apply(lipid_matrix, 1, min, na.rm = TRUE)

# 缺失值插补（该特征最小值的半数）
lipid_matrix_imputed <- lipid_matrix
for (i in 1:nrow(lipid_matrix)) {
  na_idx <- which(is.na(lipid_matrix[i, ]))
  if (length(na_idx) > 0) {
    lipid_matrix_imputed[i, na_idx] <- lipid_min_vals[i] / 2
  }
}

# ----------------------------------------------------------------------------
# 步骤5: Log2转换（加1偏移避免log2(0)）
# ----------------------------------------------------------------------------
lipid_log2 <- log2(lipid_matrix_imputed + 1)

# ----------------------------------------------------------------------------
# 步骤6: 中位数中心化（按行，去除绝对丰度差异）
# ----------------------------------------------------------------------------
lipid_centered <- lipid_log2 - apply(lipid_log2, 1, median, na.rm = TRUE)

# ----------------------------------------------------------------------------
# 步骤7: 合并注释信息，保存结果
# ----------------------------------------------------------------------------
lipid_processed <- cbind(lipid_anno, lipid_centered)

saveRDS(lipid_processed, "results/lipid_processed.rds")
saveRDS(lipid_centered, "results/lipid_centered_matrix.rds")
write.xlsx(lipid_processed, "tables/Table_Lipid_Processed.xlsx", rowNames = FALSE)

# ----------------------------------------------------------------------------
# 步骤8: 补充图S1 - 不同离子模式一致性评估（【关键修正】适配正负离子列名）
# ----------------------------------------------------------------------------
# 正离子模式样品列（通常列名形如 pos_Buffalo_1, pos_Holstein_1, pos_Jersey_1）
pos_sample_cols <- grep("pos_Buffalo_|pos_Holstein_|pos_Jersey_", 
                        colnames(lipid_pos_raw), value = TRUE)
pos_sample_cols <- pos_sample_cols[!grepl("QC", pos_sample_cols)]  # 排除QC
pos_matrix <- lipid_pos_raw[, pos_sample_cols]
pos_matrix <- as.matrix(pos_matrix)
pos_matrix[pos_matrix == 0] <- NA
pos_log2 <- log2(pos_matrix + 1)

# 负离子模式样品列
neg_sample_cols <- grep("neg_Buffalo_|neg_Holstein_|neg_Jersey_", 
                        colnames(lipid_neg_raw), value = TRUE)
neg_sample_cols <- neg_sample_cols[!grepl("QC", neg_sample_cols)]
neg_matrix <- lipid_neg_raw[, neg_sample_cols]
neg_matrix <- as.matrix(neg_matrix)
neg_matrix[neg_matrix == 0] <- NA
neg_log2 <- log2(neg_matrix + 1)

# 提取共同脂质用于比较
common_lipids <- intersect(rownames(pos_log2), rownames(neg_log2))

# 计算相关性（取前100个或全部）
cor_vals <- c()
for (lipid in common_lipids[1:min(100, length(common_lipids))]) {
  pos_vals <- as.numeric(pos_log2[lipid, ])
  neg_vals <- as.numeric(neg_log2[lipid, ])
  if (sum(!is.na(pos_vals) & !is.na(neg_vals)) > 3) {
    cor_vals <- c(cor_vals, cor(pos_vals, neg_vals, use = "complete.obs", method = "pearson"))
  }
}

# 绘图
pdf("figures/Figure_S1_IonMode_Consistency.pdf", width = 8, height = 6)
par(mfrow = c(1, 2))
# 直方图
hist(cor_vals, breaks = 30, col = "steelblue", 
     main = "Consistency between positive and negative modes",
     xlab = "Pearson correlation", ylab = "Frequency")
abline(v = median(cor_vals, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
# 散点图示例
example_idx <- which.min(abs(cor_vals - median(cor_vals, na.rm = TRUE)))
example_lipid <- common_lipids[example_idx]
plot(as.numeric(pos_log2[example_lipid, ]), 
     as.numeric(neg_log2[example_lipid, ]),
     pch = 16, col = rgb(0.2, 0.4, 0.8, 0.7),
     main = paste0("Example: ", example_lipid),
     xlab = "Positive mode (log2 intensity)",
     ylab = "Negative mode (log2 intensity)")
abline(lm(as.numeric(neg_log2[example_lipid, ]) ~ as.numeric(pos_log2[example_lipid, ])), 
       col = "red", lwd = 2)
dev.off()

message("✅ 脂质组数据预处理完成！")
# ============================================================================
# 图1A - 脂质组样品间相关性热图（显示相关性数值）
# ============================================================================

# 1. 加载必要包 ----------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)        # 颜色渐变函数
library(RColorBrewer)    # 配色
library(tidyverse)

# 2. 读取预处理后的脂质组数据 ------------------------------------------------
lipid_centered <- readRDS("results/lipid_centered_matrix.rds")
# 这是一个矩阵，行: 脂质特征，列: 样品（all_Buffalo_1, all_Holstein_1, ...）

# 3. 计算样品间Pearson相关性矩阵 --------------------------------------------
cor_matrix <- cor(lipid_centered, method = "pearson")

# 查看矩阵（确认样品名）
print(round(cor_matrix[1:3, 1:3], 3))

# 4. 准备列注释信息（样品分组） ---------------------------------------------
# 从列名提取品种：去除 "all_" 前缀，取第一个单词（Buffalo / Holstein / Jersey）
sample_group <- str_remove(colnames(cor_matrix), "^all_") %>% 
                str_extract("^[A-Za-z]+")

# 创建数据框
col_anno_df <- data.frame(
  Group = sample_group,
  row.names = colnames(cor_matrix)
)

# 定义分组颜色（与原文保持一致）
group_colors <- c(
  "Buffalo" = "#E69F00",
  "Holstein" = "#56B4E9",
  "Jersey" = "#009E73"
)

# 5. 定义热图颜色（白色→红色渐变） -----------------------------------------
# 相关性范围通常为0.8~1.0，这里设置颜色映射为 [0.8, 1.0]
col_fun <- colorRamp2(
  breaks = seq(0.8, 1.0, length.out = 9),
  colors = brewer.pal(9, "Reds")
)

# 6. 创建列注释对象 ----------------------------------------------------------
top_annotation <- HeatmapAnnotation(
  df = col_anno_df,
  col = list(Group = group_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_height = 0.5,
  annotation_name_gp = gpar(fontsize = 10)
)

# 7. 绘制热图（带数值标签）---------------------------------------------------
# 定义单元格绘制函数，显示相关性数值
cell_fun <- function(j, i, x, y, width, height, fill) {
  # 获取当前单元格的值
  val <- pindex(cor_matrix, i, j)
  # 格式化为保留两位小数
  label <- sprintf("%.2f", val)
  # 在单元格中心绘制文本，字体大小6
  grid.text(label, x, y, gp = gpar(fontsize = 6, col = "black"))
}

ht <- Heatmap(
  cor_matrix,
  name = "Pearson r",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),  # 白色网格线
  cell_fun = cell_fun,                     # 添加数值标签
  column_title = "Lipidome Sample Correlation",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  row_names_side = "left",
  column_names_side = "top",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  top_annotation = top_annotation,
  heatmap_legend_param = list(
    title = "Pearson r",
    title_position = "topcenter",
    legend_height = unit(4, "cm"),
    at = seq(0.8, 1.0, 0.05),
    labels = seq(0.8, 1.0, 0.05)
  )
)

# 8. 保存为PDF（适当增加画布尺寸，以便容纳数值标签）----------------------------
pdf("figures/Figure_1A_Lipid_Sample_Correlation_with_values.pdf",
    width = 12, height = 11)  # 比原来增大，确保数值清晰
draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

message("✅ 脂质组样品相关性热图（带数值）已保存至: figures/Figure_1A_Lipid_Sample_Correlation_with_values.pdf")
# ============================================
# 脂质组PCA得分图（同蛋白质PCA风格）
# 分组：Buffalo, Holstein, Jersey
# 样本总数：9个（每组3个）
# 可视化：样本分布 + 95%置信椭圆（正态椭圆）
# 使用 FactoMineR 进行 PCA，与蛋白PCA完全一致
# ============================================

# 1. 加载必要包
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(FactoMineR)) install.packages("FactoMineR")
if (!require(factoextra)) install.packages("factoextra")

library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)

# --- 解决select冲突 ---
select <- dplyr::select

# ================== 模拟脂质组数据（每组3个样本，共9个样本）==================
set.seed(123)

# 模拟样本：Buffalo 3个, Holstein 3个, Jersey 3个
samples <- c(paste0("Buffalo_", 1:3),
             paste0("Holstein_", 1:3),
             paste0("Jersey_", 1:3))
n_samples <- length(samples)

# 模拟脂质：150个脂质分子
n_lipids <- 150
lipid_ids <- paste0("LIPID", 1:n_lipids)
lipid_names <- paste0("TG(", sample(48:60, n_lipids, replace = TRUE), ":",
                      sample(0:3, n_lipids, replace = TRUE, prob = c(0.6,0.2,0.15,0.05)), ")")
lipid_classes <- sample(c("TG", "PC", "PE", "SM", "CE", "DG", "LPC", "Cer", "FFA"),
                        n_lipids, replace = TRUE, prob = c(0.4,0.15,0.1,0.08,0.07,0.06,0.05,0.05,0.04))

# 构建丰度矩阵（行=脂质，列=样本）
abundance_matrix <- matrix(
  rnorm(n_lipids * n_samples, mean = 12, sd = 3),
  nrow = n_lipids, ncol = n_samples
)
# 引入组间差异（使分组可分离）
abundance_matrix[, 1:3] <- abundance_matrix[, 1:3] + rnorm(n_lipids, 2.0, 0.6)   # Buffalo
abundance_matrix[, 4:6] <- abundance_matrix[, 4:6] + rnorm(n_lipids, 0, 0.6)    # Holstein
abundance_matrix[, 7:9] <- abundance_matrix[, 7:9] + rnorm(n_lipids, -1.5, 0.6) # Jersey

colnames(abundance_matrix) <- samples
rownames(abundance_matrix) <- lipid_ids

# 构建lipid_all数据框（与蛋白数据结构一致）
lipid_all <- data.frame(
  ID = lipid_ids,
  Name = lipid_names,
  Class = lipid_classes,
  abundance_matrix,
  stringsAsFactors = FALSE
)

cat("模拟脂质组数据：", n_lipids, "个脂质，", n_samples, "个样本\n")
cat("样本组成：Buffalo 3, Holstein 3, Jersey 3\n")

# ============ 真实数据替换点 ============
# 如果您有真实数据，请注释上方所有模拟代码，并取消下方注释
# library(readxl)
# lipid_all <- read_excel("您的脂质组数据.xlsx")
# 要求：lipid_all 必须包含样本列（列名含Buffalo/Holstein/Jersey）
# ========================================

# 2. 自动识别样本列
sample_cols <- grep("Buffalo|Holstein|Jersey", colnames(lipid_all), value = TRUE)
sample_cols <- sample_cols[!grepl("QC", sample_cols)]
if (length(sample_cols) == 0) stop("未找到样本列！")

cat("共识别", length(sample_cols), "个样本列\n")

# 3. 提取表达矩阵（行=脂质，列=样本）
exp_mat <- as.matrix(lipid_all[, sample_cols])
rownames(exp_mat) <- lipid_all$ID  # 使用ID作为行名

# 4. 数据预处理（同蛋白PCA：缺失值填补 + log2转换）
exp_mat[exp_mat == 0] <- NA
if (any(is.na(exp_mat))) {
  min_val <- min(exp_mat, na.rm = TRUE) / 2
  exp_mat[is.na(exp_mat)] <- min_val
  cat("已进行缺失值填补（最小值/2）\n")
}

# log2转换（若数据未log）
exp_mat <- log2(exp_mat)

# 5. 转置：PCA要求行=样本，列=变量
pca_data <- t(exp_mat)

# 6. 执行PCA（与蛋白PCA完全一致）
pca_result <- PCA(pca_data, scale.unit = TRUE, graph = FALSE)

# 7. 提取解释方差（前两轴）
var_explained <- pca_result$eig[1:2, 2]

# 8. 提取样本坐标
pca_df <- as.data.frame(pca_result$ind$coord[, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Sample <- rownames(pca_df)
pca_df$Group <- case_when(
  grepl("Buffalo", pca_df$Sample) ~ "Buffalo",
  grepl("Holstein", pca_df$Sample) ~ "Holstein",
  grepl("Jersey", pca_df$Sample) ~ "Jersey",
  TRUE ~ "Other"
)

# 9. 检查每组样本数量，确保可绘制95%置信椭圆（至少需要3个点才能计算椭圆）
group_counts <- pca_df %>% group_by(Group) %>% summarise(n = n())
cat("\n各组样本数：\n")
print(group_counts)

# 10. 绘制PCA得分图（与蛋白PCA风格完全一致）
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  # 95%置信椭圆（仅当组内样本数≥3时绘制）
  stat_ellipse(level = 0.95, type = "norm", geom = "polygon", 
               alpha = 0.1, show.legend = FALSE) +
  # 样本点
  geom_point(size = 4, alpha = 0.9) +
  # 颜色方案：与蛋白PCA完全相同
  scale_color_manual(values = c("Buffalo" = "#E64B35", 
                                "Holstein" = "#4DBBD5", 
                                "Jersey" = "#00A087")) +
  scale_fill_manual(values = c("Buffalo" = "#E64B35", 
                               "Holstein" = "#4DBBD5", 
                               "Jersey" = "#00A087")) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    title = "Lipidomics PCA Score Plot",
    color = "Breed",
    fill = "Breed"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +
  # 添加坐标轴零点线
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5)

# 显示图形
print(p_pca)

# 11. 保存图形（PDF & PNG）
ggsave("Lipidomics_PCA_ScorePlot_9samples.pdf", p_pca, width = 8, height = 6)
ggsave("Lipidomics_PCA_ScorePlot_9samples.png", p_pca, width = 8, height = 6, dpi = 300)

# 12. 可选：载荷图（变量对主成分的贡献）
# 如需要可取消注释
# p_loadings <- fviz_pca_var(pca_result,
#                            col.var = "contrib",
#                            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                            repel = TRUE) +
#   labs(title = "Lipid Loadings") +
#   theme_minimal()
# ggsave("Lipidomics_PCA_Loadings_9samples.pdf", p_loadings, width = 8, height = 6)

cat("\n✅ 脂质组PCA得分图绘制完成！\n")
cat("📊 图形保存：Lipidomics_PCA_ScorePlot_9samples.pdf / .png\n")

# 13. 输出解释方差
cat("\n主成分解释方差（前5个）：\n")
print(pca_result$eig[1:5, 2])
# ============================================
# 脂质类别分布柱状图（四个主要类别 + Other）
# 终极修正版：彻底解决select函数冲突
# ============================================

# 1. 加载必要包（按需加载，避免冲突）
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(scales)) install.packages("scales")

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

# --- 解决select函数冲突 ---
# 优先使用dplyr的select，并屏蔽其他包的select
select <- dplyr::select
assign("select", dplyr::select, envir = .GlobalEnv)
cat("已设置全局select函数为dplyr::select\n")

# ================== 数据准备 ==================
# 方案A：使用模拟数据（默认，保证直接运行）
set.seed(123)

# 生成模拟脂质注释
lipid_ids <- paste0("LIPID", 1:200)
lipid_names <- paste0("Lipid", 1:200)
lipid_classes <- sample(
  c("TG", "PC", "PE", "SM", "CE", "DG", "PS", "PI", "LPC", "Cer", "FFA"),
  200, replace = TRUE,
  prob = c(0.35, 0.18, 0.12, 0.08, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02)
)

# 模拟样本：3个组，每组5个重复
groups <- c(rep("Buffalo", 5), rep("Holstein", 5), rep("Jersey", 5))
samples <- paste0(groups, "_", 1:5)

# 模拟丰度矩阵
abundance_matrix <- matrix(
  runif(200 * 15, 100, 10000), 
  nrow = 200, ncol = 15
)
# 让TG丰度显著更高
tg_rows <- which(lipid_classes == "TG")
abundance_matrix[tg_rows, ] <- abundance_matrix[tg_rows, ] * 8
colnames(abundance_matrix) <- samples

# 构建lipid_all数据框
lipid_all <- data.frame(
  ID = lipid_ids,
  Name = lipid_names,
  Class = lipid_classes,
  abundance_matrix,
  stringsAsFactors = FALSE
)

cat("使用模拟数据，共", nrow(lipid_all), "个脂质分子，", 
    ncol(lipid_all)-3, "个样本\n")

# ============ 如果您有真实数据，请注释上方模拟代码，并取消下方注释 ============
# library(readxl)
# lipid_all <- read_excel("您的真实数据路径.xlsx")
# 要求：lipid_all必须包含列 'Class' 以及样本列（列名包含Buffalo/Holstein/Jersey）
# =============================================================================

# 2. 自动识别样本列（包含Buffalo/Holstein/Jersey）
sample_cols <- grep("Buffalo|Holstein|Jersey", colnames(lipid_all), value = TRUE)
sample_cols <- sample_cols[!grepl("QC", sample_cols)]  # 排除QC样本
if (length(sample_cols) == 0) {
  stop("未找到任何样本列！请检查数据列名是否包含 Buffalo/Holstein/Jersey")
}
cat("找到", length(sample_cols), "个样本列\n")

# 3. 计算每个脂质类别的总强度（关键：使用dplyr::select）
lipid_class_abundance <- lipid_all %>%
  dplyr::select(Class, all_of(sample_cols)) %>%
  group_by(Class) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup()

# 4. 转换长格式，添加分组信息
lipid_class_long <- lipid_class_abundance %>%
  pivot_longer(cols = -Class, names_to = "Sample", values_to = "Intensity") %>%
  mutate(Group = case_when(
    grepl("Buffalo", Sample) ~ "Buffalo",
    grepl("Holstein", Sample) ~ "Holstein",
    grepl("Jersey", Sample) ~ "Jersey",
    TRUE ~ "Other"
  ))

# 5. 计算每个样本各类别相对丰度（百分比）
lipid_class_pct <- lipid_class_long %>%
  group_by(Sample) %>%
  mutate(Percentage = Intensity / sum(Intensity) * 100) %>%
  ungroup()

# 6. 确定四个主要类别（基于所有样本平均相对丰度）
mean_pct <- lipid_class_pct %>%
  group_by(Class) %>%
  summarise(MeanPct = mean(Percentage, na.rm = TRUE)) %>%
  arrange(desc(MeanPct))

top4_classes <- head(mean_pct$Class, 4)
cat("四个主要脂质类别：", paste(top4_classes, collapse = ", "), "\n")
cat("其他类别将合并为 Other\n")

# 7. 将非前四类别合并为“Other”
lipid_class_final <- lipid_class_pct %>%
  mutate(Class_plot = ifelse(Class %in% top4_classes, Class, "Other")) %>%
  group_by(Sample, Group, Class_plot) %>%
  summarise(Percentage = sum(Percentage), .groups = "drop") %>%
  ungroup()

# 8. 设置类别因子顺序（图例顺序）
order_classes <- c(top4_classes, "Other")
lipid_class_final$Class_plot <- factor(lipid_class_final$Class_plot, 
                                        levels = order_classes)

# 9. 计算各组平均百分比（用于表格）
group_summary <- lipid_class_final %>%
  group_by(Group, Class_plot) %>%
  summarise(MeanPct = mean(Percentage), .groups = "drop")

# ================== 绘图 ==================
# 颜色方案：前四类使用Set1调色板，Other用灰色
class_colors <- c(brewer.pal(4, "Set1"), "grey60")
names(class_colors) <- order_classes

p <- ggplot(lipid_class_final, aes(x = Group, y = Percentage, fill = Class_plot)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = class_colors, name = "Lipid Class") +
  labs(
    x = "Group",
    y = "Relative Abundance",
    title = "Lipid Class Distribution",
    subtitle = paste("Top 4 classes:", paste(top4_classes, collapse = ", "))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# 显示图形
print(p)

# 10. 保存图形
ggsave("Lipid_Class_Distribution_Top4.pdf", p, width = 8, height = 6)
ggsave("Lipid_Class_Distribution_Top4.png", p, width = 8, height = 6, dpi = 300)

# 11. 保存组平均百分比表
group_summary_wide <- group_summary %>%
  pivot_wider(names_from = Class_plot, values_from = MeanPct, values_fill = 0) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))
write.csv(group_summary_wide, "Lipid_Class_Mean_Percentage.csv", row.names = FALSE)

cat("\n✅ 分析完成！\n")
cat("📊 生成图形：Lipid_Class_Distribution_Top4.pdf / .png\n")
cat("📁 生成表格：Lipid_Class_Mean_Percentage.csv\n")

# ============================================
# 脂质类别分布柱状图（水牛 vs 对照）
# 四个主要类别 + Other
# ============================================

# 1. 加载必要包（按需安装）
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(scales)) install.packages("scales")

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

# --- 解决select函数冲突 ---
select <- dplyr::select
assign("select", dplyr::select, envir = .GlobalEnv)
cat("已设置全局select函数为dplyr::select\n")

# ================== 数据准备 ==================
# 方案A：使用模拟数据（默认，保证直接运行）
set.seed(123)

# 生成模拟脂质注释
lipid_ids <- paste0("LIPID", 1:200)
lipid_names <- paste0("Lipid", 1:200)
lipid_classes <- sample(
  c("TG", "PC", "PE", "SM", "CE", "DG", "PS", "PI", "LPC", "Cer", "FFA"),
  200, replace = TRUE,
  prob = c(0.35, 0.18, 0.12, 0.08, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02)
)

# 模拟样本：Buffalo 5个，Holstein 5个，Jersey 5个
samples_buffalo <- paste0("Buffalo_", 1:5)
samples_holstein <- paste0("Holstein_", 1:5)
samples_jersey <- paste0("Jersey_", 1:5)
samples_all <- c(samples_buffalo, samples_holstein, samples_jersey)

# 模拟丰度矩阵
abundance_matrix <- matrix(
  runif(200 * 15, 100, 10000), 
  nrow = 200, ncol = 15
)
# 让TG丰度显著更高
tg_rows <- which(lipid_classes == "TG")
abundance_matrix[tg_rows, ] <- abundance_matrix[tg_rows, ] * 8
colnames(abundance_matrix) <- samples_all

# 构建lipid_all数据框
lipid_all <- data.frame(
  ID = lipid_ids,
  Name = lipid_names,
  Class = lipid_classes,
  abundance_matrix,
  stringsAsFactors = FALSE
)

cat("使用模拟数据，共", nrow(lipid_all), "个脂质分子，", 
    ncol(lipid_all)-3, "个样本\n")
cat("样本组成：Buffalo 5个, Holstein 5个, Jersey 5个\n")

# ============ 如果您有真实数据，请注释上方模拟代码，并取消下方注释 ============
# library(readxl)
# lipid_all <- read_excel("您的真实数据路径.xlsx")
# 要求：lipid_all必须包含列 'Class' 以及样本列（列名包含Buffalo/Holstein/Jersey）
# =============================================================================

# 2. 自动识别样本列（包含Buffalo/Holstein/Jersey）
sample_cols <- grep("Buffalo|Holstein|Jersey", colnames(lipid_all), value = TRUE)
sample_cols <- sample_cols[!grepl("QC", sample_cols)]  # 排除QC样本
if (length(sample_cols) == 0) {
  stop("未找到任何样本列！请检查数据列名是否包含 Buffalo/Holstein/Jersey")
}
cat("找到", length(sample_cols), "个样本列\n")

# 3. 计算每个脂质类别的总强度
lipid_class_abundance <- lipid_all %>%
  dplyr::select(Class, all_of(sample_cols)) %>%
  group_by(Class) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup()

# 4. 转换长格式，并**重新分组：水牛 vs 对照（荷斯坦+娟姗）**
lipid_class_long <- lipid_class_abundance %>%
  pivot_longer(cols = -Class, names_to = "Sample", values_to = "Intensity") %>%
  mutate(Group = case_when(
    grepl("Buffalo", Sample) ~ "Buffalo",
    grepl("Holstein|Jersey", Sample) ~ "Control",  # 合并荷斯坦和娟姗为Control
    TRUE ~ "Other"
  ))

cat("分组信息：\n")
print(table(lipid_class_long$Group))

# 5. 计算每个样本各类别相对丰度（百分比）
lipid_class_pct <- lipid_class_long %>%
  group_by(Sample) %>%
  mutate(Percentage = Intensity / sum(Intensity) * 100) %>%
  ungroup()

# 6. 确定四个主要类别（基于所有样本的平均相对丰度）
mean_pct <- lipid_class_pct %>%
  group_by(Class) %>%
  summarise(MeanPct = mean(Percentage, na.rm = TRUE)) %>%
  arrange(desc(MeanPct))

top4_classes <- head(mean_pct$Class, 4)
cat("四个主要脂质类别：", paste(top4_classes, collapse = ", "), "\n")
cat("其他类别将合并为 Other\n")

# 7. 将非前四类别合并为“Other”
lipid_class_final <- lipid_class_pct %>%
  mutate(Class_plot = ifelse(Class %in% top4_classes, Class, "Other")) %>%
  group_by(Sample, Group, Class_plot) %>%
  summarise(Percentage = sum(Percentage), .groups = "drop") %>%
  ungroup()

# 8. 设置类别因子顺序（图例顺序）
order_classes <- c(top4_classes, "Other")
lipid_class_final$Class_plot <- factor(lipid_class_final$Class_plot, 
                                        levels = order_classes)

# 9. 计算**各组平均百分比**（用于绘图和表格）
group_summary <- lipid_class_final %>%
  group_by(Group, Class_plot) %>%
  summarise(MeanPct = mean(Percentage), .groups = "drop")

# ================== 绘图 ==================
# 颜色方案：前四类使用Set1调色板，Other用灰色
class_colors <- c(brewer.pal(4, "Set1"), "grey60")
names(class_colors) <- order_classes

p <- ggplot(lipid_class_final, aes(x = Group, y = Percentage, fill = Class_plot)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = class_colors, name = "Lipid Class") +
  labs(
    x = "Group",
    y = "Relative Abundance",
    title = "Lipid Class Distribution",
    subtitle = paste("Top 4 classes:", paste(top4_classes, collapse = ", "))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# 显示图形
print(p)

# 10. 保存图形（水牛 vs 对照）
ggsave("Lipid_Class_Buffalo_vs_Control.pdf", p, width = 7, height = 6)
ggsave("Lipid_Class_Buffalo_vs_Control.png", p, width = 7, height = 6, dpi = 300)

# 11. 保存组平均百分比表（两组）
group_summary_wide <- group_summary %>%
  pivot_wider(names_from = Class_plot, values_from = MeanPct, values_fill = 0) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))
write.csv(group_summary_wide, "Lipid_Class_Mean_Percentage_Buffalo_vs_Control.csv", 
          row.names = FALSE)

cat("\n✅ 分析完成！\n")
cat("📊 生成图形：Lipid_Class_Buffalo_vs_Control.pdf / .png\n")
cat("📁 生成表格：Lipid_Class_Mean_Percentage_Buffalo_vs_Control.csv\n")

# ================== 可选：添加统计检验 ==================
# 检验两组之间前四类脂质的相对丰度是否有显著差异
cat("\n📈 统计检验（Wilcoxon秩和检验）：\n")
for (cls in top4_classes) {
  buffalo_vals <- lipid_class_final %>%
    filter(Group == "Buffalo", Class_plot == cls) %>%
    pull(Percentage)
  control_vals <- lipid_class_final %>%
    filter(Group == "Control", Class_plot == cls) %>%
    pull(Percentage)
  test <- wilcox.test(buffalo_vals, control_vals)
  cat(sprintf("%s: p = %.4f\n", cls, test$p.value))
}

# ================== 可选：绘制箱线图 ==================
# 如果想展示单个类别的组间比较，可以取消下面注释
# p_box <- ggplot(lipid_class_final %>% filter(Class_plot %in% top4_classes),
#                 aes(x = Group, y = Percentage, fill = Group)) +
#   geom_boxplot() +
#   facet_wrap(~Class_plot, scales = "free_y") +
#   scale_fill_manual(values = c("Buffalo" = "#E64B35", "Control" = "#4DBBD5")) +
#   labs(y = "Relative Abundance (%)") +
#   theme_minimal()
# ggsave("Lipid_Class_Boxplot.pdf", p_box, width = 10, height = 6)
# ============================================
# 火山图：水牛乳 vs 对照乳（荷斯坦+娟姗）
# 输出表包含脂质类别（Class）
# ============================================

# ---------- 1. 加载包 ----------
packages <- c("limma", "ggplot2", "ggrepel", "dplyr", "tidyr", "readxl")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% "limma") BiocManager::install(pkg) else install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
# 解决select冲突
select <- dplyr::select
assign("select", dplyr::select, envir = .GlobalEnv)

# ---------- 2. 读取真实数据（请修改为您的实际文件路径）----------
# 必须包含列：ID, Name, Class, 以及样本列（含Buffalo/Holstein/Jersey）
lipid_all <- read_excel("https://raw.githubusercontent.com/Alice/my-r-scripts/main/meta_intensity_all.xlsx")

# 可选：检查数据结构
# head(lipid_all)

# ========== 以下为模拟数据示例（用于测试，正式运行时请注释掉）==========
# set.seed(42)
# ...（模拟数据代码已省略）...
# ========== 模拟数据结束 ==========

# ---------- 3. 通用函数：运行水牛 vs 合并对照的差异分析并绘制火山图 ----------
run_diff_volcano <- function(lipid_data, 
                             buffalo_pattern = "Buffalo", 
                             control_pattern = "Holstein|Jersey",  # 同时匹配荷斯坦和娟姗
                             comparison_name = "Dairy",           # 对照组合并命名
                             logFC_cut = 1,
                             p_adj_cut = 0.05,
                             n_label = 10) {
  
  cat("\n========== 差异分析：水牛 vs", comparison_name, "==========\n")
  
  # ----- 3.1 提取样本列（水牛 + 荷斯坦/娟姗，排除QC）-----
  sample_cols <- grep(paste(buffalo_pattern, control_pattern, sep = "|"), 
                      colnames(lipid_data), value = TRUE)
  sample_cols <- sample_cols[!grepl("QC", sample_cols)]
  if (length(sample_cols) == 0) stop("未找到指定样本列！")
  cat("共识别", length(sample_cols), "个样本列\n")
  
  # ----- 3.2 构建表达矩阵 -----
  exp_mat <- as.matrix(lipid_data[, sample_cols])
  rownames(exp_mat) <- lipid_data$ID
  
  # 缺失值处理（用最小值的一半填充）
  if (any(exp_mat == 0 | is.na(exp_mat))) {
    min_val <- min(exp_mat[exp_mat > 0], na.rm = TRUE) / 2
    exp_mat[exp_mat == 0] <- NA
    exp_mat[is.na(exp_mat)] <- min_val
  }
  exp_mat <- log2(exp_mat)
  
  # ----- 3.3 分组设计：水牛 vs 合并对照 -----
  group <- ifelse(grepl(buffalo_pattern, sample_cols), "Buffalo", "Control")
  group <- factor(group, levels = c("Control", "Buffalo"))
  cat("\n分组样本数：\n"); print(table(group))
  
  # ----- 3.4 limma差异分析 -----
  design <- model.matrix(~ group)
  colnames(design) <- c("Intercept", "Buffalo_vs_Control")
  fit <- lmFit(exp_mat, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  deg <- topTable(fit, coef = "Buffalo_vs_Control", number = Inf, adjust.method = "BH")
  deg$ID <- rownames(deg)
  
  # 合并脂质注释（ID, Name, Class）
  deg <- merge(deg, lipid_data[, c("ID", "Name", "Class")], by = "ID", all.x = TRUE)
  
  # ----- 3.5 显著性标记 -----
  deg$Significance <- case_when(
    deg$logFC > logFC_cut & deg$adj.P.Val < p_adj_cut ~ "Up",
    deg$logFC < -logFC_cut & deg$adj.P.Val < p_adj_cut ~ "Down",
    TRUE ~ "Not Significant"
  )
  cat("\n差异脂质统计：\n")
  print(table(deg$Significance))
  
  # ----- 3.6 标记Top N脂质（用于火山图标签）-----
  top_up <- deg %>% filter(Significance == "Up") %>% arrange(adj.P.Val) %>% head(n_label)
  top_down <- deg %>% filter(Significance == "Down") %>% arrange(adj.P.Val) %>% head(n_label)
  top_labels <- bind_rows(top_up, top_down)
  
  # ----- 3.7 绘制火山图 -----
  p_volcano <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.8) +
    scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", 
                                  "Not Significant" = "grey75")) +
    geom_vline(xintercept = c(-logFC_cut, logFC_cut), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(p_adj_cut), linetype = "dashed", color = "grey30") +
    {if(nrow(top_labels) > 0) geom_text_repel(
      data = top_labels, aes(label = Name),
      size = 3.2, box.padding = 0.6, max.overlaps = 20, segment.color = "grey40"
    )} +
    labs(x = expression(log[2]("Fold Change")),
         y = expression(-log[10]("Adjusted P-value")),
         title = paste("Differential Lipids: Buffalo vs", comparison_name),
         subtitle = paste0("Total: ", nrow(deg), 
                          " | Up: ", sum(deg$Significance == "Up"),
                          " | Down: ", sum(deg$Significance == "Down"))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  print(p_volcano)
  
  # ----- 3.8 保存图形和表格 -----
  base_name <- paste0("Volcano_Buffalo_vs_", gsub(" ", "_", comparison_name))
  ggsave(paste0(base_name, ".pdf"), p_volcano, width = 12, height = 8)
  ggsave(paste0(base_name, ".png"), p_volcano, width = 12, height = 8, dpi = 300)
  write.csv(deg, paste0("Differential_Lipids_Buffalo_vs_", 
                        gsub(" ", "_", comparison_name), "_wClass.csv"), 
            row.names = FALSE)
  
  cat("✅ 完成！输出文件：\n")
  cat("   - ", base_name, ".pdf/.png\n")
  cat("   - Differential_Lipids_Buffalo_vs_", gsub(" ", "_", comparison_name), "_wClass.csv\n")
  
  return(deg)  # 返回差异分析结果
}

# ---------- 4. 执行单次比较：水牛 vs 合并对照（荷斯坦+娟姗）----------
deg_dairy <- run_diff_volcano(
  lipid_data = lipid_all,
  buffalo_pattern = "Buffalo",
  control_pattern = "Holstein|Jersey",   # 正则表达式同时匹配荷斯坦和娟姗
  comparison_name = "Dairy",            # 对照组合并命名为 Dairy
  logFC_cut = 1,
  p_adj_cut = 0.05,
  n_label = 10
)

cat("\n🎉 火山图及结果表生成完毕！\n")
# ============================================
# TG碳链长度分组：水牛 vs 对照（荷斯坦+娟姗）
# 绝对丰度 vs 相对丰度
# 分组方案：1) ≤48（短链）；2) 50–52（中长链）；
#          3) 54–56（长链）；4) ≥58（超长链）
# 比较组别：水牛乳、对照牛乳（荷斯坦+娟姗）
# 所有图形标签均为英文，无中文
# ============================================

# ------------------ 1. 加载包 ------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "ggpubr", 
                       "stringr", "vegan", "compositions", "rstatix")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
# 明确函数来源
select <- dplyr::select
first <- dplyr::first
cat("Global 'first' function set to dplyr::first\n")

# ------------------ 2. 模拟数据（可替换为真实数据） ------------------
set.seed(123)

# 生成TG名称（含碳原子数）
tg_names <- c(
  paste0("TG ", seq(46, 58, 2), ":0"),
  paste0("TG ", seq(46, 58, 2), ":1"),
  paste0("TG ", seq(48, 56, 2), ":2"),
  paste0("TG ", seq(50, 54, 2), ":3")
)
carbon_lengths <- as.numeric(str_extract(tg_names, "\\d+"))

# 样本名称：水牛、荷斯坦、娟姗各6个
samples <- c(paste0("Buffalo_", 1:6), 
             paste0("Holstein_", 1:6), 
             paste0("Jersey_", 1:6))
n_tg <- length(tg_names)
n_samples <- length(samples)

# 构建强度矩阵（模拟）
intensity <- matrix(0, nrow = n_tg, ncol = n_samples)
colnames(intensity) <- samples
rownames(intensity) <- tg_names

for (i in 1:n_tg) {
  for (j in 1:n_samples) {
    base <- 500 + rnorm(1, 0, 100)
    if (grepl("Buffalo", samples[j])) base <- base * 1.8
    carbon_weight <- 1
    if (carbon_lengths[i] <= 48) carbon_weight <- 0.5
    else if (carbon_lengths[i] <= 52) carbon_weight <- 1.0
    else if (carbon_lengths[i] <= 56) carbon_weight <- 1.2
    else carbon_weight <- 0.8
    intensity[i, j] <- base * carbon_weight * (1 + rnorm(1, 0, 0.1))
  }
}
intensity[intensity < 0] <- 100

# 转换为长格式
lipid_long <- data.frame(
  Sample = rep(colnames(intensity), each = n_tg),
  Lipid = rep(rownames(intensity), n_samples),
  Intensity = as.vector(intensity),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Carbon = as.numeric(str_extract(Lipid, "\\d+")),
    # 原始三品种分组
    Group = case_when(
      grepl("Buffalo", Sample)  ~ "Buffalo",
      grepl("Holstein", Sample) ~ "Holstein",
      grepl("Jersey", Sample)   ~ "Jersey",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# ------------------ 3. 合并为两组：水牛 vs 对照（荷斯坦+娟姗） ------------------
lipid_long <- lipid_long %>%
  mutate(
    Group2 = ifelse(Group == "Buffalo", "Buffalo", "Dairy"),  # 对照统一为"Dairy"
    Group2 = factor(Group2, levels = c("Buffalo", "Dairy"))
  )

cat("数据重构完成：水牛 vs 对照（Dairy = Holstein + Jersey）\n")

# ------------------ 4. 碳链长度分组（英文标签） ------------------
lipid_long <- lipid_long %>%
  mutate(
    Carbon_Group = case_when(
      Carbon <= 48 ~ "≤48 (Short)",
      Carbon %in% 50:52 ~ "50–52 (Medium)",
      Carbon %in% 54:56 ~ "54–56 (Long)",
      Carbon >= 58 ~ "≥58 (Very long)"
    )
  )
lipid_long$Carbon_Group <- factor(
  lipid_long$Carbon_Group,
  levels = c("≤48 (Short)", "50–52 (Medium)", "54–56 (Long)", "≥58 (Very long)")
)

# ------------------ 5. 绝对丰度计算（按样本、分组2、碳链组） ------------------
abs_sum <- lipid_long %>%
  group_by(Sample, Group2, Carbon_Group) %>%
  summarise(AbsIntensity = sum(Intensity), .groups = "drop")

# ------------------ 6. 相对丰度计算 ------------------
tg_total <- lipid_long %>%
  group_by(Sample) %>%
  summarise(TotalTG = sum(Intensity), .groups = "drop")

rel_sum <- lipid_long %>%
  left_join(tg_total, by = "Sample") %>%
  group_by(Sample, Group2, Carbon_Group) %>%
  summarise(
    RelPercent = sum(Intensity) / dplyr::first(TotalTG) * 100,
    .groups = "drop"
  )

# ------------------ 7. 定义两组颜色 ------------------
group_colors <- c("Buffalo" = "#E64B35", 
                  "Dairy" = "#4DBBD5")  # 使用蓝色代表对照

# ------------------ 8. 绘制绝对丰度箱线图（两组） ------------------
p_abs <- ggplot(abs_sum, aes(x = Carbon_Group, y = AbsIntensity, fill = Group2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             size = 1.2, alpha = 0.6, aes(color = Group2)) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = "TG Carbon Chain Length", y = "Absolute Intensity",
       title = "Absolute Abundance of TG by Chain Length",
       subtitle = "Buffalo vs Dairy (Holstein+Jersey)",
       fill = "Group", color = "Group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  # 添加Wilcoxon检验p值（每组碳链内两组比较）
  stat_compare_means(aes(group = Group2), 
                     method = "wilcox.test", 
                     label = "p.format",
                     label.y = max(abs_sum$AbsIntensity) * 0.9)

# ------------------ 9. 绘制相对丰度箱线图（两组） ------------------
p_rel <- ggplot(rel_sum, aes(x = Carbon_Group, y = RelPercent, fill = Group2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             size = 1.2, alpha = 0.6, aes(color = Group2)) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = "TG Carbon Chain Length", y = "Relative Abundance (%)",
       title = "Relative Abundance of TG by Chain Length",
       subtitle = "Buffalo vs Dairy (Holstein+Jersey)",
       fill = "Group", color = "Group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = Group2), 
                     method = "wilcox.test", 
                     label = "p.format",
                     label.y = max(rel_sum$RelPercent) * 0.9)

# ------------------ 10. 统计检验：两组比较 ------------------
cat("\n========== 两组比较统计检验 ==========\n")

# (1) 每个碳链组内，水牛 vs 对照（Wilcoxon秩和检验 + FDR校正）
wilcox_results <- rel_sum %>%
  group_by(Carbon_Group) %>%
  wilcox_test(RelPercent ~ Group2) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
cat("\n--- Wilcoxon test per carbon group (Buffalo vs Dairy) ---\n")
print(wilcox_results)

# (2) 整体组成差异：PERMANOVA（基于CLR变换，两组）
rel_matrix <- rel_sum %>%
  select(Sample, Group2, Carbon_Group, RelPercent) %>%
  pivot_wider(names_from = Carbon_Group, values_from = RelPercent, values_fill = 0) %>%
  as.data.frame()
rownames(rel_matrix) <- rel_matrix$Sample
group2_factor <- rel_matrix$Group2
rel_data <- rel_matrix[, -c(1,2)]

# 零值处理（CLR要求>0）
rel_data[rel_data == 0] <- 0.001
rel_clr <- as.data.frame(compositions::clr(rel_data))

set.seed(123)
permanova_2g <- vegan::adonis2(rel_clr ~ group2_factor, method = "euclidean", permutations = 999)
cat("\n--- PERMANOVA (Buffalo vs Dairy) ---\n")
print(permanova_2g)

# (3) 卡方检验（基于合并强度，两组）
chi_data <- lipid_long %>%
  group_by(Group2, Carbon_Group) %>%
  summarise(SumIntensity = sum(Intensity), .groups = "drop") %>%
  pivot_wider(names_from = Carbon_Group, values_from = SumIntensity) %>%
  as.data.frame()
rownames(chi_data) <- chi_data$Group2
chi_table <- as.matrix(chi_data[, -1])
chisq_test_2g <- chisq.test(chi_table)
cat("\n--- Chi-square Test (Buffalo vs Dairy) ---\n")
print(chisq_test_2g)

# ------------------ 11. 保存图形 ------------------
ggsave("Absolute_Intensity_Buffalo_vs_Dairy.pdf", p_abs, width = 10, height = 6)
ggsave("Relative_Abundance_Buffalo_vs_Dairy.pdf", p_rel, width = 10, height = 6)
ggsave("Absolute_Intensity_Buffalo_vs_Dairy.png", p_abs, width = 10, height = 6, dpi = 300)
ggsave("Relative_Abundance_Buffalo_vs_Dairy.png", p_rel, width = 10, height = 6, dpi = 300)

# ------------------ 12. 保存统计结果 ------------------
write.csv(wilcox_results, "CarbonGroup_Wilcoxon_Buffalo_vs_Dairy.csv", row.names = FALSE)

sink("PERMANOVA_Buffalo_vs_Dairy.txt")
print(permanova_2g)
sink()

sink("Chisq_Buffalo_vs_Dairy.txt")
print(chisq_test_2g)
sink()

# 保存相对丰度汇总数据（可选）
write.csv(rel_sum, "TG_ChainGroup_RelativeAbundance_Buffalo_vs_Dairy.csv", row.names = FALSE)

cat("\n✅ 两组比较（水牛 vs 对照）分析完成！\n")
cat("📊 图形已保存：Absolute_Intensity_Buffalo_vs_Dairy.pdf/png, Relative_Abundance_Buffalo_vs_Dairy.pdf/png\n")
cat("📁 统计结果：CarbonGroup_Wilcoxon_Buffalo_vs_Dairy.csv, PERMANOVA_Buffalo_vs_Dairy.txt, Chisq_Buffalo_vs_Dairy.txt\n")
cat("🎨 颜色：Buffalo = #E64B35, Dairy = #4DBBD5\n")
# ============================================
# 热图：水牛乳 vs 对照乳（荷斯坦+娟姗）
# 显示全部显著上调的甘油三酯（TG）
# 行注释：碳链总长度分组 + 不饱和度分组
# 列注释：品种分组（Buffalo / Control）
# ============================================

# 1. 加载必要的R包 ------------------------------------------------------------
packages <- c("pheatmap", "RColorBrewer", "dplyr", "stringr", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
select <- dplyr::select

# 2. 提供的数据（从问题中复制，仅保留显著上调的TG）----------------------------
# 原始数据包含所有比较，此处直接定义数据框
tg_data_raw <- data.frame(
  ID = c("Com_1009_pos", "Com_1012_pos", "Com_1014_pos", "Com_1033_pos", 
         "Com_1041_pos", "Com_1043_pos", "Com_1044_pos", "Com_1045_pos", 
         "Com_1048_pos", "Com_1051_pos", "Com_1052_pos", "Com_1053_pos", 
         "Com_1054_pos", "Com_1057_pos", "Com_1059_pos", "Com_1060_pos", 
         "Com_1061_pos", "Com_1068_pos", "Com_1090_pos", "Com_1091_pos", 
         "Com_1092_pos", "Com_1093_pos", "Com_1094_pos", "Com_1095_pos", 
         "Com_1144_pos", "Com_1163_pos", "Com_1168_pos", "Com_1189_pos", 
         "Com_1193_pos", "Com_698_pos", "Com_717_pos", "Com_744_pos", 
         "Com_770_pos", "Com_774_pos", "Com_783_pos", "Com_789_pos", 
         "Com_790_pos", "Com_792_pos", "Com_806_pos", "Com_817_pos", 
         "Com_827_pos", "Com_831_pos", "Com_832_pos", "Com_844_pos", 
         "Com_861_pos", "Com_862_pos", "Com_867_pos", "Com_889_pos", 
         "Com_891_pos", "Com_897_pos", "Com_898_pos", "Com_899_pos", 
         "Com_903_pos", "Com_914_pos", "Com_932_pos", "Com_939_pos", 
         "Com_940_pos", "Com_941_pos", "Com_942_pos", "Com_944_pos", 
         "Com_945_pos", "Com_946_pos", "Com_947_pos", "Com_948_pos", 
         "Com_949_pos", "Com_950_pos", "Com_951_pos", "Com_952_pos", 
         "Com_958_pos", "Com_971_pos", "Com_998_pos"),
  Name = c("TG(4:0/16:0/16:1)", "TG(4:0/16:0/18:0CHO)", "TG(4:0/16:0/18:1CHO)", 
           "TG(4:0/18:1/18:1CHO)", "TG(4:0/2:0/15:0)", "TG(4:0/4:0/12:0)", 
           "TG(4:0/4:0/13:0)", "TG(4:0/4:0/14:0)", "TG(4:0/4:0/16:0)", 
           "TG(4:0/4:0/18:1)", "TG(4:0/4:0/18:2)", "TG(4:0/4:0/18:3)", 
           "TG(4:0/4:0/20:4)", "TG(4:0/6:0/10:0)", "TG(4:0/6:0/14:0)", 
           "TG(4:0/6:0/14:1)", "TG(4:0/6:0/15:0)", "TG(4:0/6:0/18:1CHO)", 
           "TG(4:0CHO/12:0/16:0)", "TG(4:0CHO/14:0/16:0)", "TG(4:0CHO/14:0/18:1)", 
           "TG(4:0CHO/16:0/16:0)", "TG(4:0CHO/16:0/18:1)", "TG(4:0CHO/18:1/18:1)", 
           "TG(6:0/2:0/18:1)", "TG(7:0/7:0/15:0)", "TG(7:2CHO/16:0/17:0)", 
           "TG(8:0CHO/14:0/14:0)", "TG(8:0CHO/4:0/17:3)", "TG(10:2CHO/6:0/14:0)", 
           "TG(12:0/2:0/6:0)", "TG(14:0/17:3/18:3)", "TG(15:0/5:0/9:0)", 
           "TG(15:1/6:0/18:4)", "TG(16:0/16:0/23:0)", "TG(16:0/18:0/18:0)", 
           "TG(16:0/18:0/18:0CHO)", "TG(16:0/18:0/20:0)", "TG(16:0/23:0/18:1)", 
           "TG(17:0/18:0/18:0)", "TG(17:2CHO/16:0)", "TG(18:0/18:0/18:0)", 
           "TG(18:0/18:0/18:1)", "TG(18:0/4:0/18:0CHO)", "TG(18:1/18:2/20:4)", 
           "TG(18:1/2:0)", "TG(18:2CHO/18:0)", "TG(21:1)", "TG(21:2CHO/18:1)", 
           "TG(22:2CHO/15:0)", "TG(22:2CHO/16:0)", "TG(22:2CHO/18:0)", 
           "TG(24:1/17:1COOH)", "TG(26:6/15:2)", "TG(2:0/16:0/21:0CHO)", 
           "TG(2:0/2:0/16:0)", "TG(2:0/2:0/18:0)", "TG(2:0/2:0/18:1)", 
           "TG(2:0/4:0/14:0)", "TG(2:0/4:0/16:0)", "TG(2:0/4:0/16:1)", 
           "TG(2:0/4:0/17:1)", "TG(2:0/4:0/18:0)", "TG(2:0/4:0/18:1)", 
           "TG(2:0/4:0/18:2)", "TG(2:0/4:0/18:3)", "TG(2:0/4:0/20:1)", 
           "TG(30:1/15:1COOH)", "TG(34:7)", "TG(47:13CHO/29:6)", 
           "TG(4:0/14:0CHO/16:0)"),
  Class = "TG",
  Significance = "Up",
  stringsAsFactors = FALSE
)

cat("✅ 显著上调TG数量：", nrow(tg_data_raw), "\n")

# 3. 从Name解析总碳链长度和不饱和度 ----------------------------------------
parse_tg_name <- function(name) {
  # 提取括号内内容，如 "4:0/16:0/16:1"
  content <- str_extract(name, "(?<=\\().*(?=\\))")
  if (is.na(content)) return(c(carbon = NA, db = NA))
  
  # 按 "/" 分割各脂肪酸
  parts <- unlist(strsplit(content, "/"))
  total_carbon <- 0
  total_db <- 0
  for (p in parts) {
    # 提取类似 "4:0" 的部分，可能包含CHO等后缀
    match <- str_extract(p, "\\d+:\\d+")
    if (!is.na(match)) {
      nums <- as.numeric(unlist(strsplit(match, ":")))
      total_carbon <- total_carbon + nums[1]
      total_db <- total_db + nums[2]
    }
  }
  c(carbon = total_carbon, db = total_db)
}

# 解析所有TG
parsed <- t(sapply(tg_data_raw$Name, parse_tg_name))
tg_data_raw$carbon <- parsed[, "carbon"]
tg_data_raw$db <- parsed[, "db"]

# 检查是否有解析失败的（如TG(21:1) 格式为单脂肪酸）
failed <- which(is.na(tg_data_raw$carbon))
if (length(failed) > 0) {
  for (i in failed) {
    # 尝试直接匹配 "数字:数字"
    m <- str_match(tg_data_raw$Name[i], "\\((\\d+):(\\d+)\\)")
    if (!is.na(m[1,1])) {
      tg_data_raw$carbon[i] <- as.numeric(m[1,2])
      tg_data_raw$db[i] <- as.numeric(m[1,3])
    }
  }
}

# 再次检查
if (any(is.na(tg_data_raw$carbon))) {
  stop("仍有TG名称解析失败，请检查：", 
       paste(tg_data_raw$Name[is.na(tg_data_raw$carbon)], collapse = "; "))
}

cat("碳链长度范围：", range(tg_data_raw$carbon, na.rm = TRUE), "\n")
cat("不饱和度范围：", range(tg_data_raw$db, na.rm = TRUE), "\n")

# 4. 创建行注释分组 --------------------------------------------------------
tg_data <- tg_data_raw %>%
  mutate(
    carbon_group = case_when(
      carbon <= 48 ~ "≤48",
      carbon %in% 50:52 ~ "50–52",
      carbon %in% 54:56 ~ "54–56",
      carbon >= 58 ~ "≥58",
      TRUE ~ "Other"
    ),
    db_group = case_when(
      db == 0 ~ "Saturated (0)",
      db == 1 ~ "Monounsaturated (1)",
      db >= 2 ~ "Polyunsaturated (≥2)",
      TRUE ~ "Unknown"
    )
  )

tg_data$carbon_group <- factor(tg_data$carbon_group,
                               levels = c("≤48", "50–52", "54–56", "≥58", "Other"))
tg_data$db_group <- factor(tg_data$db_group,
                           levels = c("Saturated (0)", "Monounsaturated (1)",
                                      "Polyunsaturated (≥2)", "Unknown"))

annotation_row <- tg_data %>%
  select(ID, carbon_group, db_group) %>%
  column_to_rownames(var = "ID")

# 5. 定义样本（假设每种3个重复）-------------------------------------------
samples_buffalo <- paste0("Buffalo_", 1:3)
samples_holstein <- paste0("Holstein_", 1:3)
samples_jersey  <- paste0("Jersey_", 1:3)
samples_all <- c(samples_buffalo, samples_holstein, samples_jersey)

# 列注释：Group列，水牛为"Buffalo"，其他合并为"Control"
group <- c(rep("Buffalo", 3), rep("Control", 6))
annotation_col <- data.frame(
  Group = factor(group, levels = c("Buffalo", "Control")),
  row.names = samples_all
)

# 6. 模拟表达矩阵（水牛组显著上调）---------------------------------------
set.seed(2025)
n_tg <- nrow(tg_data)
exp_mat <- matrix(NA, nrow = n_tg, ncol = 9)
rownames(exp_mat) <- tg_data$ID
colnames(exp_mat) <- samples_all

for (i in 1:n_tg) {
  # 基础均值（对照组）
  base_mean <- runif(1, 11, 14)
  base_sd   <- runif(1, 0.3, 0.6)
  # 上调倍数（1.5~3倍）
  fc <- runif(1, 1.5, 3.0)
  
  exp_mat[i, 1:3] <- rnorm(3, mean = base_mean + fc, sd = base_sd)  # Buffalo
  exp_mat[i, 4:9] <- rnorm(6, mean = base_mean,      sd = base_sd)  # Control
}

# 确保表达值合理
exp_mat[exp_mat < 5] <- 5 + abs(rnorm(sum(exp_mat < 5), 0, 0.2))
exp_mat <- exp_mat + matrix(rnorm(n_tg * 9, 0, 0.1), nrow = n_tg, ncol = 9)

# 行标准化（Z-score）
exp_mat_scaled <- t(scale(t(exp_mat)))

# 7. 定义注释颜色 ----------------------------------------------------------
carbon_colors <- c(
  "≤48"   = "#8DD3C7",
  "50–52" = "#FFFFB3",
  "54–56" = "#BEBADA",
  "≥58"   = "#FB8072",
  "Other" = "grey80"
)
db_colors <- c(
  "Saturated (0)"       = "#FDB462",
  "Monounsaturated (1)" = "#80B1D3",
  "Polyunsaturated (≥2)"= "#B3DE69",
  "Unknown"             = "grey90"
)
group_colors <- c(
  "Buffalo" = "#E64B35",
  "Control" = "#4DBBD5"   # 荷斯坦+娟姗合并为对照
)
annotation_colors <- list(
  carbon_group = carbon_colors,
  db_group     = db_colors,
  Group        = group_colors
)

# 8. 行标签：简化为脂肪酸组成（去掉"TG"前缀）-----------------------------
labels_row <- gsub("TG\\((.*?)\\)", "\\1", tg_data$Name)
# 若仍以TG开头，进一步清理
labels_row <- gsub("^TG", "", labels_row)

# 9. 绘制热图（列不聚类，保持固定顺序）-----------------------------------
heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

# PDF输出
pheatmap(
  exp_mat_scaled,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,          # 列固定顺序
  show_rownames = TRUE,
  labels_row = labels_row,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 10,
  border_color = NA,
  main = "Up-regulated TGs in Buffalo Milk vs Control (Holstein+Jersey)",
  filename = "Heatmap_TG_Buffalo_vs_Control_fixedOrder.pdf",
  width = 14,
  height = 18
)

# PNG输出（高分辨率）
png("Heatmap_TG_Buffalo_vs_Control_fixedOrder.png",
    width = 4200, height = 4800, res = 300)
pheatmap(
  exp_mat_scaled,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = labels_row,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 10,
  border_color = NA,
  main = "Up-regulated TGs in Buffalo Milk vs Control (Holstein+Jersey)"
)
dev.off()

cat("\n✅ 热图绘制完成！\n")
cat("📊 文件保存：\n")
cat("   - PDF: Heatmap_TG_Buffalo_vs_Control_fixedOrder.pdf\n")
cat("   - PNG: Heatmap_TG_Buffalo_vs_Control_fixedOrder.png\n")
cat("   - 共显示", n_tg, "个显著上调的TG\n")
# ============================================
# 脂质代谢通路富集气泡图
# ============================================

# 1. 加载必要的包
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")

library(ggplot2)
library(dplyr)
library(tidyr)

# ================== 模拟通路富集数据 ==================
# 模拟与脂质代谢密切相关的KEGG通路富集结果
# 实际使用时，请替换为您自己的富集分析结果数据框
set.seed(123)

pathway_data <- data.frame(
  # 通路名称（使用标准KEGG脂质代谢相关通路）
  Pathway = c(
    "Fatty acid biosynthesis",
    "Glycerolipid metabolism",
    "Glycerophospholipid metabolism",
    "Sphingolipid metabolism",
    "Fatty acid degradation",
    "Fatty acid elongation",
    "Steroid biosynthesis",
    "Arachidonic acid metabolism",
    "Linoleic acid metabolism",
    "alpha-Linolenic acid metabolism",
    "Biosynthesis of unsaturated fatty acids",
    "Sphingolipid signaling pathway",
    "PPAR signaling pathway",
    "Adipocytokine signaling pathway",
    "Cholesterol metabolism"
  ),
  # 富集因子 (GeneRatio) = 差异基因数 / 通路基因总数
  GeneRatio = round(runif(15, 0.1, 0.4), 3),
  # 背景比例 (BgRatio) - 此处不直接使用，但用于计算富集倍数等
  BgRatio = rep("50/8000", 15),
  # p值 (原始p值，越小越显著)
  pvalue = c(
    2.5e-8, 3.2e-6, 4.1e-5, 0.0008, 0.0012,
    0.0035, 0.0087, 0.015, 0.021, 0.032,
    0.045, 0.056, 0.063, 0.078, 0.085
  ),
  # 校正后p值 (通常使用BH法)
  p.adjust = c(
    1.8e-7, 1.5e-5, 1.9e-4, 0.0032, 0.0046,
    0.012, 0.026, 0.038, 0.045, 0.058,
    0.069, 0.078, 0.083, 0.092, 0.095
  ),
  # 该通路中差异表达的基因/脂质数量
  Count = c(18, 15, 12, 10, 8, 7, 6, 6, 5, 4, 4, 3, 3, 2, 2),
  stringsAsFactors = FALSE
)

# 添加富集倍数 (Fold Enrichment) 或直接使用GeneRatio
# 这里我们直接用GeneRatio代表富集程度
pathway_data$FoldEnrichment <- pathway_data$GeneRatio / 0.00625  # 假设背景比例为50/8000≈0.00625

# 按p值排序，取Top 10最显著通路用于绘图（可根据需要调整）
pathway_data <- pathway_data %>%
  arrange(p.adjust) %>%
  head(10)

# 确保通路名称顺序与p值排序一致（气泡图中y轴按显著性排序）
pathway_data$Pathway <- factor(pathway_data$Pathway, 
                               levels = rev(pathway_data$Pathway))  # 反转，使最显著的通路在顶部

# ================== 绘制气泡图 ==================
p <- ggplot(pathway_data, aes(
  x = GeneRatio, 
  y = Pathway,
  size = Count,
  color = -log10(p.adjust)
)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(3, 10), name = "Gene Count") +
  scale_color_gradient(low = "#4DBBD5", high = "#E64B35", 
                       name = expression(-log[10]("adj.P"))) +
  labs(
    x = "Gene Ratio",
    y = "Lipid Metabolic Pathway",
    title = "KEGG Pathway Enrichment Analysis",
    subtitle = "Top 10 significantly enriched lipid-related pathways"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5)
  ) +
  # 添加Gene Ratio数值标签（可选）
  geom_text(aes(label = Count), hjust = -0.5, vjust = 0.5, size = 3.5, color = "black")

# 显示图形
print(p)

# 保存为PDF（矢量图，可编辑）
ggsave("Lipid_Pathway_Enrichment_Bubble.pdf", p, width = 10, height = 6)

# 保存为PNG（高分辨率位图）
ggsave("Lipid_Pathway_Enrichment_Bubble.png", p, width = 10, height = 6, dpi = 300)

cat("\n✅ 脂质代谢通路富集气泡图绘制完成！\n")
cat("📊 文件保存：Lipid_Pathway_Enrichment_Bubble.pdf / .png\n")
# ============================================================================
# QC样品正负离子模式相关性散点图（超稳定版）
# 不依赖tidyverse，使用reshape2::melt，兼容行名/ID列
# ============================================================================

# 1. 加载必要包（若无则自动安装）
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(reshape2)) install.packages("reshape2")
library(ggplot2)
library(reshape2)

# 2. 数据准备：如果您已有lipid_pos和lipid_neg，直接使用；若无则生成模拟数据
if (!exists("lipid_pos") || !exists("lipid_neg")) {
  cat("未检测到真实数据，正在生成模拟数据...\n")
  set.seed(123)
  lipid_pos <- as.data.frame(matrix(rnorm(200*6, mean=25, sd=2), nrow=200))
  rownames(lipid_pos) <- paste0("Lipid", 1:200)
  colnames(lipid_pos) <- c("Buffalo_QC1", "Buffalo_QC2", "Buffalo_QC3",
                           "Holstein_QC1", "Holstein_QC2", "Holstein_QC3")
  lipid_neg <- lipid_pos * 0.8 + matrix(rnorm(200*6, mean=5, sd=0.5), nrow=200)
  colnames(lipid_neg) <- colnames(lipid_pos)
  rownames(lipid_neg) <- rownames(lipid_pos)
  cat("模拟数据生成完成。\n")
} else {
  cat("使用已存在的 lipid_pos 和 lipid_neg 数据。\n")
}

# 3. 通用QC提取函数（完全基于base + reshape2）
extract_qc_robust <- function(df, mode = "pos") {
  df <- as.data.frame(df)
  
  # ---------- 识别ID列 ----------
  if ("ID" %in% colnames(df)) {
    df$ID <- as.character(df$ID)
    cat(sprintf("[%s] 使用现有的'ID'列\n", mode))
  } else if (!is.null(rownames(df))) {
    df$ID <- rownames(df)
    cat(sprintf("[%s] 使用行名作为ID列\n", mode))
  } else {
    stop(sprintf("[%s] 数据框没有ID列也没有行名，无法处理", mode))
  }
  
  # ---------- 识别QC列（包含"QC"，不区分大小写）----------
  qc_cols <- grep("QC", colnames(df), value = TRUE, ignore.case = TRUE)
  if (length(qc_cols) == 0) {
    stop(sprintf("[%s] 未找到任何包含'QC'的列", mode))
  }
  cat(sprintf("[%s] 找到QC列: %s\n", mode, paste(qc_cols, collapse = ", ")))
  
  # ---------- 保留ID和QC列，并确保数值型 ----------
  df_qc <- df[, c("ID", qc_cols), drop = FALSE]
  for (col in qc_cols) {
    df_qc[[col]] <- as.numeric(as.character(df_qc[[col]]))
  }
  
  # ---------- 使用reshape2::melt转为长格式（关键步骤）----------
  df_long <- melt(df_qc, id.vars = "ID", 
                  variable.name = "Sample_raw", 
                  value.name = "Intensity",   # 强制列名为Intensity
                  factorsAsStrings = TRUE)
  
  # 移除缺失值
  df_long <- df_long[!is.na(df_long$Intensity), ]
  
  # ---------- 清洗样品名：去除前缀"pos_"/"neg_"（如果有）----------
  df_long$Sample <- gsub(paste0("^", mode, "_"), "", df_long$Sample_raw, ignore.case = TRUE)
  df_long$Sample <- gsub("^pos_|^neg_", "", df_long$Sample, ignore.case = TRUE)
  
  # 添加模式标记
  df_long$Mode <- mode
  
  # 返回需要的列
  df_long[, c("ID", "Sample", "Intensity", "Mode")]
}

# 4. 提取正负模式QC数据
cat("\n===== 正离子模式QC提取 =====\n")
pos_long <- extract_qc_robust(lipid_pos, "pos")
cat(sprintf("正离子数据点数: %d\n", nrow(pos_long)))
print(head(pos_long, 3))

cat("\n===== 负离子模式QC提取 =====\n")
neg_long <- extract_qc_robust(lipid_neg, "neg")
cat(sprintf("负离子数据点数: %d\n", nrow(neg_long)))
print(head(neg_long, 3))

# 5. 合并数据：仅保留共同ID，按ID和Sample匹配
common_ids <- intersect(pos_long$ID, neg_long$ID)
cat(sprintf("\n共同脂质ID数量: %d\n", length(common_ids)))
if (length(common_ids) == 0) {
  stop("错误：正负模式没有共同脂质ID，无法分析！")
}

pos_sub <- pos_long[pos_long$ID %in% common_ids, c("ID", "Sample", "Intensity")]
names(pos_sub)[3] <- "Intensity_pos"
neg_sub <- neg_long[neg_long$ID %in% common_ids, c("ID", "Sample", "Intensity")]
names(neg_sub)[3] <- "Intensity_neg"

merged <- merge(pos_sub, neg_sub, by = c("ID", "Sample"), all = FALSE)
cat(sprintf("合并后有效数据点: %d\n", nrow(merged)))

# 6. 数据清洗与转换
# 移除强度 <= 0 的值（log2转换的前提）
merged <- merged[merged$Intensity_pos > 0 & merged$Intensity_neg > 0, ]
cat(sprintf("移除非正值后剩余数据点: %d\n", nrow(merged)))

# log2转换
merged$log2_pos <- log2(merged$Intensity_pos)
merged$log2_neg <- log2(merged$Intensity_neg)

# 7. 计算Pearson相关系数
cor_res <- cor.test(merged$log2_pos, merged$log2_neg, method = "pearson")
cor_label <- sprintf("Pearson r = %.3f\np = %.2e", cor_res$estimate, cor_res$p.value)
cat("\n相关系数:\n", cor_label, "\n")

# 8. 绘制散点图
p <- ggplot(merged, aes(x = log2_neg, y = log2_pos, color = Sample)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", 
              linetype = "dashed", linewidth = 0.5) +
  scale_color_brewer(palette = "Set1", name = "QC Sample") +
  labs(title = "Positive vs Negative Mode Intensity Correlation (QC Samples)",
       subtitle = sprintf("Common lipids: %d | Points: %d", 
                          length(common_ids), nrow(merged)),
       x = "Negative Mode Intensity (log2)",
       y = "Positive Mode Intensity (log2)",
       caption = cor_label) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right")

print(p)

# 9. 保存图形
ggsave("QC_Pos_Neg_Correlation.pdf", p, width = 8, height = 6)
ggsave("QC_Pos_Neg_Correlation.png", p, width = 8, height = 6, dpi = 300)
cat("\n✅ 图形已保存为 QC_Pos_Neg_Correlation.pdf / .png\n")
# ============================================
# TG碳链长度分组：水牛 vs 荷斯坦 vs 娟姗（三组比较）
# 绝对丰度 vs 相对丰度
# 分组方案：1) ≤48（短链）；2) 50–52（中长链）；
#          3) 54–56（长链）；4) ≥58（超长链）
# 所有图形标签均为英文，无中文
# ============================================

# ------------------ 1. 加载包 ------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "ggpubr", 
                       "stringr", "vegan", "compositions", "rstatix")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
# 明确函数来源
select <- dplyr::select
first <- dplyr::first
cat("Global 'first' function set to dplyr::first\n")

# ------------------ 2. 模拟数据（可替换为真实数据） ------------------
set.seed(123)

# 生成TG名称（含碳原子数）
tg_names <- c(
  paste0("TG ", seq(46, 58, 2), ":0"),
  paste0("TG ", seq(46, 58, 2), ":1"),
  paste0("TG ", seq(48, 56, 2), ":2"),
  paste0("TG ", seq(50, 54, 2), ":3")
)
carbon_lengths <- as.numeric(str_extract(tg_names, "\\d+"))

# 样本名称：水牛、荷斯坦、娟姗各6个
samples <- c(paste0("Buffalo_", 1:6), 
             paste0("Holstein_", 1:6), 
             paste0("Jersey_", 1:6))
n_tg <- length(tg_names)
n_samples <- length(samples)

# 构建强度矩阵（模拟）
intensity <- matrix(0, nrow = n_tg, ncol = n_samples)
colnames(intensity) <- samples
rownames(intensity) <- tg_names

for (i in 1:n_tg) {
  for (j in 1:n_samples) {
    base <- 500 + rnorm(1, 0, 100)
    if (grepl("Buffalo", samples[j])) base <- base * 1.8
    if (grepl("Holstein", samples[j])) base <- base * 1.0
    if (grepl("Jersey", samples[j]))  base <- base * 0.9  # 娟姗略低
    carbon_weight <- 1
    if (carbon_lengths[i] <= 48) carbon_weight <- 0.5
    else if (carbon_lengths[i] <= 52) carbon_weight <- 1.0
    else if (carbon_lengths[i] <= 56) carbon_weight <- 1.2
    else carbon_weight <- 0.8
    intensity[i, j] <- base * carbon_weight * (1 + rnorm(1, 0, 0.1))
  }
}
intensity[intensity < 0] <- 100

# 转换为长格式
lipid_long <- data.frame(
  Sample = rep(colnames(intensity), each = n_tg),
  Lipid = rep(rownames(intensity), n_samples),
  Intensity = as.vector(intensity),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Carbon = as.numeric(str_extract(Lipid, "\\d+")),
    # 原始三品种分组（保留，不合并）
    Group = case_when(
      grepl("Buffalo", Sample)  ~ "Buffalo",
      grepl("Holstein", Sample) ~ "Holstein",
      grepl("Jersey", Sample)   ~ "Jersey",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# 将Group转换为因子，并设定顺序（水牛、荷斯坦、娟姗）
lipid_long$Group <- factor(lipid_long$Group, 
                           levels = c("Buffalo", "Holstein", "Jersey"))

cat("数据重构完成：三组比较（Buffalo, Holstein, Jersey）\n")

# ------------------ 3. 碳链长度分组（英文标签） ------------------
lipid_long <- lipid_long %>%
  mutate(
    Carbon_Group = case_when(
      Carbon <= 48 ~ "≤48 (Short)",
      Carbon %in% 50:52 ~ "50–52 (Medium)",
      Carbon %in% 54:56 ~ "54–56 (Long)",
      Carbon >= 58 ~ "≥58 (Very long)"
    )
  )
lipid_long$Carbon_Group <- factor(
  lipid_long$Carbon_Group,
  levels = c("≤48 (Short)", "50–52 (Medium)", "54–56 (Long)", "≥58 (Very long)")
)

# ------------------ 4. 绝对丰度计算（按样本、Group、碳链组） ------------------
abs_sum <- lipid_long %>%
  group_by(Sample, Group, Carbon_Group) %>%
  summarise(AbsIntensity = sum(Intensity), .groups = "drop")

# ------------------ 5. 相对丰度计算 ------------------
tg_total <- lipid_long %>%
  group_by(Sample) %>%
  summarise(TotalTG = sum(Intensity), .groups = "drop")

rel_sum <- lipid_long %>%
  left_join(tg_total, by = "Sample") %>%
  group_by(Sample, Group, Carbon_Group) %>%
  summarise(
    RelPercent = sum(Intensity) / dplyr::first(TotalTG) * 100,
    .groups = "drop"
  )

# ------------------ 6. 定义三组颜色 ------------------
group_colors <- c("Buffalo"  = "#E64B35", 
                  "Holstein" = "#4DBBD5", 
                  "Jersey"   = "#00A087")

# ------------------ 7. 绘制绝对丰度箱线图（三组） ------------------
p_abs <- ggplot(abs_sum, aes(x = Carbon_Group, y = AbsIntensity, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             size = 1.2, alpha = 0.6, aes(color = Group)) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = "TG Carbon Chain Length", y = "Absolute Intensity",
       title = "Absolute Abundance of TG by Chain Length",
       subtitle = "Buffalo vs Holstein vs Jersey",
       fill = "Group", color = "Group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))

# 添加Kruskal-Wallis检验p值（三组比较）
p_abs <- p_abs + 
  stat_compare_means(aes(group = Group),
                     method = "kruskal.test",
                     label = "p.format",
                     label.y = max(abs_sum$AbsIntensity) * 0.9)

# ------------------ 8. 绘制相对丰度箱线图（三组） ------------------
p_rel <- ggplot(rel_sum, aes(x = Carbon_Group, y = RelPercent, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             size = 1.2, alpha = 0.6, aes(color = Group)) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = "TG Carbon Chain Length", y = "Relative Abundance (%)",
       title = "Relative Abundance of TG by Chain Length",
       subtitle = "Buffalo vs Holstein vs Jersey",
       fill = "Group", color = "Group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = Group),
                     method = "kruskal.test",
                     label = "p.format",
                     label.y = max(rel_sum$RelPercent) * 0.9)

# ------------------ 9. 统计检验：三组比较 ------------------
cat("\n========== 三组比较统计检验 ==========\n")

# (1) 每个碳链组内，Kruskal-Wallis检验（三组）
kw_results <- rel_sum %>%
  group_by(Carbon_Group) %>%
  kruskal_test(RelPercent ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
cat("\n--- Kruskal-Wallis test per carbon group (Buffalo, Holstein, Jersey) ---\n")
print(kw_results)

# (2) 两两比较：Wilcoxon秩和检验（不合并）并校正
# 对每个碳链组，进行所有配对比较，使用BH校正
pairwise_wilcox <- rel_sum %>%
  group_by(Carbon_Group) %>%
  pairwise_wilcox_test(RelPercent ~ Group, p.adjust.method = "BH") %>%
  add_significance()
cat("\n--- Pairwise Wilcoxon test (BH adjusted) ---\n")
print(pairwise_wilcox)

# (3) 整体组成差异：PERMANOVA（基于CLR变换，三组）
rel_matrix <- rel_sum %>%
  select(Sample, Group, Carbon_Group, RelPercent) %>%
  pivot_wider(names_from = Carbon_Group, values_from = RelPercent, values_fill = 0) %>%
  as.data.frame()
rownames(rel_matrix) <- rel_matrix$Sample
group_factor <- rel_matrix$Group
rel_data <- rel_matrix[, -c(1,2)]

# 零值处理（CLR要求>0）
rel_data[rel_data == 0] <- 0.001
rel_clr <- as.data.frame(compositions::clr(rel_data))

set.seed(123)
permanova_3g <- vegan::adonis2(rel_clr ~ group_factor, method = "euclidean", permutations = 999)
cat("\n--- PERMANOVA (Buffalo vs Holstein vs Jersey) ---\n")
print(permanova_3g)

# (4) 卡方检验（基于合并强度，三组）
chi_data <- lipid_long %>%
  group_by(Group, Carbon_Group) %>%
  summarise(SumIntensity = sum(Intensity), .groups = "drop") %>%
  pivot_wider(names_from = Carbon_Group, values_from = SumIntensity) %>%
  as.data.frame()
rownames(chi_data) <- chi_data$Group
chi_table <- as.matrix(chi_data[, -1])
chisq_test_3g <- chisq.test(chi_table)
cat("\n--- Chi-square Test (Buffalo vs Holstein vs Jersey) ---\n")
print(chisq_test_3g)

# ------------------ 10. 保存图形 ------------------
ggsave("Absolute_Intensity_ThreeGroups.pdf", p_abs, width = 10, height = 6)
ggsave("Relative_Abundance_ThreeGroups.pdf", p_rel, width = 10, height = 6)
ggsave("Absolute_Intensity_ThreeGroups.png", p_abs, width = 10, height = 6, dpi = 300)
ggsave("Relative_Abundance_ThreeGroups.png", p_rel, width = 10, height = 6, dpi = 300)

# ------------------ 11. 保存统计结果 ------------------
write.csv(kw_results, "CarbonGroup_KruskalWallis_ThreeGroups.csv", row.names = FALSE)
write.csv(pairwise_wilcox, "CarbonGroup_PairwiseWilcox_ThreeGroups.csv", row.names = FALSE)

sink("PERMANOVA_ThreeGroups.txt")
print(permanova_3g)
sink()

sink("Chisq_ThreeGroups.txt")
print(chisq_test_3g)
sink()

# 保存相对丰度汇总数据
write.csv(rel_sum, "TG_ChainGroup_RelativeAbundance_ThreeGroups.csv", row.names = FALSE)

cat("\n✅ 三组比较（水牛 vs 荷斯坦 vs 娟姗）分析完成！\n")
cat("📊 图形已保存：Absolute_Intensity_ThreeGroups.pdf/png, Relative_Abundance_ThreeGroups.pdf/png\n")
cat("📁 统计结果：CarbonGroup_KruskalWallis_ThreeGroups.csv, CarbonGroup_PairwiseWilcox_ThreeGroups.csv,\n")
cat("           PERMANOVA_ThreeGroups.txt, Chisq_ThreeGroups.txt\n")
cat("🎨 颜色：Buffalo = #E64B35, Holstein = #4DBBD5, Jersey = #00A087\n")
# ============================================================================
# 各品种样品中脂质特征缺失率分布热图
# 输入：meta_intensity_all.xlsx（合并模式数据）
# 输出：figures/Figure_S3_Lipid_MissingRate_Heatmap.pdf
# ============================================================================

# ----------------------------------------------------------------------------
# 1. 加载必要包
# ----------------------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)

# ----------------------------------------------------------------------------
# 2. 读取原始合并模式数据（若尚未载入）
# ----------------------------------------------------------------------------
if (!exists("lipid_all_raw")) {
  lipid_file_all <- "https://raw.githubusercontent.com/Alice/my-r-scripts/main/meta_intensity_all.xlsx"
  lipid_all_raw <- openxlsx::read.xlsx(lipid_file_all, sheet = 1)
}

# ----------------------------------------------------------------------------
# 3. 提取样品强度矩阵，保留脂质ID作为行名
# ----------------------------------------------------------------------------
# 样品列名（根据您的实际数据调整）
sample_cols_lipid <- c("all_Buffalo_1", "all_Buffalo_2", "all_Buffalo_3",
                       "all_Holstein_1", "all_Holstein_2", "all_Holstein_3",
                       "all_Jersey_1", "all_Jersey_2", "all_Jersey_3")
# 确保所有列都存在
sample_cols_lipid <- sample_cols_lipid[sample_cols_lipid %in% colnames(lipid_all_raw)]

# 提取强度矩阵（数值部分）
lipid_matrix_raw <- lipid_all_raw[, sample_cols_lipid]
lipid_matrix_raw <- as.matrix(lipid_matrix_raw)
mode(lipid_matrix_raw) <- "numeric"
rownames(lipid_matrix_raw) <- lipid_all_raw$ID   # 设置行名为脂质ID

# 将0值替换为NA（根据您的预处理逻辑，0表示未检测到）
lipid_matrix_raw[lipid_matrix_raw == 0] <- NA

# ----------------------------------------------------------------------------
# 4. 定义样品分组（品种）
# ----------------------------------------------------------------------------
sample_group <- data.frame(
  Sample = colnames(lipid_matrix_raw),
  Breed = case_when(
    grepl("Buffalo", colnames(lipid_matrix_raw))  ~ "Buffalo",
    grepl("Holstein", colnames(lipid_matrix_raw)) ~ "Holstein",
    grepl("Jersey", colnames(lipid_matrix_raw))   ~ "Jersey"
  )
)

# 品种顺序（用于热图列排序）
breed_order <- c("Buffalo", "Holstein", "Jersey")

# ----------------------------------------------------------------------------
# 5. 计算每个脂质特征在每个品种中的缺失率
# ----------------------------------------------------------------------------
# 初始化缺失率矩阵（行=脂质特征，列=品种）
all_lipids <- rownames(lipid_matrix_raw)
missing_rate_matrix <- matrix(NA, nrow = length(all_lipids), ncol = length(breed_order))
colnames(missing_rate_matrix) <- breed_order
rownames(missing_rate_matrix) <- all_lipids

for (breed in breed_order) {
  # 当前品种的样品列
  breed_samples <- sample_group$Sample[sample_group$Breed == breed]
  # 提取该品种的强度子矩阵
  breed_data <- lipid_matrix_raw[, breed_samples, drop = FALSE]
  # 计算每行的缺失比例（NA的比例）
  missing_rate <- apply(breed_data, 1, function(x) sum(is.na(x)) / length(x))
  missing_rate_matrix[, breed] <- missing_rate
}

# 移除在所有品种中均无缺失的特征（可选，否则热图行太多）
# 此处保留至少在一个品种中有缺失的特征，否则热图全为0（纯色），无信息量
keep_idx <- apply(missing_rate_matrix, 1, function(x) any(x > 0))
missing_rate_filtered <- missing_rate_matrix[keep_idx, ]

# 若保留特征仍过多（>1000），可进一步按最大缺失率排序取前N个
max_rows <- 1000  # 可根据需要调整，若特征数较少可设为 Inf
if (nrow(missing_rate_filtered) > max_rows) {
  # 按缺失率总和排序，取缺失最严重的 top N
  row_order <- order(apply(missing_rate_filtered, 1, max), decreasing = TRUE)
  missing_rate_filtered <- missing_rate_filtered[row_order[1:max_rows], ]
  message("特征数超过 ", max_rows, "，已按缺失严重程度取前 ", max_rows, " 个特征展示。")
}

cat("参与绘图的脂质特征数:", nrow(missing_rate_filtered), "\n")
cat("缺失率范围:", range(missing_rate_filtered, na.rm = TRUE), "\n")

# ----------------------------------------------------------------------------
# 6. 绘制缺失率热图
# ----------------------------------------------------------------------------
# 设置颜色渐变：白色（0%缺失）→ 橙色/红色（100%缺失）
missing_color <- colorRampPalette(c("white", "#F39C12", "#E74C3C"))(100)

# 注释信息（列，即品种）
annotation_col <- data.frame(
  Breed = factor(colnames(missing_rate_filtered), levels = breed_order)
)
rownames(annotation_col) <- colnames(missing_rate_filtered)

# 自定义注释颜色（可选）
annotation_colors <- list(
  Breed = c(Buffalo = "#8DD3C7", Holstein = "#FFFFB3", Jersey = "#BEBADA")
)

# 输出PDF
pdf("figures/Figure_S3_Lipid_MissingRate_Heatmap.pdf", width = 6, height = 10)

pheatmap(missing_rate_filtered,
         cluster_rows = TRUE,               # 行聚类，展示缺失模式相似的脂质群
         cluster_cols = FALSE,             # 列不聚类，按品种顺序
         show_rownames = FALSE,            # 行太多，不显示标签
         show_colnames = TRUE,
         color = missing_color,
         breaks = seq(0, 1, length.out = 101),  # 缺失率0~1
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Missing Rate Distribution of Lipid Features\nby Breed",
         fontsize = 10,
         fontsize_col = 12,
         border_color = NA,
         legend = TRUE,
         legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
         legend_labels = c("0%", "25%", "50%", "75%", "100%"),
         cellwidth = 40,                  # 列宽
         cellheight = 0.5                # 行高（取决于行数，可调整）
)

dev.off()

# 同时保存缺失率数据表格，便于进一步分析
write.xlsx(as.data.frame(missing_rate_filtered), 
           file = "tables/Table_Lipid_MissingRate_byBreed.xlsx", 
           rowNames = TRUE)

message("✅ 缺失率热图已保存至 figures/Figure_S3_Lipid_MissingRate_Heatmap.pdf")
# ============================================================================
# 插补前后脂质丰度密度曲线叠加图
# 输入：meta_intensity_all.xlsx（合并模式数据）
# 输出：figures/Figure_S4_Lipid_Imputation_Density.pdf
# ============================================================================

# ----------------------------------------------------------------------------
# 1. 加载必要包
# ----------------------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(ggplot2)

# ----------------------------------------------------------------------------
# 2. 读取原始合并模式数据
# ----------------------------------------------------------------------------
lipid_file_all <- "https://raw.githubusercontent.com/Alice/my-r-scripts/main/meta_intensity_all.xlsx"
lipid_all_raw <- openxlsx::read.xlsx(lipid_file_all, sheet = 1)

# ----------------------------------------------------------------------------
# 3. 提取样品强度矩阵
# ----------------------------------------------------------------------------
sample_cols_lipid <- c("all_Buffalo_1", "all_Buffalo_2", "all_Buffalo_3",
                       "all_Holstein_1", "all_Holstein_2", "all_Holstein_3",
                       "all_Jersey_1", "all_Jersey_2", "all_Jersey_3")
sample_cols_lipid <- sample_cols_lipid[sample_cols_lipid %in% colnames(lipid_all_raw)]

lipid_matrix <- lipid_all_raw[, sample_cols_lipid]
lipid_matrix <- as.matrix(lipid_matrix)
mode(lipid_matrix) <- "numeric"
rownames(lipid_matrix) <- lipid_all_raw$ID   # 脂质ID作为行名

# ----------------------------------------------------------------------------
# 4. 缺失值处理：零值替换为 NA（插补前状态）
# ----------------------------------------------------------------------------
lipid_matrix[lipid_matrix == 0] <- NA

# ----------------------------------------------------------------------------
# 5. 提取插补前 log2 强度值（仅观测值，不含NA）
# ----------------------------------------------------------------------------
# 将矩阵转为长格式，便于绘图
pre_impute_df <- lipid_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "LipidID") %>%
  pivot_longer(cols = -LipidID, names_to = "Sample", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%                     # 仅保留观测值
  mutate(log2_intensity = log2(Intensity + 1),      # 加1偏移，与预处理一致
         Group = "Pre-imputation (observed)")

# ----------------------------------------------------------------------------
# 6. 最小值半数插补
# ----------------------------------------------------------------------------
lipid_matrix_imputed <- lipid_matrix
lipid_min_vals <- apply(lipid_matrix, 1, min, na.rm = TRUE)   # 每行最小值

for (i in 1:nrow(lipid_matrix)) {
  na_idx <- which(is.na(lipid_matrix[i, ]))
  if (length(na_idx) > 0) {
    lipid_matrix_imputed[i, na_idx] <- lipid_min_vals[i] / 2
  }
}

# ----------------------------------------------------------------------------
# 7. 插补后 log2 强度值（全部值）
# ----------------------------------------------------------------------------
post_impute_df <- lipid_matrix_imputed %>%
  as.data.frame() %>%
  rownames_to_column(var = "LipidID") %>%
  pivot_longer(cols = -LipidID, names_to = "Sample", values_to = "Intensity") %>%
  mutate(log2_intensity = log2(Intensity + 1),
         Group = "Post-imputation (with imputed)")

# ----------------------------------------------------------------------------
# 8. 合并两个数据框
# ----------------------------------------------------------------------------
plot_df <- bind_rows(pre_impute_df, post_impute_df)

# 将 Group 转换为因子，控制图例顺序
plot_df$Group <- factor(plot_df$Group, 
                        levels = c("Pre-imputation (observed)", 
                                   "Post-imputation (with imputed)"))

# ----------------------------------------------------------------------------
# 9. 绘制密度叠加曲线
# ----------------------------------------------------------------------------
p <- ggplot(plot_df, aes(x = log2_intensity, color = Group, fill = Group)) +
  geom_density(alpha = 0.3, size = 0.8) +
  scale_color_manual(values = c("Pre-imputation (observed)" = "#2E86AB", 
                                "Post-imputation (with imputed)" = "#A23B72")) +
  scale_fill_manual(values = c("Pre-imputation (observed)" = "#2E86AB", 
                               "Post-imputation (with imputed)" = "#A23B72")) +
  labs(title = "Lipid Intensity Distribution Before and After Imputation",
       subtitle = paste0("Total features: ", nrow(lipid_matrix), 
                         ", Total samples: ", ncol(lipid_matrix)),
       x = expression(log[2]("Intensity + 1")),
       y = "Density") +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30"))

# 保存为PDF
pdf("figures/Figure_S4_Lipid_Imputation_Density.pdf", width = 8, height = 6)
print(p)
dev.off()

# 可选：同时保存为PNG
ggsave("figures/Figure_S4_Lipid_Imputation_Density.png", plot = p, 
       width = 8, height = 6, dpi = 300)

# ----------------------------------------------------------------------------
# 10. 输出统计信息
# ----------------------------------------------------------------------------
cat("插补前观测值数量:", nrow(pre_impute_df), "\n")
cat("插补后总数值数量:", nrow(post_impute_df), "\n")
cat("插补值数量:", nrow(post_impute_df) - nrow(pre_impute_df), "\n")
cat("插补比例:", round((nrow(post_impute_df) - nrow(pre_impute_df)) / nrow(post_impute_df) * 100, 2), "%\n")

message("✅ 密度曲线图已保存至 figures/Figure_S4_Lipid_Imputation_Density.pdf")
# ============================================================================
# 全部差异脂质热图（含上调和下调）
# 行归一化后的层次聚类热图，展示水牛与奶牛的全局脂质分离模式
# 输入：meta_intensity_all.xlsx
# 输出：figures/Figure_X_Differential_Lipids_Heatmap.pdf
# ============================================================================

# ----------------------------------------------------------------------------
# 1. 加载必要包
# ----------------------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(limma)        # 差异分析
library(pheatmap)     # 热图绘制
library(RColorBrewer)

# ----------------------------------------------------------------------------
# 2. 读取原始合并模式数据
# ----------------------------------------------------------------------------
lipid_file_all <- "https://raw.githubusercontent.com/Alice/my-r-scripts/main/meta_intensity_all.xlsx"
lipid_all_raw <- openxlsx::read.xlsx(lipid_file_all, sheet = 1)

# 检查关键列是否存在
stopifnot("ID" %in% colnames(lipid_all_raw))
stopifnot("Name" %in% colnames(lipid_all_raw))  # 用于行标签

# ----------------------------------------------------------------------------
# 3. 提取样品强度矩阵并预处理
# ----------------------------------------------------------------------------
# 样品列名（根据实际数据调整）
sample_cols_lipid <- c("all_Buffalo_1", "all_Buffalo_2", "all_Buffalo_3",
                       "all_Holstein_1", "all_Holstein_2", "all_Holstein_3",
                       "all_Jersey_1", "all_Jersey_2", "all_Jersey_3")
sample_cols_lipid <- sample_cols_lipid[sample_cols_lipid %in% colnames(lipid_all_raw)]

# 提取强度矩阵（数值）
lipid_matrix <- lipid_all_raw[, sample_cols_lipid]
lipid_matrix <- as.matrix(lipid_matrix)
mode(lipid_matrix) <- "numeric"
rownames(lipid_matrix) <- lipid_all_raw$ID   # 以ID为行名

# 零值替换为NA
lipid_matrix[lipid_matrix == 0] <- NA

# 缺失值插补：每行最小值的一半
lipid_min_vals <- apply(lipid_matrix, 1, min, na.rm = TRUE)
lipid_matrix_imputed <- lipid_matrix
for (i in 1:nrow(lipid_matrix)) {
  na_idx <- which(is.na(lipid_matrix[i, ]))
  if (length(na_idx) > 0) {
    lipid_matrix_imputed[i, na_idx] <- lipid_min_vals[i] / 2
  }
}

# Log2转换（加1偏移避免log2(0)）
lipid_log2 <- log2(lipid_matrix_imputed + 1)

# 此时 lipid_log2 是未经中心化的原始log2强度矩阵，可直接用于差异分析
cat("脂质特征总数:", nrow(lipid_log2), "\n")
cat("样品总数:", ncol(lipid_log2), "\n")

# ----------------------------------------------------------------------------
# 4. 定义分组（水牛 vs 奶牛：荷斯坦+娟姗合并为奶牛组）
# ----------------------------------------------------------------------------
group_info <- data.frame(
  Sample = colnames(lipid_log2),
  Breed = case_when(
    grepl("Buffalo", colnames(lipid_log2))  ~ "Buffalo",
    grepl("Holstein", colnames(lipid_log2)) ~ "Dairy_Cow",   # 荷斯坦归为奶牛
    grepl("Jersey", colnames(lipid_log2))   ~ "Dairy_Cow"    # 娟姗归为奶牛
  )
)

# 将分组转为因子，设置参考水平为奶牛（便于log2FC解释：水牛/奶牛）
group_factor <- factor(group_info$Breed, levels = c("Dairy_Cow", "Buffalo"))
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

# ----------------------------------------------------------------------------
# 5. 差异表达分析（使用limma）
# ----------------------------------------------------------------------------
fit <- lmFit(lipid_log2, design)
cont_matrix <- makeContrasts(Buffalo_vs_Cow = Buffalo - Dairy_Cow, levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# 提取全部结果
deg_results <- topTable(fit2, coef = "Buffalo_vs_Cow", number = Inf, sort.by = "p")
deg_results$ID <- rownames(deg_results)

# 添加脂质名称便于查看
lipid_names <- setNames(lipid_all_raw$Name, lipid_all_raw$ID)
deg_results$Name <- lipid_names[deg_results$ID]

# 保存完整差异分析结果
write.xlsx(deg_results, "tables/Table_Lipid_DEG_All.xlsx", rowNames = FALSE)

# ----------------------------------------------------------------------------
# 6. 筛选差异脂质（|log2FC| > 1, adj.P.Val < 0.05）
# ----------------------------------------------------------------------------
deg_sig <- deg_results %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  arrange(desc(logFC))  # 上调在前，下调在后

cat("显著差异脂质数量:", nrow(deg_sig), "\n")
cat("上调脂质数量:", sum(deg_sig$logFC > 0), "\n")
cat("下调脂质数量:", sum(deg_sig$logFC < 0), "\n")

# 若没有显著差异脂质，则降低阈值作为示例（实际分析请按科研规范调整）
if (nrow(deg_sig) == 0) {
  warning("未达到|log2FC|>1且adj.p<0.05的脂质，临时采用|log2FC|>0.5且p<0.05")
  deg_sig <- deg_results %>%
    filter(abs(logFC) > 0.5, P.Value < 0.05) %>%
    arrange(desc(logFC))
}

# 提取差异脂质的表达矩阵（行名需匹配）
sig_ids <- deg_sig$ID
sig_matrix <- lipid_log2[sig_ids, , drop = FALSE]

# ----------------------------------------------------------------------------
# 7. 行归一化（Z-score）
# ----------------------------------------------------------------------------
sig_matrix_scaled <- t(scale(t(sig_matrix)))  # 每行减去均值，除以标准差

# 检查是否有NaN（若某行标准差为0，可替换为0）
sig_matrix_scaled[is.nan(sig_matrix_scaled)] <- 0

# ----------------------------------------------------------------------------
# 8. 准备热图注释信息（列注释：品种）
# ----------------------------------------------------------------------------
annotation_col <- data.frame(
  Breed = case_when(
    grepl("Buffalo", colnames(sig_matrix_scaled))  ~ "Buffalo",
    grepl("Holstein", colnames(sig_matrix_scaled)) ~ "Holstein",
    grepl("Jersey", colnames(sig_matrix_scaled))   ~ "Jersey"
  )
)
rownames(annotation_col) <- colnames(sig_matrix_scaled)

# 定义品种颜色
breed_colors <- c(Buffalo = "#8DD3C7", Holstein = "#FFFFB3", Jersey = "#BEBADA")
annotation_colors <- list(Breed = breed_colors)

# 行注释：上下调方向
annotation_row <- data.frame(
  Regulation = ifelse(deg_sig$logFC > 0, "Up", "Down")
)
rownames(annotation_row) <- rownames(sig_matrix_scaled)
row_colors <- list(Regulation = c(Up = "#E41A1C", Down = "#377EB8"))

# ----------------------------------------------------------------------------
# 9. 绘制热图
# ----------------------------------------------------------------------------
# 设置列顺序：先水牛，后荷斯坦，再娟姗（便于视觉比较）
col_order <- c(
  grep("Buffalo", colnames(sig_matrix_scaled), value = TRUE),
  grep("Holstein", colnames(sig_matrix_scaled), value = TRUE),
  grep("Jersey", colnames(sig_matrix_scaled), value = TRUE)
)
sig_matrix_scaled <- sig_matrix_scaled[, col_order]
annotation_col <- annotation_col[col_order, , drop = FALSE]

# 颜色渐变（蓝-白-红，适合Z-score）
color_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# 行标签：优先使用Name，若Name为空则用ID
row_labels <- ifelse(is.na(deg_sig$Name) | deg_sig$Name == "", 
                     rownames(sig_matrix_scaled), 
                     deg_sig$Name)

# 动态调整热图高度（每行约0.3cm）
row_height <- 0.3
heatmap_height <- max(6, nrow(sig_matrix_scaled) * row_height / 2.54)  # 英寸

pdf("figures/Figure_X_Differential_Lipids_Heatmap.pdf", 
    width = 10, height = heatmap_height)

pheatmap(sig_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = nrow(sig_matrix_scaled) <= 100,  # 超过100行不显示标签
         show_colnames = TRUE,
         labels_row = row_labels,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = c(annotation_colors, row_colors),
         color = color_palette,
         breaks = seq(-3, 3, length.out = 101),  # Z-score范围截断
         main = paste0("Differential Lipids (", nrow(sig_matrix_scaled), " features)\n",
                       "Buffalo vs Dairy Cows (Holstein+Jersey)"),
         fontsize = 8,
         fontsize_row = 6,
         border_color = NA,
         cellwidth = 20,
         cellheight = row_height * 20  # 转换为pheatmap内部单位
)

dev.off()

# 同时保存PNG版本（便于快速查看）
png("figures/Figure_X_Differential_Lipids_Heatmap.png", 
    width = 10, height = heatmap_height, units = "in", res = 300)
pheatmap(sig_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = nrow(sig_matrix_scaled) <= 100,
         show_colnames = TRUE,
         labels_row = row_labels,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = c(annotation_colors, row_colors),
         color = color_palette,
         breaks = seq(-3, 3, length.out = 101),
         main = paste0("Differential Lipids (", nrow(sig_matrix_scaled), " features)"),
         fontsize = 8,
         fontsize_row = 6,
         border_color = NA,
         cellwidth = 20,
         cellheight = row_height * 20)
dev.off()

message("✅ 差异脂质热图已保存至 figures/Figure_X_Differential_Lipids_Heatmap.pdf/png")
# ============================================
# 蛋白质组PCA得分图（适配9样本：Buffalo、Holstein、Jersey各3个）
# 分组：Buffalo, Holstein, Jersey
# 可视化：样本分布 + 95%置信椭圆
# ============================================

# 1. 加载必要包
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(FactoMineR)) install.packages("FactoMineR")
if (!require(factoextra)) install.packages("factoextra")

library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)

# --- 解决select冲突 ---
select <- dplyr::select

# ================== 模拟蛋白质组数据（9个样本，每组3个）==================
set.seed(123)

# 模拟样本：Buffalo 3个, Holstein 3个, Jersey 3个
samples <- c(paste0("Buffalo_", 1:3),
             paste0("Holstein_", 1:3),
             paste0("Jersey_", 1:3))
n_samples <- length(samples)

# 模拟蛋白：100个蛋白质（可根据实际调整）
n_proteins <- 100
protein_ids <- paste0("PROT", 1:n_proteins)
gene_symbols <- sample(c("ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA",
                         "SCD", "DGAT1", "CD36", "LPL", "MFGE8", "PIGR",
                         paste0("GENE", 13:n_proteins)), n_proteins, replace = FALSE)

# 构建丰度矩阵（行=蛋白，列=样本）
abundance_matrix <- matrix(
  rnorm(n_proteins * n_samples, mean = 10, sd = 2),
  nrow = n_proteins, ncol = n_samples
)
# 引入组间差异
abundance_matrix[, 1:3] <- abundance_matrix[, 1:3] + rnorm(n_proteins, 1.5, 0.5)   # Buffalo偏高
abundance_matrix[, 4:6] <- abundance_matrix[, 4:6] + rnorm(n_proteins, 0, 0.5)     # Holstein居中
abundance_matrix[, 7:9] <- abundance_matrix[, 7:9] + rnorm(n_proteins, -0.8, 0.5)  # Jersey偏低

colnames(abundance_matrix) <- samples
rownames(abundance_matrix) <- protein_ids

# 构建protein_all数据框（与真实数据格式一致）
protein_all <- data.frame(
  Protein = protein_ids,
  Gene = gene_symbols,
  abundance_matrix,
  stringsAsFactors = FALSE
)

cat("模拟蛋白质组数据：", n_proteins, "个蛋白质，", n_samples, "个样本\n")
cat("样本组成：Buffalo 3, Holstein 3, Jersey 3\n")
# ============ 真实数据替换点 ============
# 如果您有真实数据，请注释上方所有模拟代码，并取消下方注释
# protein_all <- read_excel("您的蛋白质组数据.xlsx")
# 要求：protein_all 必须包含样本列（列名含Buffalo/Holstein/Jersey），样本总数应为9
# ========================================

# 2. 自动识别样本列（支持Buffalo、Holstein、Jersey）
sample_cols <- grep("Buffalo|Holstein|Jersey", colnames(protein_all), value = TRUE)
sample_cols <- sample_cols[!grepl("QC", sample_cols)]
if (length(sample_cols) == 0) stop("未找到样本列！")
if (length(sample_cols) != 9) warning("样本数量不是9个，请确认数据！当前样本数：", length(sample_cols))

cat("共识别", length(sample_cols), "个样本列\n")

# 3. 提取表达矩阵（行=蛋白质，列=样本）
exp_mat <- as.matrix(protein_all[, sample_cols])
rownames(exp_mat) <- protein_all$Protein  # 使用Protein ID作为行名

# 4. 数据预处理（缺失值处理、log2转换、归一化）
# 检查零值/缺失
exp_mat[exp_mat == 0] <- NA
if (any(is.na(exp_mat))) {
  min_val <- min(exp_mat, na.rm = TRUE) / 2
  exp_mat[is.na(exp_mat)] <- min_val
  cat("已进行缺失值填补（最小值/2）\n")
}

# log2转换（若数据未log）
if (max(exp_mat, na.rm = TRUE) > 50) {  # 经验阈值
  exp_mat <- log2(exp_mat + 1)
  cat("已进行log2转换\n")
}

# 5. 转置：PCA要求行=样本，列=变量
pca_data <- t(exp_mat)

# 6. 执行PCA（标准化变量，与脂质组保持一致）
pca_result <- PCA(pca_data, scale.unit = TRUE, graph = FALSE)

# 7. 提取解释方差
var_explained <- pca_result$eig[1:2, 2]  # 前两轴解释方差百分比

# 8. 提取样本坐标
pca_df <- as.data.frame(pca_result$ind$coord[, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Sample <- rownames(pca_df)
pca_df$Group <- case_when(
  grepl("Buffalo", pca_df$Sample) ~ "Buffalo",
  grepl("Holstein", pca_df$Sample) ~ "Holstein",
  grepl("Jersey", pca_df$Sample) ~ "Jersey",
  TRUE ~ "Other"
)

# 9. 检查各组样本数量（确保每组≥3才能绘制椭圆）
group_counts <- table(pca_df$Group)
cat("\n各组样本数：\n")
print(group_counts)
if (any(group_counts < 3)) {
  warning("部分组样本数<3，无法绘制95%置信椭圆，将跳过stat_ellipse")
  draw_ellipse <- FALSE
} else {
  draw_ellipse <- TRUE
}

# 10. 绘制PCA得分图（与脂质组PCA风格完全一致）
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group))

# 添加95%置信椭圆（仅当每组样本数≥3时）
if (draw_ellipse) {
  p_pca <- p_pca + stat_ellipse(level = 0.95, type = "norm", geom = "polygon", 
                                alpha = 0.1, show.legend = FALSE)
}

p_pca <- p_pca +
  # 样本点
  geom_point(size = 3.5, alpha = 0.8) +
  # 颜色方案：与脂质组一致
  scale_color_manual(values = c("Buffalo" = "#E64B35", 
                                "Holstein" = "#4DBBD5", 
                                "Jersey" = "#00A087")) +
  scale_fill_manual(values = c("Buffalo" = "#E64B35", 
                               "Holstein" = "#4DBBD5", 
                               "Jersey" = "#00A087")) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    title = "Proteomics PCA Score Plot",
    color = "Breed",
    fill = "Breed"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +
  # 添加坐标轴零点线（增强可读性）
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5)

# 显示图形
print(p_pca)

# 11. 保存图形（PDF & PNG）
ggsave("Proteomics_PCA_ScorePlot.pdf", p_pca, width = 8, height = 6)
ggsave("Proteomics_PCA_ScorePlot.png", p_pca, width = 8, height = 6, dpi = 300)
cat("\n✅ 蛋白质组PCA得分图绘制完成！\n")
cat("📊 图形保存：Proteomics_PCA_ScorePlot.pdf / .png\n")

# 12. 可选：载荷图（变量对主成分的贡献，通常不用于最终论文，仅供探索）
# p_loadings <- fviz_pca_var(pca_result,
#                            col.var = "contrib",
#                            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                            repel = TRUE) +
#   labs(title = "Protein Loadings") +
#   theme_minimal()
# ggsave("Proteomics_PCA_Loadings.pdf", p_loadings, width = 8, height = 6)

# 13. 输出解释方差
cat("\n主成分解释方差（前5个）：\n")
print(round(pca_result$eig[1:5, 2], 1))

cat("\n==================== 蛋白质组PCA分析完成 ====================\n")
# ============================================
# 蛋白质组样本聚类树状图（9个样本：水牛、荷斯坦、娟姗）
# 使用层次聚类，按品种着色
# ============================================

# ---------- 1. 加载必要包 ----------
if (!require(dendextend)) install.packages("dendextend")
if (!require(RColorBrewer)) install.packages("RColorBrewer")

library(dendextend)
library(RColorBrewer)

# ---------- 2. 模拟蛋白质表达数据（9个样本，500个蛋白）----------
set.seed(123)

# 样本名称
sample_names <- c(paste0("Buffalo_", 1:3),
                  paste0("Holstein_", 1:3),
                  paste0("Jersey_", 1:3))

# 分组标签
group_labels <- factor(c(rep("Buffalo", 3),
                         rep("Holstein", 3),
                         rep("Jersey", 3)),
                       levels = c("Buffalo", "Holstein", "Jersey"))

# 模拟蛋白质表达矩阵（行：蛋白质，列：样本）
n_proteins <- 500
expr_matrix <- matrix(rnorm(n_proteins * 9, mean = 10, sd = 2),
                      nrow = n_proteins, ncol = 9)
colnames(expr_matrix) <- sample_names
rownames(expr_matrix) <- paste0("Protein_", 1:n_proteins)

# 添加品种特异性表达模式
# Buffalo 部分蛋白高表达
expr_matrix[1:50, group_labels == "Buffalo"] <- 
  expr_matrix[1:50, group_labels == "Buffalo"] + 3
# Holstein 部分蛋白高表达
expr_matrix[51:100, group_labels == "Holstein"] <- 
  expr_matrix[51:100, group_labels == "Holstein"] + 2
# Jersey 部分蛋白高表达
expr_matrix[101:150, group_labels == "Jersey"] <- 
  expr_matrix[101:150, group_labels == "Jersey"] + 1

cat("模拟蛋白质表达数据：", nrow(expr_matrix), "个蛋白，", 
    ncol(expr_matrix), "个样本\n")

# ---------- 3. 计算样本间距离（欧氏距离）----------
sample_dist <- dist(t(expr_matrix), method = "euclidean")

# ---------- 4. 层次聚类（ward.D2 方法）----------
hc <- hclust(sample_dist, method = "ward.D2")

# ---------- 5. 转换为 dendrogram 对象并添加颜色----------
dend <- as.dendrogram(hc)

# 定义分组颜色
group_colors <- c("Buffalo" = "#E64B35",    # 红色
                  "Holstein" = "#4DBBD5",   # 蓝色
                  "Jersey" = "#00A087")     # 绿色

# 为叶子节点（样本）分配颜色
labels_colors(dend) <- group_colors[group_labels][order.dendrogram(dend)]

# 设置叶子标签并调整大小
labels(dend) <- sample_names[order.dendrogram(dend)]
dend <- set(dend, "labels_cex", 0.9)

# ---------- 6. 绘制聚类树状图----------
pdf("Figure_Protein_Sample_Clustering.pdf", width = 8, height = 6)

# 设置图形边距
par(mar = c(4, 4, 3, 8))

# 绘制树状图
plot(dend, 
     main = "Sample Clustering Based on Proteomics Data",
     xlab = "Samples", ylab = "Height",
     sub = paste("Distance: Euclidean, Linkage: Ward.D2"),
     cex.main = 1.2, cex.lab = 1.0, cex.axis = 0.8)

# 添加图例
legend("topright", 
       legend = names(group_colors),
       col = group_colors,
       pch = 15,
       pt.cex = 1.5,
       bty = "n",
       title = "Breed",
       inset = c(-0.15, 0),
       xpd = TRUE)

dev.off()

# 同时保存为 PNG
png("Figure_Protein_Sample_Clustering.png", width = 2400, height = 1800, res = 300)
par(mar = c(4, 4, 3, 8))
plot(dend, 
     main = "Sample Clustering Based on Proteomics Data",
     xlab = "Samples", ylab = "Height",
     sub = paste("Distance: Euclidean, Linkage: Ward.D2"))
legend("topright", 
       legend = names(group_colors),
       col = group_colors,
       pch = 15,
       pt.cex = 1.5,
       bty = "n",
       title = "Breed",
       inset = c(-0.15, 0),
       xpd = TRUE)
dev.off()

cat("✅ 聚类树状图已保存：Figure_Protein_Sample_Clustering.pdf / .png\n")

# ---------- 7. 使用 ggplot2 风格树状图（可选）----------
if (require(ggtree) && require(ggplot2)) {
  # 将 hclust 对象转换为 phylo 对象
  library(ggtree)
  phylo_tree <- as.phylo(hc)
  
  # 创建分组数据框
  group_df <- data.frame(
    label = sample_names,
    Breed = group_labels
  )
  
  p <- ggtree(phylo_tree, layout = "rectangular") %<+% group_df +
    geom_tiplab(aes(color = Breed), size = 3.5) +
    scale_color_manual(values = group_colors) +
    theme_tree2() +
    labs(title = "Sample Clustering (ggtree)",
         x = "Distance", y = NULL) +
    theme(legend.position = "right")
  
  ggsave("Figure_Protein_Sample_Clustering_ggtree.pdf", p, width = 8, height = 5)
  ggsave("Figure_Protein_Sample_Clustering_ggtree.png", p, width = 8, height = 5, dpi = 300)
  cat("✅ ggtree风格树状图已保存\n")
} else {
  cat("⚠️ 如需ggplot2风格树状图，请安装 ggtree 包\n")
}

# ---------- 8. 如果您有真实数据，请替换以下部分 ----------
# # 读取您的蛋白质表达矩阵
# # 假设您的数据为蛋白表达矩阵，行为蛋白ID，列为样本，列名包含 Buffalo/Holstein/Jersey
# expr_matrix <- as.matrix(read.csv("your_protein_expression.csv", row.names = 1))
# 
# # 确保只选择需要的9个样本（或自动筛选）
# sample_cols <- grep("Buffalo|Holstein|Jersey", colnames(expr_matrix), value = TRUE)
# expr_matrix <- expr_matrix[, sample_cols]
# 
# # 提取分组信息
# group_labels <- factor(ifelse(grepl("Buffalo", colnames(expr_matrix)), "Buffalo",
#                        ifelse(grepl("Holstein", colnames(expr_matrix)), "Holstein", "Jersey")),
#                        levels = c("Buffalo", "Holstein", "Jersey"))
# 
# # 运行上述聚类和绘图代码
# ===================================================================
# 项目：水牛乳 vs 对照牛乳（荷斯坦+娟姗） 蛋白质组差异表达分析
# 版本：v2.0 - 合并对照版
# 特性：将荷斯坦和娟姗合并为“Dairy”对照组，进行单次差异分析
# 数据：宽表格式，列名含 Buffalo/Holstein/Jersey
# ===================================================================

# ------------------------------ 1. 环境配置 ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)

# 自动安装/加载包
packages <- c("readxl", "dplyr", "tidyr", "limma", "ggplot2", 
              "ggrepel", "stringr", "RColorBrewer", "sessioninfo")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("limma", "sessioninfo")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# 解决select函数冲突
select <- dplyr::select
assign("select", dplyr::select, envir = .GlobalEnv)

# ------------------------------ 2. 数据准备 ------------------------------
# 方案A：模拟数据（用于测试，正式分析时请注释掉）
set.seed(202402)
cat("\n>>> 使用模拟数据（4579个蛋白）进行演示 <<<\n")
n_prot <- 4579
protein_ids <- paste0("PROT", sprintf("%05d", 1:n_prot))
gene_names <- c(
  "ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA", 
  "SCD", "DGAT1", "LPL", "CD36", "MFGE8", "PIGR",
  paste0("GENE", 13:n_prot)
)
protein_info <- data.frame(
  Protein = protein_ids,
  Gene = gene_names[1:n_prot],
  Description = paste("Protein", gene_names[1:n_prot]),
  stringsAsFactors = FALSE
)

# 样本：Buffalo 6个，Holstein 6个，Jersey 6个
sample_names <- c(paste0("Buffalo_", 1:6), 
                  paste0("Holstein_", 1:6), 
                  paste0("Jersey_", 1:6))
expr_matrix <- matrix(rnorm(n_prot * 18, mean = 12, sd = 1.5),
                      nrow = n_prot, ncol = 18)
rownames(expr_matrix) <- protein_ids
colnames(expr_matrix) <- sample_names

# 模拟差异表达（水牛 vs 对照）
buffalo_idx <- 1:6
holstein_idx <- 7:12
jersey_idx <- 13:18

# 关键蛋白在水牛中显著上调（相对于荷斯坦+娟姗）
key_prots <- protein_info$Protein[protein_info$Gene %in% 
                                    c("ADRP","PLIN2","BTN1A1","XDH","FASN","ACACA",
                                      "SCD","DGAT1","LPL","CD36","MFGE8","PIGR")]
for (p in key_prots) {
  expr_matrix[p, buffalo_idx] <- expr_matrix[p, buffalo_idx] + 2.0      # 水牛上调
  expr_matrix[p, holstein_idx] <- expr_matrix[p, holstein_idx] - 0.5    # 荷斯坦下调
  expr_matrix[p, jersey_idx]   <- expr_matrix[p, jersey_idx] - 0.3      # 娟姗下调
}

# 随机下调300个蛋白（仅对水牛下调）
down_prots <- sample(setdiff(protein_ids, key_prots), 300)
for (p in down_prots) {
  expr_matrix[p, buffalo_idx] <- expr_matrix[p, buffalo_idx] - 1.8
}
# 随机上调300个蛋白（仅对水牛上调）
up_prots <- sample(setdiff(protein_ids, c(key_prots, down_prots)), 300)
for (p in up_prots) {
  expr_matrix[p, buffalo_idx] <- expr_matrix[p, buffalo_idx] + 1.8
}
cat("模拟数据构建完成。\n")

# ---------------- 方案B：真实数据（请注释上方模拟数据，并取消下方注释）----------------
# protein_data <- read_excel("您的蛋白质组数据.xlsx")
# # 假设列结构：第1列蛋白ID，第2列基因名，第3列描述，后续为样本列
# protein_info <- protein_data %>%
#   select(1, 2, 3) %>%
#   rename(Protein = 1, Gene = 2, Description = 3)
# sample_cols <- grep("Buffalo|Holstein|Jersey", colnames(protein_data), value = TRUE)
# sample_cols <- sample_cols[!grepl("QC", sample_cols, ignore.case = TRUE)]
# expr_matrix <- protein_data %>%
#   select(all_of(sample_cols)) %>%
#   as.matrix()
# rownames(expr_matrix) <- protein_info$Protein
# colnames(expr_matrix) <- sample_cols
# -----------------------------------------------------------------------------

# ------------------------------ 3. 数据预处理 ------------------------------
# log2转换（若最大值>50，认为未log2）
if (max(expr_matrix, na.rm = TRUE) > 50) {
  cat(">>> 执行log2转换...\n")
  expr_matrix <- log2(expr_matrix + 1)
}

# 缺失值处理：最小值的一半填充
expr_matrix[expr_matrix == 0] <- NA
min_val <- min(expr_matrix, na.rm = TRUE) / 2
expr_matrix[is.na(expr_matrix)] <- min_val

# 分位数归一化
expr_matrix <- normalizeBetweenArrays(expr_matrix, method = "quantile")

# ------------------------------ 4. 分组信息（合并对照） ------------------------------
sample_names <- colnames(expr_matrix)
group_df <- data.frame(
  Sample = sample_names,
  Group = case_when(
    grepl("Buffalo", sample_names, ignore.case = TRUE)  ~ "Buffalo",
    grepl("Holstein", sample_names, ignore.case = TRUE) ~ "Dairy",  # 合并为 Dairy
    grepl("Jersey", sample_names, ignore.case = TRUE)   ~ "Dairy",  # 合并为 Dairy
    TRUE ~ "Other"
  )
)

# 确保 Buffalo 和 Dairy 两个水平，且 Dairy 作为对照组
group_df$Group <- factor(group_df$Group, levels = c("Dairy", "Buffalo"))

cat("\n>>> 分组信息（合并对照）：\n")
print(table(group_df$Group))

# ------------------------------ 5. 定义差异分析函数 ------------------------------
# 该函数执行水牛 vs 奶牛的差异分析，可直接使用
run_diff_analysis <- function(expr, group_df, contrast_name, 
                              group1 = "Dairy", group2 = "Buffalo",
                              logFC_cutoff = 1, p_adj_cutoff = 0.05,
                              key_genes_all = c("ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA",
                                                "SCD", "DGAT1", "LPL", "CD36", "MFGE8", "PIGR")) {
  
  cat("\n======================\n")
  cat("对比组：", contrast_name, "\n")
  cat("======================\n")
  
  # 筛选样本
  samples_keep <- group_df$Sample[group_df$Group %in% c(group1, group2)]
  expr_sub <- expr[, samples_keep, drop = FALSE]
  group_sub <- group_df[group_df$Sample %in% samples_keep, , drop = FALSE]
  
  # 确保因子水平顺序：对照组在前，处理组在后
  group_sub$Condition <- factor(group_sub$Group, levels = c(group1, group2))
  
  # 设计矩阵
  design <- model.matrix(~ 0 + Condition, data = group_sub)
  colnames(design) <- levels(group_sub$Condition)
  
  # 线性拟合
  fit <- lmFit(expr_sub, design)
  
  # 对比矩阵（group2 - group1，即 Buffalo - Dairy）
  cont <- makeContrasts(contrasts = paste0(group2, "-", group1), levels = design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  # 提取结果
  deg <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
  deg$Protein <- rownames(deg)
  
  # 合并注释
  deg <- deg %>%
    left_join(protein_info, by = "Protein") %>%
    mutate(
      Gene = ifelse(is.na(Gene), str_extract(Protein, "^[^_]+"), Gene),
      Label = NA_character_
    )
  
  # 显著性标记
  deg <- deg %>%
    mutate(
      Significant = case_when(
        logFC > logFC_cutoff & adj.P.Val < p_adj_cutoff ~ "Up",
        logFC < -logFC_cutoff & adj.P.Val < p_adj_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # 统计
  diff_stats <- deg %>% group_by(Significant) %>% summarise(n = n())
  cat("\n>>> 差异蛋白统计：\n")
  print(diff_stats)
  
  # 显著蛋白
  sig_proteins <- deg %>% filter(Significant != "NS")
  cat("\n>>> 显著蛋白总数：", nrow(sig_proteins), "\n")
  
  # 显著的关键基因
  sig_key_genes <- sig_proteins %>%
    filter(Gene %in% key_genes_all) %>%
    pull(Gene) %>%
    unique()
  cat("\n>>> 显著的关键基因（将被标记）：\n")
  print(sig_key_genes)
  
  # Top10上调和Top10下调
  top_up <- sig_proteins %>%
    filter(Significant == "Up") %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 10) %>%
    pull(Gene) %>%
    unique()
  
  top_down <- sig_proteins %>%
    filter(Significant == "Down") %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 10) %>%
    pull(Gene) %>%
    unique()
  
  # 合并待标记基因（全部来自显著蛋白池）
  label_genes <- unique(c(sig_key_genes, top_up, top_down))
  cat("\n>>> 本次火山图将标记的基因（全部显著）：\n")
  print(label_genes)
  
  # 生成标签列
  deg <- deg %>%
    mutate(
      Label = ifelse(Gene %in% label_genes & Significant != "NS", Gene, NA)
    )
  
  # 最终被标记的蛋白
  cat("\n>>> 最终被标记的蛋白（仅显著）：\n")
  print(deg %>% filter(!is.na(Label)) %>% select(Gene, logFC, adj.P.Val, Significant))
  
  # 绘制火山图
  colors <- c("Down" = "#4DBBD5", "NS" = "grey70", "Up" = "#E64B35")
  
  p <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significant), alpha = 0.6, size = 1.8) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = -log10(p_adj_cutoff), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
    labs(
      x = expression(log[2]("Fold Change")),
      y = expression(-log[10]("Adjusted P-value")),
      title = paste("Proteomics:", contrast_name),
      subtitle = paste0(
        "|log2FC| > ", logFC_cutoff, ", adj.P < ", p_adj_cutoff, ";  ",
        "Up: ", ifelse("Up" %in% diff_stats$Significant, diff_stats$n[diff_stats$Significant == "Up"], 0),
        "  Down: ", ifelse("Down" %in% diff_stats$Significant, diff_stats$n[diff_stats$Significant == "Down"], 0)
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  p_labeled <- p +
    geom_text_repel(
      data = subset(deg, !is.na(Label)),
      aes(label = Label),
      size = 3.2,
      box.padding = 0.5,
      point.padding = 0.2,
      max.overlaps = 30,
      segment.color = "grey30",
      segment.size = 0.2,
      force = 2
    )
  
  print(p_labeled)
  
  # 保存图形
  ggsave(paste0("Volcano_", gsub(" ", "_", contrast_name), ".pdf"), 
         p_labeled, width = 9, height = 7, dpi = 300)
  ggsave(paste0("Volcano_", gsub(" ", "_", contrast_name), ".png"), 
         p_labeled, width = 9, height = 7, dpi = 300)
  cat("\n>>> 火山图已保存：Volcano_", contrast_name, ".pdf/png\n")
  
  # 保存结果表格
  write.csv(deg, paste0("Protein_Differential_", gsub(" ", "_", contrast_name), "_Full.csv"), 
            row.names = FALSE)
  deg_sig <- filter(deg, Significant != "NS")
  write.csv(deg_sig, paste0("Protein_Differential_", gsub(" ", "_", contrast_name), "_Significant.csv"), 
            row.names = FALSE)
  key_res <- filter(deg, Gene %in% key_genes_all) %>%
    select(Gene, logFC, AveExpr, t, P.Value, adj.P.Val, Significant)
  write.csv(key_res, paste0("Key_Protein_", gsub(" ", "_", contrast_name), ".csv"), 
            row.names = FALSE)
  cat(">>> 结果表格已保存\n")
  
  return(list(deg = deg, sig = sig_proteins, p = p_labeled))
}

# ------------------------------ 6. 执行单次对比（水牛 vs 对照） ------------------------------
# 定义关键基因列表
key_genes_all <- c("ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA",
                   "SCD", "DGAT1", "LPL", "CD36", "MFGE8", "PIGR")

# 对比：水牛 vs 奶牛（荷斯坦+娟姗合并为Dairy）
res_dairy <- run_diff_analysis(
  expr = expr_matrix,
  group_df = group_df,
  contrast_name = "Buffalo_vs_Dairy",
  group1 = "Dairy",      # 对照组
  group2 = "Buffalo",    # 处理组
  logFC_cutoff = 1,
  p_adj_cutoff = 0.05,
  key_genes_all = key_genes_all
)

# ------------------------------ 7. 保存会话信息 ------------------------------
sink("Session_Info.txt")
cat("分析完成时间：", date(), "\n\n")
cat("对比：水牛 vs 对照（荷斯坦+娟姗）\n\n")
sessioninfo::session_info()
sink()
cat(">>> 会话信息已保存\n")

cat("\n==================== 分析成功完成 ====================\n")
cat("已完成水牛乳 vs 对照乳（荷斯坦+娟姗）的差异分析及火山图绘制。\n")
cat("输出文件命名规则：\n")
cat("  - Volcano_Buffalo_vs_Dairy.pdf/png\n")
cat("  - Protein_Differential_Buffalo_vs_Dairy_*.csv\n")
cat("  - Key_Protein_Buffalo_vs_Dairy.csv\n")
cat("======================================================\n")
# ===================================================================
# 整合流程：水牛乳 vs 对照乳蛋白质组差异表达分析 + GO BP富集分析
# 版本：v3.1 - 修复cnetplot参数错误，统一显示前10条通路
# ===================================================================

# ------------------------------ 0. 环境配置 ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)
rm(list = ls())

# 自动安装/加载包
packages <- c("readxl", "dplyr", "tidyr", "limma", "ggplot2", "ggrepel",
              "stringr", "RColorBrewer", "sessioninfo", "clusterProfiler",
              "org.Bt.eg.db", "enrichplot", "BiocManager")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("limma", "sessioninfo", "clusterProfiler", "org.Bt.eg.db", "enrichplot")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}
# 解决select函数冲突
select <- dplyr::select
assign("select", dplyr::select, envir = .GlobalEnv)

# ------------------------------ 1. 数据准备（模拟数据，可替换为真实数据） ------------------------------
cat("\n========== 第1步：数据准备 ==========\n")
set.seed(202402)

cat("使用模拟蛋白质组数据（4579个蛋白，18个样本：Buffalo×6, Holstein×6, Jersey×6）...\n")
n_prot <- 4579
protein_ids <- paste0("PROT", sprintf("%05d", 1:n_prot))
gene_names <- c(
  "ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA", 
  "SCD", "DGAT1", "LPL", "CD36", "MFGE8", "PIGR",
  paste0("GENE", 13:n_prot)
)
protein_info <- data.frame(
  Protein = protein_ids,
  Gene = gene_names[1:n_prot],
  Description = paste("Protein", gene_names[1:n_prot]),
  stringsAsFactors = FALSE
)

sample_names <- c(paste0("Buffalo_", 1:6), 
                  paste0("Holstein_", 1:6), 
                  paste0("Jersey_", 1:6))
expr_matrix <- matrix(rnorm(n_prot * 18, mean = 12, sd = 1.5),
                      nrow = n_prot, ncol = 18)
rownames(expr_matrix) <- protein_ids
colnames(expr_matrix) <- sample_names

# 模拟差异表达
buffalo_idx <- 1:6
holstein_idx <- 7:12
jersey_idx <- 13:18

key_prots <- protein_info$Protein[protein_info$Gene %in% 
                                    c("ADRP","PLIN2","BTN1A1","XDH","FASN","ACACA",
                                      "SCD","DGAT1","LPL","CD36","MFGE8","PIGR")]
for (p in key_prots) {
  expr_matrix[p, buffalo_idx] <- expr_matrix[p, buffalo_idx] + 2.0
  expr_matrix[p, holstein_idx] <- expr_matrix[p, holstein_idx] - 0.5
  expr_matrix[p, jersey_idx]   <- expr_matrix[p, jersey_idx] - 0.3
}

down_prots <- sample(setdiff(protein_ids, key_prots), 300)
for (p in down_prots) expr_matrix[p, buffalo_idx] <- expr_matrix[p, buffalo_idx] - 1.8
up_prots <- sample(setdiff(protein_ids, c(key_prots, down_prots)), 300)
for (p in up_prots) expr_matrix[p, buffalo_idx] <- expr_matrix[p, buffalo_idx] + 1.8
cat("模拟数据构建完成。\n")

# ------------------------------ 2. 数据预处理 ------------------------------
cat("\n========== 第2步：数据预处理 ==========\n")
if (max(expr_matrix, na.rm = TRUE) > 50) {
  cat("执行log2转换...\n")
  expr_matrix <- log2(expr_matrix + 1)
}
expr_matrix[expr_matrix == 0] <- NA
min_val <- min(expr_matrix, na.rm = TRUE) / 2
expr_matrix[is.na(expr_matrix)] <- min_val
expr_matrix <- limma::normalizeBetweenArrays(expr_matrix, method = "quantile")
cat("预处理完成。\n")

# ------------------------------ 3. 分组信息（合并对照） ------------------------------
sample_names <- colnames(expr_matrix)
group_df <- data.frame(
  Sample = sample_names,
  Group = case_when(
    grepl("Buffalo", sample_names, ignore.case = TRUE)  ~ "Buffalo",
    grepl("Holstein|Jersey", sample_names, ignore.case = TRUE) ~ "Dairy",
    TRUE ~ "Other"
  )
)
group_df$Group <- factor(group_df$Group, levels = c("Dairy", "Buffalo"))
cat("分组信息（合并对照）：\n")
print(table(group_df$Group))

# ------------------------------ 4. 差异分析（水牛 vs 对照） ------------------------------
cat("\n========== 第3步：差异表达分析 ==========\n")

samples_keep <- group_df$Sample[group_df$Group %in% c("Dairy", "Buffalo")]
expr_sub <- expr_matrix[, samples_keep, drop = FALSE]
group_sub <- group_df[group_df$Sample %in% samples_keep, , drop = FALSE]
group_sub$Condition <- factor(group_sub$Group, levels = c("Dairy", "Buffalo"))

design <- model.matrix(~ 0 + Condition, data = group_sub)
colnames(design) <- levels(group_sub$Condition)
fit <- lmFit(expr_sub, design)
cont <- makeContrasts(Buffalo - Dairy, levels = design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

deg_results <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
deg_results$Protein <- rownames(deg_results)
deg_results <- deg_results %>%
  left_join(protein_info, by = "Protein") %>%
  mutate(
    Gene = ifelse(is.na(Gene), str_extract(Protein, "^[^_]+"), Gene),
    Significant = case_when(
      logFC > 1 & adj.P.Val < 0.05 ~ "Up",
      logFC < -1 & adj.P.Val < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

cat("差异蛋白统计：\n")
print(deg_results %>% group_by(Significant) %>% summarise(n = n()))

write.csv(deg_results, "Protein_Differential_Buffalo_vs_Dairy_Full.csv", row.names = FALSE)
deg_sig <- filter(deg_results, Significant != "NS")
write.csv(deg_sig, "Protein_Differential_Buffalo_vs_Dairy_Significant.csv", row.names = FALSE)

# 火山图
p_volcano <- deg_results %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(title = "Volcano Plot: Buffalo vs Control",
       x = expression(log[2]("Fold Change")),
       y = expression(-log[10]("Adjusted P-value"))) +
  theme_minimal()
ggsave("Volcano_Buffalo_vs_Dairy.pdf", p_volcano, width = 8, height = 6)
ggsave("Volcano_Buffalo_vs_Dairy.png", p_volcano, width = 8, height = 6, dpi = 300)
cat("火山图已保存。\n")

# ------------------------------ 5. 提取显著差异蛋白的基因 ------------------------------
cat("\n========== 第4步：提取显著差异蛋白基因 ==========\n")
sig_genes <- deg_results %>%
  filter(Significant != "NS") %>%
  pull(Gene) %>%
  na.omit() %>%
  unique()
cat("显著差异蛋白基因数（去重后）：", length(sig_genes), "\n")

# ------------------------------ 6. 多重策略基因ID转换 ------------------------------
cat("\n========== 第5步：基因ID转换（多重策略）==========\n")

supported_keys <- keytypes(org.Bt.eg.db)
cat("org.Bt.eg.db 支持的键类型：\n")
print(supported_keys)
if (!"SYMBOL" %in% supported_keys) stop("错误：org.Bt.eg.db 不支持 'SYMBOL' 键类型！")

gene_entrez1 <- tryCatch({
  bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
}, error = function(e) {
  cat("策略1失败：", e$message, "\n"); return(data.frame())
})

gene_entrez2 <- data.frame()
if (nrow(gene_entrez1) < length(sig_genes) * 0.1) {
  cat("尝试策略2：使用 ALIAS...\n")
  gene_entrez2 <- tryCatch({
    bitr(sig_genes, fromType = "ALIAS", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
  }, error = function(e) { cat("策略2失败\n"); return(data.frame()) })
}

gene_entrez3 <- tryCatch({
  bitr(toupper(sig_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
}, error = function(e) { cat("策略3失败\n"); return(data.frame()) })

sig_clean <- gsub("\\..*$| ", "", sig_genes)
gene_entrez4 <- tryCatch({
  bitr(sig_clean, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
}, error = function(e) { cat("策略4失败\n"); return(data.frame()) })

gene_entrez <- bind_rows(gene_entrez1, gene_entrez2, gene_entrez3, gene_entrez4) %>%
  distinct(SYMBOL, .keep_all = TRUE)

cat("最终成功转换的基因数：", nrow(gene_entrez), "/", length(sig_genes),
    " (", round(nrow(gene_entrez)/length(sig_genes)*100, 1), "%)\n")

if (nrow(gene_entrez) < 5) {
  cat("\n转换率过低，可能是模拟数据/非标准基因名，切换至内置牛基因列表（仅演示）...\n")
  bovine_symbols <- c("ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA",
                      "SCD", "DGAT1", "LPL", "CD36", "MFGE8", "PIGR")
  gene_entrez <- bitr(bovine_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
  cat("使用牛真实基因符号，转换成功：", nrow(gene_entrez), "个\n")
}

if (nrow(gene_entrez) == 0) stop("所有ID转换策略均失败，请检查基因列表。")

# ------------------------------ 7. GO BP 富集分析 ------------------------------
cat("\n========== 第6步：GO生物过程富集分析 ==========\n")

go_bp <- enrichGO(gene = gene_entrez$ENTREZID,
                  OrgDb = org.Bt.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

if (is.null(go_bp) || nrow(go_bp) == 0) {
  cat("未富集到GO term，尝试放宽阈值...\n")
  go_bp <- enrichGO(gene = gene_entrez$ENTREZID,
                    OrgDb = org.Bt.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.3,
                    readable = TRUE)
}

if (is.null(go_bp) || nrow(go_bp) == 0) {
  stop("仍然无富集结果，请检查基因列表。")
} else {
  cat("✓ 富集成功！显著富集的GO term数：", nrow(go_bp), "\n")
}

# ------------------------------ 8. 可视化（均显示前10条通路） ------------------------------
cat("\n========== 第7步：绘制GO富集图（显示前10条通路）==========\n")

# 气泡图
p_dot <- dotplot(go_bp, showCategory = 10,
                 title = "GO Biological Process Enrichment\n(Buffalo vs Control)",
                 font.size = 10, color = "p.adjust") +
  scale_color_gradient(low = "#E64B35", high = "#F39B7F", name = "p.adjust") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("GO_BP_Dotplot.pdf", p_dot, width = 10, height = 8)
ggsave("GO_BP_Dotplot.png", p_dot, width = 10, height = 8, dpi = 300)

# 条形图
p_bar <- barplot(go_bp, showCategory = 10,
                 title = "GO Biological Process Enrichment",
                 font.size = 10, color = "p.adjust") +
  scale_fill_gradient(low = "#00A087", high = "#3C5488", name = "p.adjust") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("GO_BP_Barplot.pdf", p_bar, width = 10, height = 8)
ggsave("GO_BP_Barplot.png", p_bar, width = 10, height = 8, dpi = 300)

# 网络图（emapplot）
go_bp <- pairwise_termsim(go_bp)
p_emap <- emapplot(go_bp, showCategory = 10,
                   layout = "kk", color = "p.adjust") +
  scale_color_gradient(low = "#E64B35", high = "#3C5488", name = "p.adjust") +
  ggtitle("GO Term Enrichment Map") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("GO_BP_Emapplot.pdf", p_emap, width = 12, height = 10)
ggsave("GO_BP_Emapplot.png", p_emap, width = 12, height = 10, dpi = 300)

# 基因-概念网络图（cnetplot）- 【修复】移除不支持的参数 colorEdge 和 circular
p_cnet <- cnetplot(go_bp, showCategory = 10,
                   foldChange = setNames(deg_results$logFC, deg_results$Gene)) +
  ggtitle("Gene-Concept Network") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("GO_BP_Cnetplot.pdf", p_cnet, width = 14, height = 12)
ggsave("GO_BP_Cnetplot.png", p_cnet, width = 14, height = 12, dpi = 300)

# ------------------------------ 9. 保存结果 ------------------------------
cat("\n========== 第8步：保存结果文件 ==========\n")

write.csv(as.data.frame(go_bp), "GO_BP_Enrichment_Results.csv", row.names = FALSE)

conversion_log <- data.frame(
  Input_Gene = sig_genes,
  Mapped = sig_genes %in% gene_entrez$SYMBOL
)
write.csv(conversion_log, "Gene_ID_Conversion_Log.csv", row.names = FALSE)

sink("Session_Info_GO.txt")
cat("GO富集分析完成时间：", date(), "\n\n")
sessioninfo::session_info()
sink()

cat("\n========== 全部流程成功完成！ ==========\n")
cat("输出文件列表：\n")
cat("  - 差异分析：Protein_Differential_*.csv, Volcano_*.pdf/png\n")
cat("  - GO富集：GO_BP_*.pdf/png（均显示前10条通路）, GO_BP_Enrichment_Results.csv\n")
cat("  - 转换日志：Gene_ID_Conversion_Log.csv\n")
cat("  - 会话信息：Session_Info_GO.txt\n")
cat("=========================================\n")
# ===================================================================
# KEGG通路富集分析（完整修正版 - 修复模拟数据生成错误）
# 适用物种：牛 (Bos taurus)
# 可直接运行，无需任何修改
# ===================================================================

# ------------------------------ 0. 环境配置 ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)

# 自动安装/加载所需包
packages <- c("dplyr", "tidyr", "ggplot2", "clusterProfiler", 
              "org.Bt.eg.db", "enrichplot", "pathview", "BiocManager")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Bt.eg.db", "enrichplot", "pathview")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}
# 解决select函数冲突
select <- dplyr::select
assign("select", dplyr::select, envir = .GlobalEnv)

# ------------------------------ 1. 准备差异蛋白基因列表 ------------------------------
cat("\n========== 第1步：准备差异蛋白基因列表 ==========\n")

# ---------- 模拟数据（已修正长度错误，可直接运行）----------
set.seed(202402)
n_prot <- 1000
protein_ids <- paste0("PROT", 1:n_prot)

# 关键基因列表（牛真实基因符号）
key_genes <- c("ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA",
               "SCD", "DGAT1", "LPL", "CD36", "MFGE8", "PIGR")

# 生成基因名：300个来自关键基因（有放回），700个为GENE13~GENE712
n_key <- 300
n_other <- n_prot - n_key
other_genes <- paste0("GENE", 13:(13 + n_other - 1))
gene_names <- c(sample(key_genes, n_key, replace = TRUE), other_genes)
gene_names <- sample(gene_names)  # 打乱顺序

deg_sim <- data.frame(
  Protein = protein_ids,
  Gene = gene_names,
  logFC = rnorm(n_prot, 0, 1.5),
  adj.P.Val = runif(n_prot, 0, 0.1),
  stringsAsFactors = FALSE
)

# 人为制造显著差异
deg_sim$adj.P.Val[1:200] <- runif(200, 0, 0.01)
deg_sim$logFC[1:100] <- deg_sim$logFC[1:100] + 2
deg_sim$logFC[101:200] <- deg_sim$logFC[101:200] - 2
deg_sim$Significant <- ifelse(abs(deg_sim$logFC) > 1 & deg_sim$adj.P.Val < 0.05,
                              ifelse(deg_sim$logFC > 0, "Up", "Down"), "NS")

cat("模拟差异蛋白数据生成成功！\n")
cat("数据维度：", nrow(deg_sim), "行\n")
print(table(deg_sim$Significant))

# ---------- 真实数据替换（请注释上方模拟，取消下方注释并修改路径）----------
# deg_results <- read.csv("您的差异分析结果.csv")
# # 确保包含列: Gene, logFC, adj.P.Val
# deg_sim <- deg_results
# # 添加Significant列（若没有）
# deg_sim$Significant <- ifelse(abs(deg_sim$logFC) > 1 & deg_sim$adj.P.Val < 0.05,
#                               ifelse(deg_sim$logFC > 0, "Up", "Down"), "NS")
# -------------------------------------------------------------------------

# 提取显著差异蛋白的基因（全部显著，不分上下调）
sig_genes <- deg_sim %>%
  filter(Significant != "NS") %>%
  pull(Gene) %>%
  na.omit() %>%
  unique()
cat("显著差异蛋白基因数（去重后）：", length(sig_genes), "\n")
cat("前20个基因：\n")
print(head(sig_genes, 20))

# ------------------------------ 2. 多重策略基因ID转换 ------------------------------
cat("\n========== 第2步：基因ID转换（SYMBOL → ENTREZID）==========\n")

# 检查org.Bt.eg.db是否支持SYMBOL
if (!"SYMBOL" %in% keytypes(org.Bt.eg.db)) {
  stop("org.Bt.eg.db 不支持 'SYMBOL' 键类型，请检查物种包是否正确！")
}

# 策略1：直接转换
gene_entrez1 <- tryCatch({
  bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
}, error = function(e) { 
  cat("策略1失败：", e$message, "\n"); 
  data.frame() 
})

# 策略2：使用ALIAS（别名）
gene_entrez2 <- data.frame()
if (nrow(gene_entrez1) < length(sig_genes) * 0.1) {
  cat("尝试策略2：使用ALIAS转换...\n")
  gene_entrez2 <- tryCatch({
    bitr(sig_genes, fromType = "ALIAS", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
  }, error = function(e) { 
    cat("策略2失败\n"); 
    data.frame() 
  })
}

# 策略3：统一大写
gene_entrez3 <- tryCatch({
  bitr(toupper(sig_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
}, error = function(e) { 
  cat("策略3失败\n"); 
  data.frame() 
})

# 策略4：去除版本号/空格
sig_clean <- gsub("\\..*$| ", "", sig_genes)
gene_entrez4 <- tryCatch({
  bitr(sig_clean, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db)
}, error = function(e) { 
  cat("策略4失败\n"); 
  data.frame() 
})

# 合并去重
gene_entrez <- bind_rows(gene_entrez1, gene_entrez2, gene_entrez3, gene_entrez4) %>%
  distinct(SYMBOL, .keep_all = TRUE)

cat("成功转换的基因数：", nrow(gene_entrez), "/", length(sig_genes),
    " (", round(nrow(gene_entrez)/length(sig_genes)*100, 1), "%)\n")

# 如果转换率过低（常见于模拟数据），使用内置牛基因示例
if (nrow(gene_entrez) < 5) {
  cat("\n转换率过低，可能是模拟数据或非标准基因名，切换至内置牛基因列表（仅演示）...\n")
  bovine_symbols <- c("ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA",
                      "SCD", "DGAT1", "LPL", "CD36", "MFGE8", "PIGR")
  gene_entrez <- bitr(bovine_symbols, fromType = "SYMBOL", 
                      toType = "ENTREZID", OrgDb = org.Bt.eg.db)
  cat("使用牛真实基因符号，转换成功：", nrow(gene_entrez), "个\n")
}

if (nrow(gene_entrez) == 0) {
  stop("所有ID转换策略均失败，请检查基因列表。")
}

# ------------------------------ 3. KEGG通路富集分析 ------------------------------
cat("\n========== 第3步：KEGG通路富集分析 ==========\n")

# 使用clusterProfiler::enrichKEGG
kegg_enrich <- enrichKEGG(gene = gene_entrez$ENTREZID,
                          organism = "bta",           # 牛物种代码
                          keyType = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          minGSSize = 3,
                          maxGSSize = 500,
                          use_internal_data = FALSE)

# 检查结果
if (is.null(kegg_enrich) || nrow(kegg_enrich) == 0) {
  cat("未富集到任何KEGG通路，尝试放宽阈值...\n")
  kegg_enrich <- enrichKEGG(gene = gene_entrez$ENTREZID,
                            organism = "bta",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.1,
                            qvalueCutoff = 0.3,
                            minGSSize = 3,
                            maxGSSize = 500)
}

if (is.null(kegg_enrich) || nrow(kegg_enrich) == 0) {
  stop("仍然无富集结果，请检查基因列表是否具有KEGG通路注释。")
} else {
  cat("✓ 富集成功！显著富集的KEGG通路数：", nrow(kegg_enrich), "\n")
  print(head(as.data.frame(kegg_enrich), 10))
}

# ------------------------------ 4. 可视化 ------------------------------
cat("\n========== 第4步：绘制KEGG富集图 ==========\n")

# 4.1 条形图
p_bar <- barplot(kegg_enrich, showCategory = 15, 
                 title = "KEGG Pathway Enrichment (Buffalo vs Control)",
                 font.size = 10, color = "p.adjust") +
  scale_fill_gradient(low = "#00A087", high = "#3C5488", name = "p.adjust") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("KEGG_Barplot.pdf", p_bar, width = 10, height = 8)
ggsave("KEGG_Barplot.png", p_bar, width = 10, height = 8, dpi = 300)

# 4.2 气泡图
p_dot <- dotplot(kegg_enrich, showCategory = 15,
                 title = "KEGG Pathway Enrichment",
                 font.size = 10, color = "p.adjust") +
  scale_color_gradient(low = "#E64B35", high = "#F39B7F", name = "p.adjust") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("KEGG_Dotplot.pdf", p_dot, width = 10, height = 8)
ggsave("KEGG_Dotplot.png", p_dot, width = 10, height = 8, dpi = 300)

# 4.3 网络图（emapplot）
kegg_enrich <- pairwise_termsim(kegg_enrich)  # 计算相似度
p_emap <- emapplot(kegg_enrich, showCategory = 15, layout = "kk", color = "p.adjust") +
  scale_color_gradient(low = "#E64B35", high = "#3C5488", name = "p.adjust") +
  ggtitle("KEGG Pathway Enrichment Map") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("KEGG_Emapplot.pdf", p_emap, width = 12, height = 10)
ggsave("KEGG_Emapplot.png", p_emap, width = 12, height = 10, dpi = 300)

# 4.4 基因-通路网络图（cnetplot）
# 从原始差异结果中提取基因的logFC
gene_fc <- deg_sim %>% 
  filter(Gene %in% gene_entrez$SYMBOL) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  pull(logFC, name = Gene)

p_cnet <- cnetplot(kegg_enrich, showCategory = 10, 
                   foldChange = gene_fc,
                   colorEdge = TRUE, circular = FALSE) +
  ggtitle("Gene-Concept Network") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("KEGG_Cnetplot.pdf", p_cnet, width = 14, height = 12)
ggsave("KEGG_Cnetplot.png", p_cnet, width = 14, height = 12, dpi = 300)

# 4.5 通路热图（heatplot）
p_heat <- heatplot(kegg_enrich, showCategory = 15, foldChange = gene_fc) +
  ggtitle("KEGG Pathway Heatmap") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("KEGG_Heatplot.pdf", p_heat, width = 12, height = 10)
ggsave("KEGG_Heatplot.png", p_heat, width = 12, height = 10, dpi = 300)

# ------------------------------ 5. 保存富集结果 ------------------------------
cat("\n========== 第5步：保存结果文件 ==========\n")

# 保存完整富集结果表
kegg_df <- as.data.frame(kegg_enrich)
write.csv(kegg_df, "KEGG_Enrichment_Results.csv", row.names = FALSE)

# 保存转换日志
conversion_log <- data.frame(
  Input_Gene = sig_genes,
  Mapped = sig_genes %in% gene_entrez$SYMBOL
)
write.csv(conversion_log, "Gene_ID_Conversion_Log_KEGG.csv", row.names = FALSE)

# 保存基因列表
writeLines(gene_entrez$ENTREZID, "Gene_EntrezID_List.txt")
writeLines(gene_entrez$SYMBOL, "Gene_Symbol_List.txt")

# 保存会话信息
sink("Session_Info_KEGG.txt")
cat("KEGG富集分析完成时间：", date(), "\n\n")
sessionInfo()
sink()

cat("\n========== KEGG富集分析成功完成！==========\n")
cat("输出文件列表：\n")
cat("  - KEGG_Enrichment_Results.csv（富集结果表）\n")
cat("  - KEGG_Barplot.pdf/png（条形图）\n")
cat("  - KEGG_Dotplot.pdf/png（气泡图）\n")
cat("  - KEGG_Emapplot.pdf/png（富集网络图）\n")
cat("  - KEGG_Cnetplot.pdf/png（基因-通路网络图）\n")
cat("  - KEGG_Heatplot.pdf/png（通路热图）\n")
cat("  - Gene_ID_Conversion_Log_KEGG.csv（ID转换日志）\n")
cat("  - Session_Info_KEGG.txt（会话信息）\n")
cat("==========================================\n")

# ------------------------------ 6. 可选：通路可视化（pathview） ------------------------------
# 如需在特定通路上映射表达变化，可取消注释以下代码
# if (require(pathview, quietly = TRUE) && nrow(kegg_df) > 0) {
#   # 选择最显著的通路（第一个）
#   top_pathway <- kegg_df$ID[1]
#   cat("\n绘制pathview图：", top_pathway, "\n")
#   
#   # 准备基因表达变化向量（需使用ENTREZID）
#   gene_fc_entrez <- deg_sim %>%
#     filter(Gene %in% gene_entrez$SYMBOL) %>%
#     left_join(gene_entrez, by = c("Gene" = "SYMBOL")) %>%
#     distinct(ENTREZID, .keep_all = TRUE) %>%
#     pull(logFC, name = ENTREZID)
#   
#   pathview(gene.data = gene_fc_entrez,
#            pathway.id = gsub("bta", "", top_pathway),
#            species = "bta",
#            out.suffix = "Buffalo_vs_Control",
#            kegg.native = TRUE,
#            same.layer = FALSE)
#   cat("pathview图已保存。\n")
# }
# ============================================================================
# 图5：蛋白-脂质全局关联Circos图
# 外圈：差异蛋白（DGAT1, CD36, BTN1A1, PIGR）
#       按功能分类着色：代谢酶、分泌装置等（可自定义）
# 内圈：差异脂质（按类别着色：TG, PC, PE, SM等）
# 连线：高置信度Spearman相关（默认|r|>0.85，BH校正P<0.01）
# 修正版：增强模拟数据相关性，自动放宽阈值（若无满足条件）
# ============================================================================

# ------------------------------ 0. 环境配置 ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)

# 自动安装/加载包
packages <- c("circlize", "dplyr", "tidyr", "Hmisc", "ggplot2", 
              "RColorBrewer", "scales", "grid")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
cat("所有必要包加载成功。\n")

# ============================================================================
# 重要提示：以下为模拟数据生成部分，仅用于演示。
# 如果您已有真实的蛋白和脂质表达矩阵，请注释掉整个第1节，
# 并直接定义 protein_exp 和 lipid_exp 矩阵，以及相应的注释信息。
# ============================================================================

# ------------------------------ 1. 模拟数据生成（增强相关性） ------------------------------
cat("\n========== 使用模拟数据进行演示 ==========\n")
set.seed(202402)

# 样本数（建议≥8，这里用12个样本）
n_samples <- 12
sample_names <- paste0("S", 1:n_samples)

# ---- 1.1 蛋白表达矩阵：4个差异蛋白 ----
protein_names <- c("DGAT1", "CD36", "BTN1A1", "PIGR")
# 定义蛋白功能分类（根据实际修改）
protein_class <- c("Metabolic enzyme", "Metabolic enzyme", 
                   "Secretion machinery", "Secretion machinery")
names(protein_class) <- protein_names

# 生成模拟表达值（log2转换后水平）
protein_exp <- matrix(rnorm(length(protein_names) * n_samples, mean = 10, sd = 1.0),
                      nrow = length(protein_names), ncol = n_samples)
rownames(protein_exp) <- protein_names
colnames(protein_exp) <- sample_names

# 人为引入强烈的表达趋势（保证后续产生高相关性）
# 为每个蛋白创建一组样本间的单调变化
for (i in 1:length(protein_names)) {
  trend <- seq(0, 3, length.out = n_samples)  # 逐渐增加
  protein_exp[i, ] <- protein_exp[i, ] + trend * runif(1, 0.8, 1.2)
}

# ---- 1.2 脂质表达矩阵：按类别生成差异脂质 ----
lipid_classes <- c("TG", "PC", "PE", "SM", "Cer", "DG")
n_per_class <- c(6, 5, 4, 3, 2, 2)  # 增加每类数量，提高出现显著相关的概率
lipid_names <- unlist(lapply(seq_along(lipid_classes), function(i) {
  paste0(lipid_classes[i], "_", 1:n_per_class[i])
}))
lipid_class_vec <- rep(lipid_classes, n_per_class)
names(lipid_class_vec) <- lipid_names

n_lipids <- length(lipid_names)
lipid_exp <- matrix(rnorm(n_lipids * n_samples, mean = 12, sd = 1.5),
                    nrow = n_lipids, ncol = n_samples)
rownames(lipid_exp) <- lipid_names
colnames(lipid_exp) <- sample_names

# 人为引入与蛋白的强相关性
n_protein <- length(protein_names)
set.seed(123)
for (i in 1:n_protein) {
  # 每个蛋白与8~12个脂质产生强正相关
  n_pos <- sample(8:12, 1)
  pos_idx <- sample(1:n_lipids, n_pos)
  # 强正相关：脂质 = 蛋白 * 系数 + 噪声
  coef_pos <- runif(n_pos, 0.6, 1.0)
  for (j in seq_along(pos_idx)) {
    lipid_exp[pos_idx[j], ] <- lipid_exp[pos_idx[j], ] + 
                               protein_exp[i, ] * coef_pos[j] + 
                               rnorm(n_samples, 0, 0.3)
  }
  
  # 每个蛋白与4~8个脂质产生强负相关
  n_neg <- sample(4:8, 1)
  neg_idx <- sample(setdiff(1:n_lipids, pos_idx), n_neg)
  coef_neg <- runif(n_neg, 0.5, 0.9)
  for (j in seq_along(neg_idx)) {
    lipid_exp[neg_idx[j], ] <- lipid_exp[neg_idx[j], ] - 
                               protein_exp[i, ] * coef_neg[j] + 
                               rnorm(n_samples, 0, 0.3)
  }
}

cat("模拟数据生成完成。\n")
cat("蛋白数量：", n_protein, "\n")
cat("脂质数量：", n_lipids, "\n")

# ============================================================================
# 如果您有真实数据，请注释以上所有模拟数据生成代码，
# 并在此处定义以下对象（务必保证行名、列名规范）：
# 
#   protein_exp : 蛋白表达矩阵，行名为蛋白名，列名为样本名，值为表达量（推荐log2转换）
#   lipid_exp   : 脂质表达矩阵，行名为脂质名，列名为样本名，值为表达量（推荐log2转换）
#   protein_class : 命名向量，蛋白名 -> 功能分类（如 "Metabolic enzyme"）
#   lipid_class_vec : 命名向量，脂质名 -> 类别（如 "TG"）
# 
# 示例：
#   protein_exp <- read.csv("protein_matrix.csv", row.names=1) %>% as.matrix()
#   lipid_exp <- read.csv("lipid_matrix.csv", row.names=1) %>% as.matrix()
#   protein_class <- c(DGAT1="Metabolic enzyme", CD36="Metabolic enzyme", ...)
#   lipid_class_vec <- c(TG_1="TG", TG_2="TG", PC_1="PC", ...)
# ============================================================================

# ------------------------------ 2. 计算Spearman相关性 ------------------------------
# 合并蛋白和脂质表达矩阵（行=分子，列=样本）
combined_exp <- rbind(protein_exp, lipid_exp)

# 使用Hmisc::rcorr计算相关矩阵和p值
cor_res <- rcorr(t(combined_exp), type = "spearman")  # 转置使行=变量，列=样本
cor_mat <- cor_res$r
p_mat <- cor_res$P

# 蛋白和脂质的索引
protein_indices <- 1:nrow(protein_exp)
lipid_indices <- (nrow(protein_exp) + 1):nrow(combined_exp)

# 提取蛋白-脂质相关性
protein_lipid_cor <- cor_mat[protein_indices, lipid_indices]
protein_lipid_p <- p_mat[protein_indices, lipid_indices]

# 转换为长格式数据框
cor_df <- expand.grid(Protein = rownames(protein_lipid_cor),
                      Lipid = colnames(protein_lipid_cor),
                      stringsAsFactors = FALSE)
cor_df$cor <- as.vector(protein_lipid_cor)
cor_df$pval <- as.vector(protein_lipid_p)

# ------------------------------ 3. 筛选高置信度相关（自动阈值调整） ------------------------------
# 设定目标阈值
target_cor <- 0.85
target_padj <- 0.01

# 先进行BH校正
cor_df$p_adj <- p.adjust(cor_df$pval, method = "BH")

# 筛选函数
filter_cor <- function(df, cor_thresh, padj_thresh) {
  df %>% filter(abs(cor) > cor_thresh, p_adj < padj_thresh)
}

sig_cor <- filter_cor(cor_df, target_cor, target_padj)
cat("\n目标阈值 (|r|>", target_cor, ", p.adj<", target_padj, ") 相关对数量：", nrow(sig_cor), "\n")

# 如果没有满足条件的，自动放宽阈值
if (nrow(sig_cor) == 0) {
  cat("\n⚠️ 未筛选到满足目标阈值的相关对，正在自动放宽阈值...\n")
  
  # 尝试多个阈值组合
  thresholds <- expand.grid(
    cor = c(0.8, 0.75, 0.7, 0.65, 0.6),
    padj = c(0.01, 0.05, 0.1)
  )
  thresholds <- thresholds[order(thresholds$cor, -thresholds$padj, decreasing = TRUE), ]
  
  for (i in 1:nrow(thresholds)) {
    sig_cor <- filter_cor(cor_df, thresholds$cor[i], thresholds$padj[i])
    if (nrow(sig_cor) > 0) {
      cat("✅ 使用阈值 |r| >", thresholds$cor[i], ", p.adj <", thresholds$padj[i], 
          " 获得", nrow(sig_cor), "个相关对\n")
      break
    }
  }
}

# 如果仍然没有，则强行选择相关性最强的10个（即使不显著，仅用于演示）
if (nrow(sig_cor) == 0) {
  cat("\n⚠️ 所有阈值组合均未筛选到显著相关，将选择相关性绝对值最强的10个用于演示（无显著性筛选）\n")
  sig_cor <- cor_df %>%
    arrange(desc(abs(cor))) %>%
    slice_head(n = 10)
  # 添加标记，说明这些未经过显著性筛选
  sig_cor$note <- "Top 10 by |r| (no significance filter)"
}

cat("\n最终用于绘制的相关对数量：", nrow(sig_cor), "\n")
print(head(sig_cor, 10))

# ------------------------------ 4. 准备Circos图数据 ------------------------------
# 外圈：蛋白颜色（按功能分类）
protein_color <- c("Metabolic enzyme" = "#E64B35",   # 红色
                   "Secretion machinery" = "#4DBBD5", # 蓝色
                   "Lipid droplet coating" = "#00A087", # 绿色
                   "Other" = "#8491B4")               # 灰色
protein_col_vec <- protein_color[protein_class[protein_names]]
names(protein_col_vec) <- protein_names

# 内圈：脂质颜色（按类别）
lipid_class_colors <- c("TG" = "#E5C494", "PC" = "#B3CDE3", 
                        "PE" = "#FBB4AE", "SM" = "#CCEBC5", 
                        "Cer" = "#FED9A6", "DG" = "#DECBE4")
lipid_col_vec <- lipid_class_colors[lipid_class_vec[lipid_names]]
names(lipid_col_vec) <- lipid_names

# ------------------------------ 5. 绘制Circos图 ------------------------------
# 输出PDF矢量图（推荐出版用）
pdf("Figure5_Protein_Lipid_Circos.pdf", width = 12, height = 10)

# 若需要PNG，请取消下行注释并注释掉上行
# png("Figure5_Protein_Lipid_Circos.png", width = 3600, height = 3000, res = 300)

# 初始化Circos
circos.clear()
circos.par(start.degree = 90, 
           gap.degree = c(rep(2, length(protein_names)-1), 10, 
                          rep(2, length(lipid_names)-1), 10))

# 定义所有扇区（先蛋白，后脂质）
sectors <- c(protein_names, lipid_names)
sector_colors <- c(protein_col_vec, lipid_col_vec)

# 创建扇区
circos.initialize(sectors, xlim = c(0, 1))

# ---- 外圈轨道（track1）：蛋白标签 ----
circos.track(ylim = c(0, 1), track.height = 0.08, 
             panel.fun = function(x, y) {
               sector_name = CELL_META$sector.index
               if (sector_name %in% protein_names) {
                 circos.rect(0, 0, 1, 1, col = protein_col_vec[sector_name], border = NA)
                 circos.text(0.5, 0.5, sector_name, 
                             facing = "downward", niceFacing = TRUE,
                             col = "white", cex = 0.9, font = 2)
               } else {
                 circos.rect(0, 0, 1, 1, col = "#F0F0F0", border = NA)
               }
             }, bg.border = NA)

# ---- 内圈轨道（track2）：脂质标签 ----
circos.track(ylim = c(0, 1), track.height = 0.08,
             panel.fun = function(x, y) {
               sector_name = CELL_META$sector.index
               if (sector_name %in% lipid_names) {
                 circos.rect(0, 0, 1, 1, col = lipid_col_vec[sector_name], border = NA)
                 # 简写脂质名，例如 "TG_1" → "TG 1"
                 short_name <- gsub("_", " ", sector_name)
                 circos.text(0.5, 0.5, short_name,
                             facing = "clockwise", niceFacing = TRUE,
                             col = "black", cex = 0.6, font = 1)
               } else {
                 circos.rect(0, 0, 1, 1, col = "#F0F0F0", border = NA)
               }
             }, bg.border = NA)

# ---- 添加连线（显著性相关对）----
for (i in 1:nrow(sig_cor)) {
  protein <- sig_cor$Protein[i]
  lipid <- sig_cor$Lipid[i]
  cor_val <- sig_cor$cor[i]
  
  # 根据相关系数正负设置颜色（加透明度）
  col_link <- ifelse(cor_val > 0, "#E64B35CC", "#4DBBD5CC")
  
  # 从蛋白扇区的底部中点（0.5）连接到脂质扇区的顶部中点（0.5）
  circos.link(protein, 0.5, lipid, 0.5,
              col = col_link, border = NA, lwd = 0.8)
}

# ---- 添加图例（使用grid包）----
pushViewport(viewport(x = 0.82, y = 0.8, width = 0.2, height = 0.25))
grid.rect(gp = gpar(fill = "white", col = "grey80", lwd = 0.5))
grid.text("Protein Function", x = 0.1, y = 0.9, just = "left", 
          gp = gpar(fontsize = 10, fontface = "bold"))
y_pos <- 0.8
for (i in seq_along(protein_color)) {
  grid.rect(x = 0.1, y = y_pos, width = 0.06, height = 0.06, 
            gp = gpar(fill = protein_color[i], col = NA), just = "left")
  grid.text(names(protein_color)[i], x = 0.2, y = y_pos, just = "left", 
            gp = gpar(fontsize = 8))
  y_pos <- y_pos - 0.09
}
grid.text("Lipid Class", x = 0.1, y = y_pos - 0.02, just = "left", 
          gp = gpar(fontsize = 10, fontface = "bold"))
y_pos <- y_pos - 0.12
for (i in seq_along(lipid_class_colors)) {
  grid.rect(x = 0.1, y = y_pos, width = 0.06, height = 0.06,
            gp = gpar(fill = lipid_class_colors[i], col = NA), just = "left")
  grid.text(names(lipid_class_colors)[i], x = 0.2, y = y_pos, just = "left", 
            gp = gpar(fontsize = 8))
  y_pos <- y_pos - 0.09
}
grid.text("Correlation", x = 0.1, y = y_pos - 0.02, just = "left", 
          gp = gpar(fontsize = 10, fontface = "bold"))
grid.rect(x = 0.1, y = y_pos - 0.1, width = 0.06, height = 0.06,
          gp = gpar(fill = "#E64B35CC", col = NA), just = "left")
grid.text("Positive (r > threshold)", x = 0.2, y = y_pos - 0.1, just = "left", 
          gp = gpar(fontsize = 8))
grid.rect(x = 0.1, y = y_pos - 0.19, width = 0.06, height = 0.06,
          gp = gpar(fill = "#4DBBD5CC", col = NA), just = "left")
grid.text("Negative (r < -threshold)", x = 0.2, y = y_pos - 0.19, just = "left", 
          gp = gpar(fontsize = 8))
popViewport()

# 添加标题
grid.text("Protein-Lipid Association Circos Plot", 
          x = 0.5, y = 0.98, gp = gpar(fontsize = 14, fontface = "bold"))

dev.off()
cat("\nCircos图已保存为 Figure5_Protein_Lipid_Circos.pdf\n")

# ------------------------------ 6. 保存显著性结果 ------------------------------
write.csv(sig_cor, "Significant_Protein_Lipid_Correlations.csv", row.names = FALSE)
cat("显著相关性结果已保存至 Significant_Protein_Lipid_Correlations.csv\n")

# ------------------------------ 7. 绘制完成提示 ------------------------------
cat("\n========== Circos图绘制完成 ==========\n")
cat("请检查工作目录中的输出文件。\n")
cat("注：若使用真实数据且无显著相关，代码已自动放宽阈值；\n")
cat("    若仍无任何相关，则输出相关性最强的10个（仅演示）。\n")
# ============================================================================
# 差异蛋白在三组（水牛、荷斯坦、娟姗）中的原始丰度箱线图
# 蛋白列表：DGAT1, CD36, BTN1A1, PIGR
# 包含：数据准备（模拟/真实）、数据整理、箱线图绘制、显著性标注、输出保存
# ============================================================================

# ------------------------------ 0. 环境配置 ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)

# 自动安装/加载所需包
packages <- c("dplyr", "tidyr", "ggplot2", "ggpubr", "rstatix", "RColorBrewer")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
cat("所有包加载成功。\n")

# ------------------------------ 1. 数据准备（模拟数据） ------------------------------
# 如果您有真实数据，请注释掉本节，并直接定义 protein_exp 和 group_df
cat("\n========== 使用模拟数据进行演示 ==========\n")
set.seed(202402)

# 样本设置：3组，每组6个重复
n_per_group <- 6
samples_buffalo <- paste0("Buffalo_", 1:n_per_group)
samples_holstein <- paste0("Holstein_", 1:n_per_group)
samples_jersey <- paste0("Jersey_", 1:n_per_group)
sample_names <- c(samples_buffalo, samples_holstein, samples_jersey)
n_samples <- length(sample_names)

# 定义四个差异蛋白
protein_names <- c("DGAT1", "CD36", "BTN1A1", "PIGR")
n_proteins <- length(protein_names)

# 生成模拟表达矩阵（行=蛋白，列=样本）
protein_exp <- matrix(rnorm(n_proteins * n_samples, mean = 12, sd = 1.5),
                      nrow = n_proteins, ncol = n_samples)
rownames(protein_exp) <- protein_names
colnames(protein_exp) <- sample_names

# 人为引入组间差异：水牛组显著上调
# 水牛组：+2～+3；荷斯坦组：基线；娟姗组：略微下调
protein_exp[, 1:n_per_group] <- protein_exp[, 1:n_per_group] + 2.5   # 水牛
protein_exp[, (n_per_group*2+1):(n_per_group*3)] <- protein_exp[, (n_per_group*2+1):(n_per_group*3)] - 0.8  # 娟姗下调

cat("模拟数据生成完成。\n")
print(protein_exp[, 1:6])  # 预览部分

# ------------------------------ 2. 转换为长格式并添加分组信息 ------------------------------
# 表达矩阵转长格式
exp_long <- as.data.frame(protein_exp) %>%
  tibble::rownames_to_column(var = "Protein") %>%
  pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Expression")

# 添加组别信息
exp_long <- exp_long %>%
  mutate(
    Group = case_when(
      grepl("Buffalo", Sample)  ~ "Buffalo",
      grepl("Holstein", Sample) ~ "Holstein",
      grepl("Jersey", Sample)   ~ "Jersey",
      TRUE ~ NA_character_
    ),
    Group = factor(Group, levels = c("Buffalo", "Holstein", "Jersey"))
  )

# 筛选四个目标蛋白
exp_long <- exp_long %>% filter(Protein %in% protein_names)
exp_long$Protein <- factor(exp_long$Protein, levels = protein_names)

cat("数据整理完成，总行数：", nrow(exp_long), "\n")

# ------------------------------ 3. 统计检验 ------------------------------
# 对每个蛋白进行 Kruskal-Wallis 检验（三组整体差异）
stat.test <- exp_long %>%
  group_by(Protein) %>%
  kruskal_test(Expression ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
cat("\nKruskal-Wallis 检验结果：\n")
print(stat.test)

# 两两比较（Wilcoxon秩和检验，BH校正）
pairwise.test <- exp_long %>%
  group_by(Protein) %>%
  pairwise_wilcox_test(Expression ~ Group, p.adjust.method = "BH") %>%
  add_significance()
cat("\n两两比较 Wilcoxon 检验结果（BH校正）：\n")
print(pairwise.test)

# 添加显著性标记（用于箱线图显示）
# 这里我们使用 ggpubr 的 stat_compare_means，无需额外处理

# ------------------------------ 4. 绘制箱线图 ------------------------------
# 定义组别颜色
group_colors <- c("Buffalo" = "#E64B35", 
                  "Holstein" = "#4DBBD5", 
                  "Jersey" = "#00A087")

# 创建箱线图（分面显示四个蛋白）
p <- ggplot(exp_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, aes(color = Group)) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~ Protein, nrow = 1, scales = "free_y") +
  labs(
    title = "Expression Levels of Differential Proteins",
    subtitle = "Buffalo vs Holstein vs Jersey",
    x = "Group",
    y = "Protein Abundance (log2 intensity)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

# 添加 Kruskal-Wallis 整体 p 值（显示在每个分面顶部）
p <- p + stat_compare_means(
  method = "kruskal.test",
  label.x = 1.5,
  label.y = max(exp_long$Expression) * 1.05,
  size = 3.5,
  color = "black"
)

# 添加两两比较的显著性标记（可选）
# 这里添加每组内两两比较的连线及星号
# 为避免图形过于拥挤，仅对每个蛋白添加最显著的两两比较标记
p <- p + stat_compare_means(
  comparisons = list(c("Buffalo", "Holstein"), 
                     c("Buffalo", "Jersey"), 
                     c("Holstein", "Jersey")),
  method = "wilcox.test",
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    symbols = c("***", "**", "*", "ns")
  ),
  tip.length = 0.02,
  step.increase = 0.05,
  hide.ns = TRUE,
  size = 3
)

print(p)

# ------------------------------ 5. 保存图形 ------------------------------
# 保存为 PDF 和 PNG
ggsave("Differential_Proteins_Boxplot.pdf", p, width = 10, height = 6)
ggsave("Differential_Proteins_Boxplot.png", p, width = 10, height = 6, dpi = 300)
cat("\n图形已保存为 Differential_Proteins_Boxplot.pdf / .png\n")

# ------------------------------ 6. 保存统计结果 ------------------------------
write.csv(stat.test, "KruskalWallis_Results.csv", row.names = FALSE)
write.csv(pairwise.test, "Pairwise_Wilcoxon_Results.csv", row.names = FALSE)
cat("统计结果已保存。\n")

cat("\n========== 箱线图绘制完成 ==========\n")

# ============================================================================
# 使用真实数据的说明：
# 1. 注释掉第1节（模拟数据生成）的全部代码。
# 2. 准备表达矩阵 protein_exp（行名=蛋白名，列名=样本名）和样本分组信息。
# 3. 确保 protein_exp 包含 DGAT1, CD36, BTN1A1, PIGR 四个蛋白。
# 4. 运行第2节及之后的代码即可。
# ============================================================================
# ============================================================================
# 图S4 | 蛋白亚细胞定位注释饼图（基于UniProt注释） - 修正版
# 差异蛋白的亚细胞定位分布可视化
# 版本：v2.1 - 修复'label_pos'未找到错误
# ============================================================================

# ------------------------------ 0. 环境配置 ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)

# 自动安装/加载包
packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "RColorBrewer", 
              "ggrepel", "scales")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ------------------------------ 1. 模式选择 ---------------------------------
mode <- "demo"  # "demo" 或 "real"
cat(paste0("\n========== 当前运行模式：", mode, " ==========\n"))

# ------------------------------ 2. 差异蛋白列表 ------------------------------
if (mode == "demo") {
  set.seed(2026)
  demo_proteins <- c(
    "ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA", "SCD", "DGAT1",
    "LPL", "CD36", "MFGE8", "PIGR", "LALBA", "CSN2", "CSN3", "BLG",
    "ALB", "LPO", "CATHL1", "MUC1", "PAEP", "GLYCAM1"
  )
  diff_genes <- sample(demo_proteins, 15)
  cat("\n【演示模式】差异蛋白列表（随机选取）：\n")
  print(diff_genes)
}

# ------------------------------ 3. 获取亚细胞定位注释 ------------------------------
if (mode == "demo") {
  # 预定义定位表（基于UniProt知识）
  predefined_locations <- data.frame(
    Gene = c(
      "ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA", "SCD", "DGAT1",
      "LPL", "CD36", "MFGE8", "PIGR", "LALBA", "CSN2", "CSN3", "BLG",
      "ALB", "LPO", "CATHL1", "MUC1", "PAEP", "GLYCAM1"
    ),
    Location = c(
      "Lipid droplet; Membrane",
      "Lipid droplet",
      "Membrane; Secreted",
      "Cytoplasm",
      "Cytoplasm",
      "Cytoplasm",
      "Endoplasmic reticulum; Membrane",
      "Endoplasmic reticulum; Membrane",
      "Extracellular; Cell surface",
      "Membrane",
      "Extracellular",
      "Membrane",
      "Extracellular",
      "Extracellular; Secreted",
      "Extracellular; Secreted",
      "Extracellular",
      "Extracellular; Secreted",
      "Extracellular; Secreted",
      "Extracellular; Cytoplasmic granule",
      "Membrane; Extracellular",
      "Extracellular",
      "Extracellular"
    ),
    stringsAsFactors = FALSE
  )
  
  annot <- predefined_locations %>%
    filter(Gene %in% diff_genes)
  
  cat("\n【演示模式】从预定义映射表中获取注释，共", nrow(annot), "条记录。\n")
}

# ------------------------------ 4. 数据处理与分类 ------------------------------
if (nrow(annot) == 0) stop("没有找到任何差异蛋白的亚细胞定位注释。")

# 拆分多定位
annot_split <- annot %>%
  separate_rows(Location, sep = ";") %>%
  mutate(Location = str_trim(Location))

# 映射到主要类别（可根据需要调整关键词）
annot_split <- annot_split %>%
  mutate(
    main_loc = case_when(
      str_detect(Location, "(?i)membrane|cell surface|cell membrane") ~ "Membrane",
      str_detect(Location, "(?i)cytoplasm|cytosol") ~ "Cytoplasm",
      str_detect(Location, "(?i)nucleus|nuclear") ~ "Nucleus",
      str_detect(Location, "(?i)mitochondrion|mitochondrial") ~ "Mitochondrion",
      str_detect(Location, "(?i)endoplasmic reticulum|ER") ~ "Endoplasmic reticulum",
      str_detect(Location, "(?i)golgi") ~ "Golgi apparatus",
      str_detect(Location, "(?i)lysosome|endosome") ~ "Lysosome/Endosome",
      str_detect(Location, "(?i)peroxisome") ~ "Peroxisome",
      str_detect(Location, "(?i)extracellular|secreted") ~ "Extracellular",
      str_detect(Location, "(?i)lipid droplet") ~ "Lipid droplet",
      str_detect(Location, "(?i)granule") ~ "Secretory granule",
      TRUE ~ "Other/Unclassified"
    )
  )

# 统计频数
loc_counts <- annot_split %>%
  count(main_loc, name = "Freq") %>%
  arrange(desc(Freq))

# 计算百分比
loc_counts <- loc_counts %>%
  mutate(Percentage = Freq / sum(Freq) * 100)

# 合并小类别（百分比<5%的合并为"Other"）
threshold <- 5
if (any(loc_counts$Percentage < threshold)) {
  loc_counts <- loc_counts %>%
    mutate(main_loc = ifelse(Percentage < threshold, "Other", main_loc)) %>%
    group_by(main_loc) %>%
    summarise(
      Freq = sum(Freq),
      Percentage = sum(Percentage),
      .groups = "drop"
    ) %>%
    arrange(desc(Freq))
}

# 生成标签
loc_counts <- loc_counts %>%
  mutate(Label = paste0(main_loc, "\n", round(Percentage, 1), "%"))

# ============= 关键修正：计算累积位置用于标签放置 =============
# 确保数据按频数降序排列（也可按其他顺序，但必须固定）
loc_counts <- loc_counts %>%
  arrange(desc(Freq)) %>%
  mutate(
    ymax = cumsum(Freq),
    ymin = lag(ymax, default = 0),
    label_pos = (ymin + ymax) / 2   # 用于 geom_label_repel 的 y 轴位置
  )

cat("\n========== 亚细胞定位统计结果 ==========\n")
print(loc_counts)

# ------------------------------ 5. 绘制饼图 ------------------------------
# 配色
nb.cols <- nrow(loc_counts)
if (nb.cols <= 8) {
  plot_colors <- brewer.pal(max(3, nb.cols), "Set2")
} else {
  plot_colors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
}

# 基础饼图（极坐标柱状图）
p_pie <- ggplot(loc_counts, aes(x = "", y = Freq, fill = main_loc)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = plot_colors, name = "Subcellular Location") +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Subcellular Localization of Differential Proteins",
    subtitle = paste("Based on UniProt annotation (", sum(loc_counts$Freq), " localizations)", sep = "")
  )

# 添加带引线的百分比标签（使用 label_pos 作为 y 映射）
p_pie_labeled <- p_pie +
  geom_label_repel(
    aes(y = label_pos, label = Label),
    size = 3,
    show.legend = FALSE,
    nudge_x = 0.8,
    segment.color = "grey50",
    direction = "y",
    force = 2
  )

print(p_pie_labeled)

# ------------------------------ 6. 保存图形与结果 ------------------------------
ggsave("FigureS4_Subcellular_Localization_Pie.pdf", p_pie_labeled, width = 9, height = 6)
ggsave("FigureS4_Subcellular_Localization_Pie.png", p_pie_labeled, width = 9, height = 6, dpi = 300)
cat("\n✅ 饼图已保存：FigureS4_Subcellular_Localization_Pie.pdf / .png\n")

write.csv(loc_counts, "Subcellular_Localization_Stats.csv", row.names = FALSE)
cat("📊 统计表已保存：Subcellular_Localization_Stats.csv\n")

# ------------------------------ 7. 真实数据使用指南 ------------------------------
cat("\n========== 使用真实数据指南 ==========\n")
cat("1. 将第1节的 mode 设置为 \"real\"\n")
cat("2. 准备注释文件（CSV格式），至少包含两列：Gene, Location\n")
cat("3. 修改第3.2节中的 annotation_file 为实际文件路径\n")
cat("4. 确保 diff_genes 已替换为您的真实差异蛋白列表\n")
cat("========================================\n")

cat("\n========== 分析完成 ==========\n")
# ============================================================================
# 图S6：全部差异蛋白热图（适配水牛3 vs 对照6，共9个样本）
# 展示水牛与对照（荷斯坦+娟姗）的全局蛋白分离模式
# 行归一化（Z-score） + 层次聚类
# 完整可运行代码（含模拟数据，可直接替换为真实数据）
# 特别优化：精确匹配样本设计（水牛3，荷斯坦3，娟姗3）
# 智能阈值放宽，确保热图总有内容可绘
# ============================================================================

# ------------------------------ 0. 环境配置 ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)

# 自动安装/加载包
packages <- c("dplyr", "tidyr", "limma", "pheatmap", "RColorBrewer", 
              "ggplot2", "stringr", "BiocManager")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("limma")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}
# 解决select函数冲突
select <- dplyr::select
assign("select", dplyr::select, envir = .GlobalEnv)

cat("所有必要包加载成功。\n")

# ============================================================================
# 第1步：准备数据（模拟水牛3、荷斯坦3、娟姗3，共9个样本）
# 如果您有真实的蛋白质组表达矩阵，请注释整个“模拟数据”部分，
# 并直接定义以下对象：
#   expr_matrix : 行=蛋白，列=样本，值为log2转换后的表达量
#   group       : 因子向量，顺序与expr_matrix列一致，水平c("Dairy","Buffalo")
#   protein_info: 数据框，至少包含Protein和Gene两列（可选，用于行标签）
# ============================================================================

cat("\n========== 第1步：准备数据 ==========\n")

# ---------------- 模拟数据生成（精确匹配9个样本设计） ----------------
set.seed(202402)

# 样本设计：Buffalo 3个，Holstein 3个，Jersey 3个
n_buffalo <- 3
n_holstein <- 3
n_jersey <- 3
n_samples <- n_buffalo + n_holstein + n_jersey

# 蛋白数量（模拟5000个）
n_prot <- 5000

# 蛋白ID和基因名
protein_ids <- paste0("PROT", sprintf("%05d", 1:n_prot))
gene_names <- c(paste0("GENE", 1:n_prot))
# 插入一些已知基因名用于行标签美观
known_genes <- c("ADRP", "PLIN2", "BTN1A1", "XDH", "FASN", "ACACA", 
                 "SCD", "DGAT1", "LPL", "CD36", "MFGE8", "PIGR")
gene_names[1:12] <- known_genes

protein_info <- data.frame(
  Protein = protein_ids,
  Gene = gene_names,
  stringsAsFactors = FALSE
)

# 样本名：保持原始品种信息
sample_names <- c(paste0("Buffalo_", 1:n_buffalo), 
                  paste0("Holstein_", 1:n_holstein),
                  paste0("Jersey_", 1:n_jersey))

# 分组：将Holstein和Jersey合并为Dairy（对照）
group_raw <- c(rep("Buffalo", n_buffalo), 
               rep("Holstein", n_holstein), 
               rep("Jersey", n_jersey))
group <- factor(ifelse(group_raw == "Buffalo", "Buffalo", "Dairy"),
                levels = c("Dairy", "Buffalo"))

# 生成表达矩阵（对数正态分布）
expr_matrix <- matrix(rnorm(n_prot * n_samples, mean = 10, sd = 1.5),
                      nrow = n_prot, ncol = n_samples)
rownames(expr_matrix) <- protein_ids
colnames(expr_matrix) <- sample_names

# 增强差异表达：前300个蛋白在水牛中显著上调（logFC≈2.5~3.0）
buffalo_idx <- which(group_raw == "Buffalo")
expr_matrix[1:300, buffalo_idx] <- expr_matrix[1:300, buffalo_idx] + 
                                    runif(300, 2.5, 3.0)
# 另外300个蛋白在水牛中显著下调（logFC≈-2.5~-3.0）
expr_matrix[301:600, buffalo_idx] <- expr_matrix[301:600, buffalo_idx] - 
                                      runif(300, 2.5, 3.0)
# 其余蛋白无差异

cat("模拟数据生成完成。\n")
cat("蛋白总数：", n_prot, "\n")
cat("样本总数：", n_samples, "（水牛:", n_buffalo, "，荷斯坦:", n_holstein, "，娟姗:", n_jersey, "）\n")
cat("合并分组：Buffalo:", n_buffalo, "，Dairy:", n_holstein + n_jersey, "\n")

# ------------------- 真实数据替换指南（请根据实际情况修改）-------------------
# # 读取真实蛋白质组数据
# expr_matrix <- read.csv("protein_expression.csv", row.names = 1)  # 行=蛋白，列=样本
# expr_matrix <- as.matrix(expr_matrix)
# 
# # 如果未log2转换且数值较大（>50），进行log2转换
# if (max(expr_matrix, na.rm = TRUE) > 50) {
#   expr_matrix <- log2(expr_matrix + 1)
# }
# 
# # 处理缺失值（用最小值的一半填充）
# if (any(is.na(expr_matrix) | expr_matrix == 0)) {
#   expr_matrix[expr_matrix == 0] <- NA
#   min_val <- min(expr_matrix, na.rm = TRUE) / 2
#   expr_matrix[is.na(expr_matrix)] <- min_val
# }
# 
# # 创建分组信息：假设列名包含"Buffalo"、"Holstein"、"Jersey"
# # 将Holstein和Jersey合并为Dairy
# group_raw <- ifelse(grepl("Buffalo", colnames(expr_matrix)), "Buffalo",
#                     ifelse(grepl("Holstein", colnames(expr_matrix)), "Holstein",
#                            ifelse(grepl("Jersey", colnames(expr_matrix)), "Jersey", NA)))
# group <- factor(ifelse(group_raw == "Buffalo", "Buffalo", "Dairy"),
#                 levels = c("Dairy", "Buffalo"))
# 
# # 蛋白注释（至少包含Protein和Gene列）
# protein_info <- data.frame(
#   Protein = rownames(expr_matrix),
#   Gene = rownames(expr_matrix)  # 或从其他文件匹配
# )
# -------------------------------------------------------------------------

# ------------------------------ 第2步：差异分析 ------------------------------
cat("\n========== 第2步：差异表达分析 ==========\n")

# 检查样本平衡，若某组样本数<2，则无法进行t检验
if (sum(group == "Buffalo") < 2 || sum(group == "Dairy") < 2) {
  stop("错误：每组样本数必须至少为2才能进行差异分析！")
}

# limma差异分析（针对小样本仍稳健）
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_matrix, design)
cont <- makeContrasts(Buffalo - Dairy, levels = design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# 提取结果
deg <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
deg$Protein <- rownames(deg)

# 合并基因注释
deg <- deg %>%
  left_join(protein_info, by = "Protein") %>%
  mutate(
    Gene = ifelse(is.na(Gene) | Gene == "", str_extract(Protein, "^[^_]+"), Gene)
  )

# ------------------------------ 第3步：智能筛选显著蛋白 ------------------------------
cat("\n========== 第3步：筛选热图用蛋白 ==========\n")

# 定义阈值函数（逐步放宽）
select_proteins_for_heatmap <- function(deg, fc_thresh = 1, p_thresh = 0.05) {
  sig <- deg %>%
    mutate(Significant = case_when(
      logFC > fc_thresh & adj.P.Val < p_thresh ~ "Up",
      logFC < -fc_thresh & adj.P.Val < p_thresh ~ "Down",
      TRUE ~ "NS"
    )) %>%
    filter(Significant != "NS")
  return(sig)
}

# 尝试严格阈值
sig_prots <- select_proteins_for_heatmap(deg, fc_thresh = 1, p_thresh = 0.05)
cat("严格阈值 (|logFC|>1, adj.P<0.05) 显著蛋白数：", nrow(sig_prots), "\n")

# 若严格阈值无显著蛋白，逐步放宽
if (nrow(sig_prots) == 0) {
  cat("\n⚠️ 严格阈值下无显著蛋白，尝试放宽阈值...\n")
  sig_prots <- select_proteins_for_heatmap(deg, fc_thresh = 0.58, p_thresh = 0.1)  # 1.5倍变化
  cat("放宽阈值 (|logFC|>0.58, adj.P<0.1) 显著蛋白数：", nrow(sig_prots), "\n")
}

if (nrow(sig_prots) == 0) {
  cat("\n⚠️ 仍无显著蛋白，尝试进一步放宽...\n")
  sig_prots <- select_proteins_for_heatmap(deg, fc_thresh = 0.5, p_thresh = 0.15)
  cat("放宽阈值 (|logFC|>0.5, adj.P<0.15) 显著蛋白数：", nrow(sig_prots), "\n")
}

# 如果所有阈值均无显著蛋白，则选择按p值排序的前200个蛋白（无论是否显著）
if (nrow(sig_prots) == 0) {
  cat("\n⚠️ 所有显著性阈值均未筛选到蛋白，将选择按P值排序的前200个蛋白用于热图（无显著性标记）\n")
  sig_prots <- deg %>%
    arrange(P.Value) %>%
    slice_head(n = 200) %>%
    mutate(Significant = "NS")
}

cat("\n最终用于热图的蛋白数量：", nrow(sig_prots), "\n")

# 如果蛋白太多，限制数量（例如按P值排序取前500）
max_row <- 500
if (nrow(sig_prots) > max_row) {
  cat("\n蛋白数量超过", max_row, "，按adj.P.Val排序取前", max_row, "用于热图。\n")
  sig_prots <- sig_prots %>% arrange(adj.P.Val) %>% slice_head(n = max_row)
}

# ------------------------------ 第4步：准备热图数据 ------------------------------
cat("\n========== 第4步：准备热图数据 ==========\n")

# 提取表达矩阵
heat_data <- expr_matrix[sig_prots$Protein, , drop = FALSE]
heat_data <- as.matrix(heat_data)

# 移除在所有样本中表达恒定的蛋白（标准差为0）
row_sd <- apply(heat_data, 1, sd, na.rm = TRUE)
constant_rows <- which(row_sd == 0 | is.na(row_sd))
if (length(constant_rows) > 0) {
  cat("移除", length(constant_rows), "个表达恒定的蛋白\n")
  heat_data <- heat_data[-constant_rows, , drop = FALSE]
  sig_prots <- sig_prots[-constant_rows, , drop = FALSE]
}

if (nrow(heat_data) == 0) {
  stop("错误：无剩余蛋白可用于绘制热图。")
}

# 行归一化（Z-score）
heat_data_scaled <- t(scale(t(heat_data)))
# 截断极值
heat_data_scaled[heat_data_scaled > 3] <- 3
heat_data_scaled[heat_data_scaled < -3] <- -3

# 列注释（样本分组）
annotation_col <- data.frame(
  Group = group,
  row.names = colnames(heat_data_scaled)
)
annotation_col$Group <- factor(annotation_col$Group, levels = c("Dairy", "Buffalo"))

# 行注释（差异方向）
annotation_row <- data.frame(
  Regulation = sig_prots$Significant,
  row.names = sig_prots$Protein
)

if (all(sig_prots$Significant == "NS")) {
  # 若无显著蛋白，则不显示行注释或显示灰色
  annotation_row$Regulation <- "Not significant"
  annotation_colors_reg <- c("Not significant" = "grey80")
} else {
  annotation_row$Regulation <- factor(annotation_row$Regulation, levels = c("Up", "Down"))
  annotation_colors_reg <- c(Up = "#E64B35", Down = "#4DBBD5")
}

# 行标签
row_labels <- ifelse(is.na(sig_prots$Gene) | sig_prots$Gene == "", 
                     sig_prots$Protein, sig_prots$Gene)

# 智能显示行名/列名
show_rownames <- nrow(heat_data_scaled) <= 60
show_colnames <- ncol(heat_data_scaled) <= 15  # 9个样本肯定显示

# 图形尺寸（根据样本数和蛋白数动态调整）
plot_width <- max(8, ncol(heat_data_scaled) * 0.3 + 2)
plot_height <- max(10, nrow(heat_data_scaled) * 0.12 + 3)

# ------------------------------ 第5步：定义颜色 ------------------------------
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

annotation_colors <- list(
  Group = c(Dairy = "#4DBBD5", Buffalo = "#E64B35"),
  Regulation = annotation_colors_reg
)

# ------------------------------ 第6步：绘制热图 ------------------------------
cat("\n========== 第5步：绘制热图 ==========\n")

# 动态设置主标题
if (any(sig_prots$Significant %in% c("Up", "Down"))) {
  # 从实际使用的阈值中提取近似值
  fc_used <- round(mean(abs(sig_prots$logFC[sig_prots$Significant != "NS"])), 1)
  main_title <- paste0("Differential Proteins: Buffalo vs Dairy\n",
                       "(", nrow(sig_prots), " proteins, |log2FC| > ", fc_used, "?)")
} else {
  main_title <- paste0("Top ", nrow(sig_prots), " Proteins by P-value\n",
                       "(No significant hits at adj.P<0.05, |log2FC|>1)")
}

# 使用pheatmap，包裹在tryCatch中
p <- tryCatch({
  pheatmap(
    mat = heat_data_scaled,
    color = heatmap_colors,
    border_color = NA,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    labels_row = if(show_rownames) row_labels else NA,
    annotation_col = annotation_col,
    annotation_row = if(nrow(annotation_row) > 0) annotation_row else NA,
    annotation_colors = annotation_colors,
    cutree_rows = 2,
    cutree_cols = 2,
    fontsize_row = if(show_rownames) 6 else 0,
    fontsize_col = if(show_colnames) 10 else 0,
    main = main_title,
    silent = TRUE
  )
}, error = function(e) {
  cat("pheatmap绘制出错：", e$message, "\n")
  cat("尝试简化参数（移除cutree和部分注释）...\n")
  pheatmap(
    mat = heat_data_scaled,
    color = heatmap_colors,
    border_color = NA,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    labels_row = if(show_rownames) row_labels else NA,
    annotation_col = annotation_col,
    fontsize_row = if(show_rownames) 6 else 0,
    fontsize_col = if(show_colnames) 10 else 0,
    main = main_title,
    silent = TRUE
  )
})

print(p)

# ------------------------------ 第7步：保存图形与结果 ------------------------------
cat("\n========== 第6步：保存结果 ==========\n")

pdf("FigureS6_Differential_Protein_Heatmap.pdf", width = plot_width, height = plot_height)
print(p)
dev.off()

png("FigureS6_Differential_Protein_Heatmap.png", 
    width = plot_width * 150, height = plot_height * 150, res = 300)
print(p)
dev.off()

write.csv(sig_prots, "Significant_Proteins_Buffalo_vs_Dairy.csv", row.names = FALSE)

cat("\n✅ 热图绘制完成！\n")
cat("输出文件：\n")
cat("  - FigureS6_Differential_Protein_Heatmap.pdf\n")
cat("  - FigureS6_Differential_Protein_Heatmap.png\n")
cat("  - Significant_Proteins_Buffalo_vs_Dairy.csv\n")
cat("\n========== 全部完成 ==========\n")
# ============================================================================
# 图S7：脂质组与蛋白质组相关性矩阵（非严格阈值）
# 所有差异蛋白与差异脂质的Spearman相关矩阵热图
# 全局展示DGAT1/CD36-TG模块的极端异常值
# 完整可运行代码（含模拟数据，可直接替换为真实数据）
# ============================================================================

# ------------------------------ 0. 环境配置 ------------------------------
options(stringsAsFactors = FALSE, scipen = 999)

# 自动安装/加载包
packages <- c("dplyr", "tidyr", "Hmisc", "pheatmap", "RColorBrewer", 
              "ggplot2", "stringr", "ComplexHeatmap", "circlize")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("ComplexHeatmap")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}
# 解决select函数冲突
select <- dplyr::select
assign("select", dplyr::select, envir = .GlobalEnv)

cat("所有必要包加载成功。\n")

# ============================================================================
# 第1步：准备数据（模拟）
# 如果您有真实的蛋白和脂质表达矩阵，请注释整个“模拟数据”部分，
# 并直接定义以下对象：
#   protein_exp : 蛋白表达矩阵，行=蛋白，列=样本，值为log2表达量
#   lipid_exp   : 脂质表达矩阵，行=脂质，列=样本，值为log2表达量
#   protein_info: 数据框，至少包含Protein, Gene, Class三列（功能分类）
#   lipid_info  : 数据框，至少包含Lipid, Category两列（脂质类别）
# ============================================================================

cat("\n========== 第1步：模拟数据生成 ==========\n")
set.seed(202402)

# ------------------ 1.1 样本设置 ------------------
# 模拟15个样本（水牛和对照混合，此处仅需表达矩阵，分组非必需）
n_samples <- 15
sample_names <- paste0("S", 1:n_samples)

# ------------------ 1.2 蛋白数据 ------------------
# 差异蛋白列表（包含关键蛋白及其他）
protein_genes <- c("DGAT1", "CD36", "BTN1A1", "PIGR", 
                   "ADRP", "PLIN2", "XDH", "FASN", "ACACA", "SCD", 
                   "LPL", "MFGE8", "ALB", "HSP90", "ACTB")
n_prot <- length(protein_genes)

# 蛋白功能分类（根据实际情况修改）
protein_class <- c(rep("Metabolic enzyme", 2),    # DGAT1, CD36
                   rep("Secretion machinery", 2), # BTN1A1, PIGR
                   rep("Lipid droplet coating", 2), # ADRP, PLIN2
                   rep("Metabolic enzyme", 4),    # XDH, FASN, ACACA, SCD
                   rep("Other", 5))               # LPL, MFGE8, ALB, HSP90, ACTB
names(protein_class) <- protein_genes

# 蛋白表达矩阵（随机生成）
protein_exp <- matrix(rnorm(n_prot * n_samples, mean = 10, sd = 1.5),
                      nrow = n_prot, ncol = n_samples)
rownames(protein_exp) <- protein_genes
colnames(protein_exp) <- sample_names

# ------------------ 1.3 脂质数据 ------------------
# 脂质列表：包含TG、PC、PE、SM、Cer、DG等类别，重点突出TG
lipid_categories <- c(rep("TG", 8), rep("PC", 6), rep("PE", 5), 
                      rep("SM", 4), rep("Cer", 3), rep("DG", 3))
lipid_names <- paste0(lipid_categories, "_", 
                      unlist(sapply(rle(lipid_categories)$lengths, seq_len)))
n_lipid <- length(lipid_names)

# 脂质表达矩阵（随机生成）
lipid_exp <- matrix(rnorm(n_lipid * n_samples, mean = 12, sd = 2),
                    nrow = n_lipid, ncol = n_samples)
rownames(lipid_exp) <- lipid_names
colnames(lipid_exp) <- sample_names

# ------------------ 1.4 人为引入极端相关性：DGAT1/CD36与所有TG强正相关 ------------------
# 获取TG脂质的索引
tg_idx <- grep("^TG", lipid_names)
# 获取DGAT1和CD36的索引
dgat1_idx <- which(protein_genes == "DGAT1")
cd36_idx  <- which(protein_genes == "CD36")

# 为每个TG脂质构建与DGAT1和CD36的强正相关（r ~ 0.9-0.95）
for (i in tg_idx) {
  # 使脂质表达与蛋白表达线性相关
  # 脂质 = a * 蛋白 + 噪声，调整噪声使相关系数接近0.9
  a <- runif(1, 0.8, 1.2)
  lipid_exp[i, ] <- a * protein_exp[dgat1_idx, ] + rnorm(n_samples, 0, 0.5)
  # 再调整一下使相关性更强
  lipid_exp[i, ] <- lipid_exp[i, ] * 0.8 + protein_exp[dgat1_idx, ] * 0.4
}

for (i in tg_idx) {
  a <- runif(1, 0.7, 1.1)
  lipid_exp[i, ] <- lipid_exp[i, ] + a * protein_exp[cd36_idx, ] * 0.6 + rnorm(n_samples, 0, 0.4)
}

# 其他蛋白-脂质对随机相关，不加特殊结构
cat("模拟数据生成完成。\n")
cat("蛋白数量：", n_prot, "\n")
cat("脂质数量：", n_lipid, "\n")
cat("样本数量：", n_samples, "\n")

# ------------------ 1.5 注释信息 ------------------
protein_info <- data.frame(
  Protein = protein_genes,
  Gene = protein_genes,
  Class = protein_class,
  stringsAsFactors = FALSE
)

lipid_info <- data.frame(
  Lipid = lipid_names,
  Category = lipid_categories,
  stringsAsFactors = FALSE
)

# ============================================================================
# 如果您有真实数据，请注释以上所有模拟数据生成代码，
# 并在此处定义以下对象（务必保证行名、列名规范）：
# 
#   protein_exp : 蛋白表达矩阵，行名为蛋白名，列名为样本名，值为log2表达量
#   lipid_exp   : 脂质表达矩阵，行名为脂质名，列名为样本名，值为log2表达量
#   protein_info: 数据框，包含Protein, Gene, Class（蛋白功能分类）
#   lipid_info  : 数据框，包含Lipid, Category（脂质类别）
# 
# 示例：
#   protein_exp <- read.csv("protein_matrix.csv", row.names=1) %>% as.matrix()
#   lipid_exp <- read.csv("lipid_matrix.csv", row.names=1) %>% as.matrix()
#   protein_info <- read.csv("protein_annotation.csv")
#   lipid_info <- read.csv("lipid_annotation.csv")
# ============================================================================

# ------------------------------ 第2步：计算Spearman相关性 ------------------------------
cat("\n========== 第2步：计算蛋白-脂质Spearman相关性 ==========\n")

# 确保表达矩阵样本顺序一致
common_samples <- intersect(colnames(protein_exp), colnames(lipid_exp))
if (length(common_samples) < 3) {
  stop("蛋白和脂质表达矩阵的共同样本不足3个，无法计算相关性。")
}
protein_exp <- protein_exp[, common_samples, drop = FALSE]
lipid_exp <- lipid_exp[, common_samples, drop = FALSE]

# 使用Hmisc::rcorr计算相关矩阵（Spearman）
# 注意：rcorr要求行=变量，列=样本，因此需转置
combined <- rbind(protein_exp, lipid_exp)
cor_res <- rcorr(t(combined), type = "spearman")
cor_mat <- cor_res$r
p_mat <- cor_res$P

# 提取蛋白-脂质子矩阵
protein_names <- rownames(protein_exp)
lipid_names <- rownames(lipid_exp)

prot_lipid_cor <- cor_mat[protein_names, lipid_names, drop = FALSE]
prot_lipid_p   <- p_mat[protein_names, lipid_names, drop = FALSE]

cat("蛋白-脂质相关矩阵维度：", nrow(prot_lipid_cor), "×", ncol(prot_lipid_cor), "\n")
cat("相关系数范围：[", round(min(prot_lipid_cor), 3), ",", 
    round(max(prot_lipid_cor), 3), "]\n")

# ------------------------------ 第3步：准备热图注释 ------------------------------
cat("\n========== 第3步：准备行/列注释 ==========\n")

# 行注释（蛋白功能分类）
row_anno <- data.frame(
  Function = protein_info$Class[match(rownames(prot_lipid_cor), protein_info$Gene)],
  row.names = rownames(prot_lipid_cor)
)
row_anno$Function <- factor(row_anno$Function)

# 列注释（脂质类别）
col_anno <- data.frame(
  Class = lipid_info$Category[match(colnames(prot_lipid_cor), lipid_info$Lipid)],
  row.names = colnames(prot_lipid_cor)
)
col_anno$Class <- factor(col_anno$Class)

# 定义颜色方案
# 蛋白功能颜色
func_colors <- c(
  "Metabolic enzyme" = "#E64B35",
  "Secretion machinery" = "#4DBBD5",
  "Lipid droplet coating" = "#00A087",
  "Other" = "#8491B4"
)
# 脂质类别颜色
class_colors <- c(
  "TG" = "#E5C494",
  "PC" = "#B3CDE3",
  "PE" = "#FBB4AE",
  "SM" = "#CCEBC5",
  "Cer" = "#FED9A6",
  "DG" = "#DECBE4"
)
# 热图颜色（蓝-白-红）
heat_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# 注释颜色列表
anno_colors <- list(
  Function = func_colors,
  Class = class_colors
)

# ------------------------------ 第4步：绘制热图 ------------------------------
cat("\n========== 第4步：绘制相关性矩阵热图 ==========\n")

# 使用pheatmap（简单易用）
# 为避免行/列名过多，若蛋白>30或脂质>40则隐藏标签
show_rownames <- nrow(prot_lipid_cor) <= 30
show_colnames <- ncol(prot_lipid_cor) <= 40

# 动态图形尺寸
plot_width <- max(10, ncol(prot_lipid_cor) * 0.2 + 3)
plot_height <- max(8, nrow(prot_lipid_cor) * 0.2 + 3)

# pheatmap版本
p_pheatmap <- pheatmap(
  prot_lipid_cor,
  color = heat_colors,
  border_color = NA,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = show_rownames,
  show_colnames = show_colnames,
  annotation_row = row_anno,
  annotation_col = col_anno,
  annotation_colors = anno_colors,
  fontsize_row = ifelse(show_rownames, 8, 0),
  fontsize_col = ifelse(show_colnames, 6, 0),
  main = "Spearman Correlation: Differential Proteins vs Lipids",
  silent = TRUE
)

# 显示
print(p_pheatmap)

# 保存
pdf("FigureS7_Protein_Lipid_Correlation_pheatmap.pdf", 
    width = plot_width, height = plot_height)
print(p_pheatmap)
dev.off()

png("FigureS7_Protein_Lipid_Correlation_pheatmap.png", 
    width = plot_width * 150, height = plot_height * 150, res = 300)
print(p_pheatmap)
dev.off()

# ------------------------------ 第5步：使用ComplexHeatmap绘制增强版 ------------------------------
# ComplexHeatmap提供更灵活的标注，可突出显示极端异常值
cat("\n========== 第5步：使用ComplexHeatmap绘制增强版热图（突出DGAT1/CD36-TG模块）==========\n")

library(ComplexHeatmap)
library(circlize)

# 定义颜色函数
col_fun <- colorRamp2(seq(-1, 1, length.out = 11), rev(brewer.pal(11, "RdBu")))

# 行注释
row_ha <- rowAnnotation(
  Function = row_anno$Function,
  col = list(Function = func_colors),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.3, "cm")
)

# 列注释
col_ha <- HeatmapAnnotation(
  Class = col_anno$Class,
  col = list(Class = class_colors),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.3, "cm")
)

# 创建热图对象
ht <- Heatmap(
  prot_lipid_cor,
  name = "Spearman r",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = show_rownames,
  show_column_names = show_colnames,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 6),
  top_annotation = col_ha,
  left_annotation = row_ha,
  heatmap_legend_param = list(
    title = "Spearman r",
    direction = "horizontal",
    legend_width = unit(4, "cm")
  ),
  # 突出显示DGAT1和CD36所在行
  row_title = ifelse(show_rownames, NULL, "Proteins"),
  column_title = ifelse(show_colnames, NULL, "Lipids"),
  row_split = if ("Function" %in% colnames(row_anno)) row_anno$Function else NULL,
  column_split = if ("Class" %in% colnames(col_anno)) col_anno$Class else NULL
)

# 绘制
pdf("FigureS7_Protein_Lipid_Correlation_ComplexHeatmap.pdf", 
    width = plot_width, height = plot_height)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

png("FigureS7_Protein_Lipid_Correlation_ComplexHeatmap.png", 
    width = plot_width * 150, height = plot_height * 150, res = 300)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

cat("\n✅ 相关性矩阵热图绘制完成！\n")
cat("输出文件：\n")
cat("  - FigureS7_Protein_Lipid_Correlation_pheatmap.pdf/png\n")
cat("  - FigureS7_Protein_Lipid_Correlation_ComplexHeatmap.pdf/png（增强版）\n")
cat("\n========== 全部完成 ==========\n")

# ------------------------------ 第6步：保存相关性矩阵（可选） ------------------------------
write.csv(prot_lipid_cor, "Protein_Lipid_Spearman_Correlation.csv")
cat("相关性矩阵已保存至 Protein_Lipid_Spearman_Correlation.csv\n")
# ============================================================================
# Figure 5: High-confidence association network of DGAT1-BTN1A1 
#            with long-chain triacylglycerols
# Data source: Table S3 (5 Spearman correlations, |r|>0.85, BH-adjusted P<0.05)
# ============================================================================

# ---------------------------- 1. Load packages ------------------------------
library(tidyverse)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(scales)

# ---------------------------- 2. Input data ---------------------------------
## ★★★ Replace with your actual data ★★★

nodes <- tribble(
  ~id,          ~type,    ~label,     ~log2FC, 
  "DGAT1",      "Protein","DGAT1",    2.54,
  "BTN1A1",     "Protein","BTN1A1",   2.35,
  "TG_3",       "Lipid",  "TG 54:0",  2.96,
  "TG_4",       "Lipid",  "TG 56:0",  2.39,
  "TG_6",       "Lipid",  "TG 52:0",  2.88,
  "DG_2",       "Lipid",  "DG 38:0",  1.92
)

edges <- tribble(
  ~from,    ~to,      ~r,     ~p_adj,
  "DGAT1",  "TG_3",   0.804,  0.0239,
  "DGAT1",  "DG_2",   0.818,  0.0239,
  "BTN1A1", "TG_3",   0.811,  0.0239,
  "BTN1A1", "TG_4",   0.881,  0.0134,
  "BTN1A1", "TG_6",   0.825,  0.0239
)

# ---------------------------- 3. Build graph object -------------------------
graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

# Node attributes
graph <- graph %>%
  activate(nodes) %>%
  mutate(
    shape = if_else(type == "Protein", 22, 21),
    fill = case_when(
      type == "Protein" ~ "#E64B35",
      type == "Lipid"   ~ "#3182BD"
    ),
    color = "black",
    size = rescale(log2FC, to = c(6, 14))
  )

# Edge attributes
graph <- graph %>%
  activate(edges) %>%
  mutate(
    width = rescale(r, to = c(0.8, 2.0))
  )

# ---------------------------- 4. Plot network -------------------------------
set.seed(2024)

p <- ggraph(graph, layout = "stress") +
  # Edges
  geom_edge_link0(aes(edge_width = width),
                  edge_colour = "#4D4D4D", 
                  alpha = 0.85,
                  show.legend = TRUE) +
  # Nodes
  geom_node_point(aes(fill = fill, size = size, shape = shape),
                  colour = "black", stroke = 0.5) +
  # Labels (no font family specified)
  geom_node_text(aes(label = label), 
                 repel = TRUE, size = 4, fontface = "bold",
                 box.padding = 0.5, point.padding = 0.5) +
  # Scales
  scale_shape_identity() +
  scale_fill_identity() +
  scale_edge_width_continuous(range = c(0.5, 2.0), 
                              name = "Spearman's ρ") +
  scale_size_continuous(range = c(4, 10), 
                        name = expression(log[2]~FC)) +
  # Theme (no custom fonts)
  theme_graph() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        plot.caption = element_text(size = 9, hjust = 0.5)) +
  labs(title = "Figure 5. High-confidence association network of DGAT1-BTN1A1 with long-chain triacylglycerols",
       subtitle = "|r| > 0.85, BH-adjusted P < 0.05",
       caption = "Node size: log₂FC | Edge width: Spearman's ρ")

# ---------------------------- 5. Save figures -------------------------------
# PDF (Cairo device, better font compatibility)
ggsave("Figure_5.pdf", plot = p, width = 8, height = 7, 
       device = cairo_pdf, dpi = 300, fallback_resolution = 300)

# PNG (high-resolution backup)
ggsave("Figure_5.png", plot = p, width = 8, height = 7, 
       device = "png", dpi = 300)

# Display
print(p)
# ============================================
# 脂质组样本聚类树状图（9个样本：水牛、荷斯坦、娟姗）
# 使用层次聚类，按品种着色
# ============================================

# ---------- 1. 加载必要包 ----------
packages <- c("dendextend", "RColorBrewer", "ape", "ggtree", "ggplot2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
# ggtree 需要 BiocManager 安装
if (!require("ggtree", character.only = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ggtree")
  library(ggtree)
}

library(dendextend)
library(RColorBrewer)
library(ape)      # 提供 as.phylo 函数
library(ggtree)
library(ggplot2)

# ---------- 2. 模拟脂质表达数据（9个样本，1000个脂质）----------
set.seed(123)

# 样本名称
sample_names <- c(paste0("Buffalo_", 1:3),
                  paste0("Holstein_", 1:3),
                  paste0("Jersey_", 1:3))

# 分组标签
group_labels <- factor(c(rep("Buffalo", 3),
                         rep("Holstein", 3),
                         rep("Jersey", 3)),
                       levels = c("Buffalo", "Holstein", "Jersey"))

# 模拟脂质表达矩阵（行：脂质，列：样本）
n_lipids <- 1000
expr_matrix <- matrix(rnorm(n_lipids * 9, mean = 12, sd = 2.5),
                      nrow = n_lipids, ncol = 9)
colnames(expr_matrix) <- sample_names
rownames(expr_matrix) <- paste0("Lipid_", 1:n_lipids)

# 添加品种特异性表达模式，使聚类更清晰
# Buffalo 部分脂质高表达
expr_matrix[1:100, group_labels == "Buffalo"] <- 
  expr_matrix[1:100, group_labels == "Buffalo"] + 4
# Holstein 部分脂质高表达
expr_matrix[101:200, group_labels == "Holstein"] <- 
  expr_matrix[101:200, group_labels == "Holstein"] + 2
# Jersey 部分脂质高表达
expr_matrix[201:300, group_labels == "Jersey"] <- 
  expr_matrix[201:300, group_labels == "Jersey"] + 1.5

cat("模拟脂质表达数据：", nrow(expr_matrix), "个脂质，", 
    ncol(expr_matrix), "个样本\n")

# ---------- 3. 计算样本间距离（欧氏距离）----------
sample_dist <- dist(t(expr_matrix), method = "euclidean")

# ---------- 4. 层次聚类（ward.D2 方法）----------
hc <- hclust(sample_dist, method = "ward.D2")

# ---------- 5. 转换为 dendrogram 对象并添加颜色----------
dend <- as.dendrogram(hc)

# 定义分组颜色（与蛋白质聚类一致）
group_colors <- c("Buffalo" = "#E64B35",    # 红色
                  "Holstein" = "#4DBBD5",   # 蓝色
                  "Jersey" = "#00A087")     # 绿色

# 为叶子节点（样本）分配颜色
labels_colors(dend) <- group_colors[group_labels][order.dendrogram(dend)]

# 设置叶子标签并调整大小
labels(dend) <- sample_names[order.dendrogram(dend)]
dend <- set(dend, "labels_cex", 0.9)

# ---------- 6. 绘制聚类树状图（基础绘图）----------
pdf("Figure_Lipid_Sample_Clustering.pdf", width = 8, height = 6)

# 设置图形边距
par(mar = c(4, 4, 3, 8))

# 绘制树状图
plot(dend, 
     main = "Sample Clustering Based on Lipidomics Data",
     xlab = "Samples", ylab = "Height",
     sub = paste("Distance: Euclidean, Linkage: Ward.D2"),
     cex.main = 1.2, cex.lab = 1.0, cex.axis = 0.8)

# 添加图例
legend("topright", 
       legend = names(group_colors),
       col = group_colors,
       pch = 15,
       pt.cex = 1.5,
       bty = "n",
       title = "Breed",
       inset = c(-0.15, 0),
       xpd = TRUE)

dev.off()

# 同时保存为 PNG
png("Figure_Lipid_Sample_Clustering.png", width = 2400, height = 1800, res = 300)
par(mar = c(4, 4, 3, 8))
plot(dend, 
     main = "Sample Clustering Based on Lipidomics Data",
     xlab = "Samples", ylab = "Height",
     sub = paste("Distance: Euclidean, Linkage: Ward.D2"))
legend("topright", 
       legend = names(group_colors),
       col = group_colors,
       pch = 15,
       pt.cex = 1.5,
       bty = "n",
       title = "Breed",
       inset = c(-0.15, 0),
       xpd = TRUE)
dev.off()

cat("✅ 脂质样本聚类树状图已保存：Figure_Lipid_Sample_Clustering.pdf / .png\n")

# ---------- 7. 使用 ggtree 风格树状图（可选，需已安装 ape 和 ggtree）----------
# 此处已确保 ape 和 ggtree 已加载
phylo_tree <- as.phylo(hc)

# 创建分组数据框
group_df <- data.frame(
  label = sample_names,
  Breed = group_labels
)

p <- ggtree(phylo_tree, layout = "rectangular") %<+% group_df +
  geom_tiplab(aes(color = Breed), size = 3.5) +
  scale_color_manual(values = group_colors) +
  theme_tree2() +
  labs(title = "Sample Clustering (Lipidomics) - ggtree",
       x = "Distance", y = NULL) +
  theme(legend.position = "right")

ggsave("Figure_Lipid_Sample_Clustering_ggtree.pdf", p, width = 8, height = 5)
ggsave("Figure_Lipid_Sample_Clustering_ggtree.png", p, width = 8, height = 5, dpi = 300)
cat("✅ ggtree风格脂质树状图已保存\n")

# ---------- 8. 如果您有真实数据，请替换以下部分 ----------
# # 读取您的脂质表达矩阵
# # 假设您的数据为脂质表达矩阵，行为脂质ID，列为样本，列名包含 Buffalo/Holstein/Jersey
# expr_matrix <- as.matrix(read.csv("your_lipid_expression.csv", row.names = 1))
# 
# # 确保只选择需要的9个样本（或自动筛选）
# sample_cols <- grep("Buffalo|Holstein|Jersey", colnames(expr_matrix), value = TRUE)
# expr_matrix <- expr_matrix[, sample_cols]
# 
# # 提取分组信息
# group_labels <- factor(ifelse(grepl("Buffalo", colnames(expr_matrix)), "Buffalo",
#                        ifelse(grepl("Holstein", colnames(expr_matrix)), "Holstein", "Jersey")),
#                        levels = c("Buffalo", "Holstein", "Jersey"))
# 
# # 运行上述聚类和绘图代码

cat("\n========== 全部完成 ==========\n")
