library(ggplot2)
library(dplyr)
library(ggsignif)
library(stats)
options(warn = -1)
# 1. 数据导入与预处理
df <- read.csv("oral_ratio_data.csv", header = TRUE, encoding = "UTF-8")
df_clean <- df %>% filter(Oral_ratio != 0)
df_clean$Species <- factor(df_clean$Species, levels = c("GP", "RP", "MD", "RB"))
# 2. Kruskal-Wallis 全局检验
kw_res <- kruskal.test(Oral_ratio ~ Species, data = df_clean)
cat("==== Kruskal-Wallis 全局检验结果 ====\n")
print(kw_res)
# 3. 定义全部6组两两对比
all_comp <- list(
  c("GP", "RP"),
  c("GP", "MD"),
  c("GP", "RB"),
  c("RP", "MD"),
  c("RP", "RB"),
  c("MD", "RB"))
comp_names <- c("GP vs RP","GP vs MD","GP vs RB","RP vs MD","RP vs RB","MD vs RB")
# 4. 批量计算Wilcoxon检验 + FDR校正P值，输出完整P值表格（附表用）
p_raw <- c()
for (pair in all_comp) {
  sub1 = df_clean$Oral_ratio[df_clean$Species == pair[1]]
  sub2 = df_clean$Oral_ratio[df_clean$Species == pair[2]]
  wt = wilcox.test(sub1, sub2)
  p_raw = c(p_raw, wt$p.value)}
p_fdr <- p.adjust(p_raw, method = "fdr")
# 构建P值结果表
p_table <- data.frame(
  Comparison = comp_names,
  Raw_P = p_raw,
  FDR_adjust_P = p_fdr)
cat("\n==== 两两比较原始P值 & FDR校正P值（补充表）====\n")
print(p_table, row.names = F)
# 5. 筛选显著组别（P < 0.05） + 手动精准赋值星号
# 对应顺序：GP-RP(ns), GP-MD(****), GP-RB(****), RP-MD(****), RP-RB(****), MD-RB(ns)
sig_comp <- all_comp[c(2,3,4,5)]
sig_star <- c("****","****","****","****")
# 对应高度（只保留4条显著线的高度）
max_y <- max(df_clean$Oral_ratio)
y_pos_sig <- c(
  max_y * 1.13,
  max_y * 1.08,
  max_y * 1.03,
  max_y * 0.98)
# 6. 绘图：只标****，无任何P值、无显著组不显示
p <- ggplot(df_clean, aes(x = Species, y = Oral_ratio, fill = Species)) +
  geom_boxplot(
    width = 0.6,
    color = "black",
    size = 0.5,
    outlier.shape = NA,
    alpha = 0.8) +
  scale_fill_manual(values = c(
    "GP" = "#F1AB78",
    "RP" = "#80C7A7",
    "MD" = "#CBCFE9",
    "RB" = "#C884B6")) +
  scale_y_continuous(labels = function(x) x*100, expand = c(0.05,0.05)) +
  geom_signif(
    comparisons = sig_comp,
    annotations = sig_star,
    y_position = y_pos_sig,
    tip_length = 0.01,
    size = 0.5,
    textsize = 8,
    vjust = 0.6,
    color = "black") +
  labs(x = "", y = "Oral Microbe Source Ratio (%)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line("black", 0.5),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 13, face = "bold"),
    legend.position = "none",
    plot.margin = margin(10,10,10,10))
print(p)
# 保存图片
ggsave(
  "oral_ratio_boxplot_final_staronly.pdf",
  plot = p, width = 6, height = 5, device = "pdf", dpi = 300)
# 导出完整P值附表
write.csv(p_table, "pairwise_pvalue_FDR.csv", row.names = F, fileEncoding = "UTF-8")
options(warn = 0)
