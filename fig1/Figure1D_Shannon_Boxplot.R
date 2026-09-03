# 安装并加载所需包
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("readxl")) install.packages("readxl")
if (!require("ggsignif")) install.packages("ggsignif")
library(ggplot2)
library(readxl)
library(ggsignif)
options(warn = -1)
# ===================== 1. 数据读取与分组设置 =====================
df <- read_excel("D:/999-sj/xin/2-1.Alpha-index/alpha1.xlsx")
# X轴分组固定顺序
group_order <- c("GPO", "RPO", "MDO", "RBO",  # 口腔4组
                 "GPF", "RPF", "MDF", "RBF") # 粪便4组
# 自定义填充配色
group_colors <- c(
  "GPO" = "#F1AB78", "RPO" = "#80C7A7", "MDO" = "#CBCFE9", "RBO" = "#C884B6",
  "GPF" = "#3D73B9", "RPF" = "#F3A8A8", "MDF" = "#C7E3F6", "RBF" = "#E1C4A6")
df$`Sample ID` <- factor(df$`Sample ID`, levels = group_order)
# ===================== 2. Kruskal-Wallis 全局多组检验 =====================
kw_result <- kruskal.test(shannon_entropy ~ `Sample ID`, data = df)
cat("================ Kruskal-Wallis 全局检验结果 ================\n")
print(kw_result)
# ===================== 3. 全部12组两两对比列表 =====================
all_compare_list <- list(
  # 口腔内部6组
  c("GPO", "RPO"),
  c("GPO", "MDO"),
  c("GPO", "RBO"),
  c("RPO", "MDO"),
  c("RPO", "RBO"),
  c("MDO", "RBO"), # FDR P>0.05，无显著，绘图剔除
  # 粪便内部6组
  c("GPF", "RPF"),
  c("GPF", "MDF"),
  c("GPF", "RBF"),
  c("RPF", "MDF"),
  c("RPF", "RBF"),
  c("MDF", "RBF"))
comp_names <- sapply(all_compare_list, function(x) paste0(x[1], " vs ", x[2]))
# ===================== 4. 批量计算原始P、FDR校正P，导出完整附表 =====================
raw_P <- c()
for (pair in all_compare_list) {
  grp1 <- df$shannon_entropy[df$`Sample ID` == pair[1]]
  grp2 <- df$shannon_entropy[df$`Sample ID` == pair[2]]
  wt_res <- wilcox.test(grp1, grp2)
  raw_P <- c(raw_P, wt_res$p.value)}
fdr_P <- p.adjust(raw_P, method = "fdr")
# 构建12组完整P值表格
p_table <- data.frame(
  Comparison = comp_names,
  Raw_Pvalue = raw_P,
  FDR_adjust_Pvalue = fdr_P)
cat("\n================ 12组两两比较P值汇总表 ================\n")
print(p_table, row.names = FALSE)
# 导出CSV补充表
write.csv(
  p_table,
  file = "D:/999-sj/xin/2-1.Alpha-index/Shannon_12pairwise_P_table.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8")
cat("\n完整P值表格已导出至目标文件夹\n")
# ===================== 5. 筛选显著对比（FDR P < 0.05） =====================
sig_index <- which(fdr_P < 0.05)
sig_compare_list <- all_compare_list[sig_index]
# 11条显著性横线高度（口腔5条 + 粪便6条，按需微调数值）
sig_y_height <- c(
  12.0, 11.6, 11.2, 10.8, 10.4, # 口腔5组
  12.0, 11.6, 11.2, 10.8, 10.4, 10.0 # 粪便6组)
# 固定对应星号（严格匹配你的P值结果，无任何P数字）
sig_star_labels <- c(
  "****", # GPO vs RPO
  "****", # GPO vs MDO
  "****", # GPO vs RBO
  "****", # RPO vs MDO
  "****", # RPO vs RBO
  "**",   # GPF vs RPF
  "****", # GPF vs MDF
  "****", # GPF vs RBF
  "****", # RPF vs MDF
  "****", # RPF vs RBF
  "****"  # MDF vs RBF)
# ===================== 6. 绘图：仅显示星号，完全不出现P值 =====================
p <- ggplot(df, aes(x = `Sample ID`, y = shannon_entropy, fill = `Sample ID`)) +
  geom_boxplot(width = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = group_colors) +
  labs(x = "", y = "shannon_entropy") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    legend.position = "none") +
  geom_signif(
    comparisons = sig_compare_list,
    annotations = sig_star_labels, # 强制仅展示星号，屏蔽P值
    y_position = sig_y_height,
    tip_length = 0.01,
    size = 0.5,
    textsize = 15,
    vjust = 0.8,
    color = "black")
# 预览画布
print(p)
# ===================== 7. 导出高清PDF主图 =====================
ggsave(
  filename = "Shannon_Boxplot_Only_Signif_Star.pdf",
  plot = p,
  path = "D:/999-sj/xin/2-1.Alpha-index",
  width = 9.5,
  height = 8,
  device = "pdf",
  dpi = 300)
options(warn = 0)
