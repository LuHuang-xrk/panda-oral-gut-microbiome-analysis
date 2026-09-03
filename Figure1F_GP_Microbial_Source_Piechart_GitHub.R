# 加载必要的包
library(ggplot2)  # 绘图核心包
library(dplyr)    # 数据处理包
# 设置文件路径（建议：用file.path()增强跨系统兼容性）
gpo_file <- file.path("sink_predictions_GPO_contributions.txt")
unknown_file <- file.path("GP sink_predictions_Unknown_contributions.txt")
output_path <- file.path("GP_pie_chart.pdf")
# 读取数据文件（优化：增加错误处理，避免文件不存在时报错）
if (!file.exists(gpo_file)) stop(paste("文件不存在：", gpo_file))
if (!file.exists(unknown_file)) stop(paste("文件不存在：", unknown_file))
gpo_data <- read.delim(gpo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
unknown_data <- read.delim(unknown_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# 计算两类的总和（核心逻辑：提取所有数值列并求和，忽略NA）
gpo_total <- sum(unlist(gpo_data[sapply(gpo_data, is.numeric)]), na.rm = TRUE)
unknown_total <- sum(unlist(unknown_data[sapply(unknown_data, is.numeric)]), na.rm = TRUE)
# 构建饼图数据框
pie_data <- data.frame(
  category = c("oral", "unknown"),
  value = c(gpo_total, unknown_total),
  stringsAsFactors = FALSE  # 优化：显式关闭因子转换，避免后续绘图警告)
# 计算占比（保留2位小数，仅用于打印）
pie_data <- pie_data %>%
  mutate(
    percentage = round(value / sum(value) * 100, 2))
# 绘制饼图
p <- ggplot(pie_data, aes(x = "", y = value, fill = category)) +
  geom_col(width = 1, color = NA) +
  coord_polar("y", start = 0) +
  # ===================== 自定义配色 =====================
scale_fill_manual(
  values = c("oral" = "#C7E3F6", "unknown" = "#CBCFE9"), 
  labels = c("oral", "unknown")) +
  # ======================================================
labs(title = "GP", fill = "") +
  theme_minimal() +
  theme(
    # 隐藏坐标轴
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # 标题样式（已放大）
    plot.title = element_text(
      hjust = 0.5,        
      size = 30,          
      face = "bold",      
      margin = margin(b = -20)),
    # ===================== 调整：图例远离饼图 =====================
    legend.position = c(1, 0.9),  # 数值变大 = 向右移动，远离饼图
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.3, "cm"),
    legend.text = element_text(size = 27),
    legend.title = element_blank(),
    # ============================================================
    # 去除网格线和背景
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10, unit = "mm")  # 右边距加大，防止图例被裁切)
# ================ 在 RStudio 中显示图表 ================
print(p)  
# ========================================================
# 保存饼图为PDF
ggsave(
  filename = output_path,
  plot = p,
  width = 8,
  height = 8,
  device = "pdf",
  dpi = 300)
# 打印占比信息
cat(sprintf("口腔(oral)占比：%.2f%%\n", pie_data$percentage[1]))
cat(sprintf("未知(unknown)占比：%.2f%%\n", pie_data$percentage[2]))
