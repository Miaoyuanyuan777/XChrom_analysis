# /1_all_mild_severe_rna/tf_activity_mono_RevSigmoid/cMono_sm_out/cMono_sm-diffout_add-cMono-corr.csv
# cMono中severe vs mild的mean diff以及cMono中的activity与表达corr

# 加载库
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

set.seed(123)
results_df <- read_csv("cMono_sm-diffout_add-cMono-corr.csv")

# 数据准备
results_df <- results_df %>%
  mutate(
    significant = (pval_adj < 0.05) & (abs(mean_difference) > 0.006),
    neg_log_padj = `-log10(padj)`,
    neg_log_corr_padj = -log10(corr_padj+1e-300),
    Correlation = as.numeric(Correlation)
  )
results_df <- results_df %>%
  mutate(
    corrp_size = ifelse(neg_log_corr_padj < 1.3,1.3,ifelse(neg_log_corr_padj > 50,50,neg_log_corr_padj))
  )

# 颜色映射范围设置
vmax <- max(results_df$Correlation, na.rm = TRUE)
vmin <- min(results_df$Correlation, na.rm = TRUE)

# 主绘图
p <- ggplot(results_df, aes(x = mean_difference, y = neg_log_padj)) +
  geom_point(
    aes(
      color = Correlation,
      size = corrp_size
    ),
    alpha = 0.7
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0,
    limits = c(vmin, vmax),
    name = "TF activity_expression\ncorrelation(cMono)"
  )+
  # scale_color_gradientn(
  #   colours = viridisLite::viridis(256, #option = "rocket",
  #                                  begin = 0.1,   # 调整起始位置，增强低端颜色
  #                                  end = 0.9),
  #   limits = c(vmin, vmax),
  #   name = "TF activity_expression\ncorrelation(cMono)"
  # ) +
  scale_size_continuous(
    range = c(0.5,1.5),  # 对应Python的sizes参数
    breaks = c(1.3, 10, 50),  # 自定义断点
    name = "Correlation:-log10(Adjusted-Pval)"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    x = "Activity mean difference (cMono:severe vs mild)",
    y = "-log10(Adjusted-Pval)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(linewidth = 0.2),
    axis.line = element_line(linewidth = 0.5)
  )
p
# 添加文本标签
significant_data <- subset(results_df, significant)
p1 <- p + geom_text_repel(
  data = significant_data,
  aes(label = TF,color = Correlation),
  size = 3,
  box.padding = 0.5,
  point.padding = 0.4,
  min.segment.length = 0.3,
  segment.color = "grey50",
  segment.size = 0.05,
  direction = "both",
  max.overlaps = Inf,
  force = 30,
  fontface = "bold",       # 加粗字体
  bg.color = "white",      # 设置文本背景颜色
  bg.r = 0.15 
)
p1
# 坐标轴格式调整
p2 <- p1 + scale_x_continuous(
  trans = scales::pseudo_log_trans(base = 3, sigma = 0.001),
  breaks = scales::pretty_breaks(n = 3)
)


p2
# 输出图形
ggsave(
  filename = file.path("MA_cMono-sm_addCorr_R2.pdf"),
  plot = p2,
  width = 10,
  height = 4,
  dpi = 300
)
