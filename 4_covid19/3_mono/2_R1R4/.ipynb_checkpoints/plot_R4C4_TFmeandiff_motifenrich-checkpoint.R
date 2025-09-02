## motif富集用bubble图，mean diff用bar,颜色映射统一是显著性
library(scales)  # 用于归一化函数
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(tidyverse)
library(cowplot)
library(scales)

r4 <- read_csv("./R4_other_diff/DEGresults.csv")
motif_r4 <- read_csv("./R4_other_diff/C4UP_diff005p_enrichmotif.csv")


# # 筛选显著结果
r4_sig <- r4 %>%
  filter(pval_adj < 0.05, abs(mean_difference) > 0.02)

merge_data <- function(rna_data, motif_data, cluster){
  inner_join(rna_data %>% 
               mutate(cluster = cluster),motif_data %>% 
               rename(tf = TF_prefix),by = "tf") %>%
    mutate(sig_stars = cut(mlog10padj,
                           breaks = c(-Inf, 1.3, Inf),
                           labels = c(' ',"*")))}

df_r4 <- merge_data(r4_sig, motif_r4, "R4")

## ##################################
# ----这里是根据R1/R4的mean diff 排序------
## ##################################
r4_diffsort <- df_r4 %>%
  # filter(tf %in% common_tfs) %>%
  group_by(tf) %>%
  arrange(desc(mean_difference)) %>% 
  pull(tf)

r41 <- df_r4 %>% 
  # filter(tf %in% common_tfs) %>% 
  mutate(tf = factor(tf, levels = r4_diffsort)) %>% 
  mutate(mlog10_pval_adj = -log10(pval_adj + 1e-300))

# 创建图形 --------------------------------------------------------
# (1) 上方均值差异条形图
p_mean_diff <- ggplot(r41, aes(x = tf)) +
  geom_col(aes(y = mean_difference, fill = mlog10_pval_adj),
           width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid",
             color = "grey", linewidth = 0.5) +
  scale_fill_gradient(
    low = "grey90", high = "red",
    limits = c(0, ceiling(max(r41$mlog10_pval_adj))),
    name = "-log10(Adjusted-Pval)",
    guide = guide_colorbar(      # 控制图例的核心参数
      direction = "horizontal",  # 水平方向
      title.position = "top",    # 标题位置
      barwidth = unit(4, "cm")   # 色条长度
    )) +
  labs(y = "Activity mean difference\n(R4 vs Other)", title = "COVID-19:R4(cMono)") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )
p_mean_diff

# (2) 下方motif富集气泡图

# 创建新列处理颜色映射逻辑（保留原数据完整性）
r41 <- r41 %>%
  mutate(fill_value = if_else(mlog10padj >= 1.3, mlog10padj, 0))

# 绘图代码
p_motif <- ggplot(r41, aes(x = tf)) +
  geom_point(aes(y = 0, fill = fill_value),
             size = 7, 
             shape = 21,
             alpha = 0.8) +
  # geom_text(aes(y = 0, label = sig_stars),
  #           size = 7, 
  #           color = "red",
  #           vjust = 0.5,   # 垂直居中
  #           hjust = 0.5) + # 水平居中
  scale_fill_gradientn(
    # 颜色设置：前两个颜色为灰色，后续依次为蓝色到黑色
    colours = c("grey80", "blue", "black"),
    # values 参数根据 limits 对照，1.3 对应归一化后的比例为 1.3/max_val
    values = scales::rescale(c(0, 1.3, ceiling(max(r41$mlog10padj)))),
    # breaks 只显示三个刻度：0、1.3 和最大值
    breaks = c(0,1.3,ceiling(max(r41$mlog10padj))),
    # 明确设定填充值的取值范围
    limits = c(0, ceiling(max(r41$mlog10padj))),
    name = "-log10(Adjusted-Pval)",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(4, "cm"))
  ) +
  labs(x = "Transcription factor (TF)") +
  scale_y_continuous(
    limits = c(-0.001, 0.001),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1.2,
      margin = margin(t = -6)
    ),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(-2,0,0,0), "cm"),
    aspect.ratio = 0.1,
    legend.position = "right"
  )

# 绘制图形
p_motif


# 组合图形 --------------------------------------------------------
aligned_plots <- align_plots(p_mean_diff, p_motif, align = "v")
final_plot <- plot_grid(
  aligned_plots[[1]], 
  aligned_plots[[2]],
  ncol = 1, 
  rel_heights = c(0.2, 0.4)  # 调整上下部分高度比例
)

final_plot

ggsave("R4C4_bubble-bar2-1.pdf", final_plot, 
       width = 16, height = 8, dpi = 300)
