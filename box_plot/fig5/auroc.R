library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
df <- read_csv("auroc.csv")

## human物种的human model中有一个离群点，
# 在geom_boxplot作箱型图时指定outlier.shape = NA不画它，否则在画geom抖动点时会重复

# make model factor
df$Model <- factor(df$Model, levels = c("human_model", "mouse_model"))
df$Species<-factor(df$Species,levels=c('human','mouse','macaque','marmoset'))

##############----------auROC--------#########
filtered_df <- df %>% filter(Metric == "per cell auROC")  
## per cell auROC
## per peak auROC


#########----------filt human & mouse which are not paired-------##########
filtered_df <- filtered_df %>% 
  group_by(Species, Samples) %>%
  filter(n() == 2)
################################################################################
# 1. 为了标注，我们确定每个 species 中数值的上限，用于放置水平线和文字
y_pos <- filtered_df %>% 
  group_by(Species) %>% 
  summarise(max_val = max(Value)) %>% 
  mutate(y_position = max_val + 0.02)  # 根据实际数据调整偏移量

annotation_df <- filtered_df %>%
  pivot_wider(names_from = Model, values_from = Value) %>%
  mutate(diff =  mouse_model - human_model) %>%
  group_by(Species) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    p_value = wilcox.test(human_model, mouse_model, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_value_label = sprintf("p = %.3f", p_value)
  )

# 合并标注数据
annotation_df <- left_join(annotation_df, y_pos, by = "Species")

# # 3. 绘制 boxplot，并添加水平线和标注
# ggplot(filtered_df, aes(x = Species, y = Value,group = interaction(Species,Model),fill=Model)) +
#   geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.9,outlier.shape = NA)+ #不画离群点，否则会和抖动点重复
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
#               size = 2) +
#   # 手动设置填充色
#   scale_fill_manual(values = c("human_model" = "orange", "mouse_model" = "#7B4F94")) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black"),
#     axis.title = element_text(size = 12, color = "black", face = "bold"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
#   ) +
#   labs(title = "", x = "", y = "auROC") +
#   geom_text(data = annotation_df[annotation_df$p_value<0.05,],inherit.aes = FALSE,
#             aes(x = Species, y = y_position, 
#                 label = paste0(p_value_label, "\n95%CI = [", round(lower_ci,2), ",", round(upper_ci,2), "]")),
#             size = 3, color = "black")


# 绘图
ggplot(filtered_df, aes(x = Model, y = Value, fill = Model)) +
  # 箱线图
  geom_boxplot(width = 0.6, alpha = 0.9, outlier.shape = NA) + 
  # 散点图（抖动点）
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) + 
  # 分面（按Species分组）
  facet_wrap(~ Species, scales = "free_x", strip.position = "bottom",nrow=1) + 
  # 颜色设置
  scale_fill_manual(values = c("human_model" = "orange", 
                               "mouse_model" = "#7B4F94")) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  # 隐藏x轴标签
    axis.ticks.x = element_blank(), # 隐藏x轴刻度
    axis.text.y = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold", color = "black"),  # 分面标签样式
    strip.background = element_blank(),  # 分面标签背景透明
    panel.spacing = unit(0.5, "lines")  # 增加分面之间的间距
  ) +
  labs(title = "", x = "", y = "per cell auROC") +
  geom_text(
    data = annotation_df[annotation_df$p_value<0.05,],
    aes(x = 1.5, y = y_position,  # x=1.5 将注释居中于分面面板
        label = paste0(p_value_label, "\n95%CI = [", 
                       round(lower_ci, 2), ",", 
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4,
    color = "black",hjust = 0.5  # 文字水平居中
  )


ggsave("percellauROC1.pdf", plot = last_plot(), 
       width = 9, height = 4, units = "in",dpi=300)
