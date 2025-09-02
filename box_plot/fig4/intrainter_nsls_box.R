library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

df <- read_csv("intrainter_nsls.csv")

# 将 Metric 转换为因子并指定顺序
df$Metric <- factor(df$Metric, levels = c("ns(10)", "ns(50)", "ns(100)", 
                                          "ls(10)", "ls(50)", "ls(100)"))
df$Scenario<-factor(df$Scenario,levels=c('raw_atac','inter','intra'))
df$Dataset<-factor(df$Dataset,levels = c('s1d2','s2d1','s3d10','s4d1','s2d4'))


# 绘制分面图，使用不同的纵坐标尺度并去除没有数据的分面轴
df <- df %>%
  mutate(MetricType = ifelse(grepl("ns", Metric), "neighbor score", "label score"))%>%
  mutate(MetricType = factor(MetricType, levels = c("neighbor score", "label score")))

# # 进行配对 t 检验,只计算了intra和inter的配对t，没有算raw atac的
# t_test_results <- df %>%
#   spread(Scenario, Value) %>%
#   group_by(Metric) %>%
#   summarise(
#     t_test_p_value = t.test(intra, inter, paired = TRUE)$p.value,
#     max_value = max(c(max(intra, na.rm = TRUE), max(inter, na.rm = TRUE)))  
#   )
# # 创建一个新的数据框用于标注位置
# p_value_labels <- t_test_results %>%
#   mutate(label = paste0("p = ", signif(t_test_p_value, digits = 3))) %>%
#   mutate(y_position = max_value + 0.05)  # 设置 p 值标注的 y 坐标位置（标注位置略高于最大值）

y_pos <- df %>% 
  group_by(Metric) %>% 
  summarise(max_val = max(Value))

annotation_df <- df %>%
  pivot_wider(names_from = Scenario, values_from = Value) %>%
  mutate(diff =  intra - inter) %>%
  group_by(Metric) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    p_value = t.test(intra, inter, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_value_label=ifelse(p_value < 0.05,
                              paste0("p = ", signif(p_value, digits = 3)
                              ),NA_character_))%>%
  mutate(MetricType = factor(case_when(
    grepl("ls", Metric) ~ "label score",
    grepl("ns", Metric) ~ "neighbor score"),
    levels = c("neighbor score", "label score")))

# 合并标注数据
annotation_df <- left_join(annotation_df, y_pos, by = "Metric")
annotation_df$y_position <- ifelse(annotation_df$p_value < 0.05,
                                   annotation_df$max_val + 0.04, 
                                   NA_real_) 

ggplot(df, aes(x = Metric, y = Value,)) +
  geom_boxplot(aes(group = interaction(Metric, Scenario), fill = Scenario), 
               alpha = 1, outlier.shape = NA) +  # 绘制箱线图
  geom_jitter(
    aes(group = interaction(Metric, Scenario), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.3,   # 控制抖动幅度
      dodge.width = 0.6     # 控制不同 Scenario 的间距
    ),size = 3,alpha = 0.8,color = "black",stroke = 1)+
  scale_fill_manual(values = c('raw_atac'='grey',"intra" = "#aed09c", "inter" = "#FFD700")) +  # 自定义箱线图颜色
  # scale_color_manual(values = c("s1d2" = "#2ca02c", "s2d1" = "#d62728", "s2d4" = "purple", 
  #                               "s4d1" = "#8c564b", "s3d10" = "#FF6699")) +  # 自定义点的颜色
  scale_shape_manual(
    values = c(
      "s1d2" = 1,    # 空心圆
      "s2d1" = 2,    # 空心三角
      "s2d4" = 5,    # 空心菱形
      "s4d1" = 0,    # 空心方框
      "s3d10" = 6    # 空心倒三角
    ))+
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  # 增加 x 轴 Metric 之间的间距
  guides(shape = guide_legend(order = 1), fill = guide_legend(order = 2)) +  # 图例顺序
  facet_wrap(~MetricType, scales = "free")+
  theme_set(theme_bw())+
  # theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black", margin = margin(t =2)),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold",margin = margin(b = 10)),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3),color = "black", size = 0.5),  # 设置分面标签背景
    strip.text = element_text(size = 12, face = "bold", color = "black"),  # 设置分面标签的文字样式
    strip.text.x = element_text(size = 12, face = "bold",color = "black"),  # 设置分面标签的文本样式
  ) +
  labs(title = "Metrics Across Scenario", x = "", y = "Value") +
  geom_text(
    data = annotation_df[annotation_df$p_value<0.05,],
    aes(x = Metric, y = y_position,group = MetricType,  # x=1.5 将注释居中于分面面板
        label = paste0(p_value_label, "\n95%CI = [", 
                       round(lower_ci, 2), ",", 
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4,
    color = "black",hjust = 0.5)

ggsave("intrainter_nsls_box1.pdf", plot = last_plot(), 
       width = 10, height = 6, units = "in",dpi=300)
