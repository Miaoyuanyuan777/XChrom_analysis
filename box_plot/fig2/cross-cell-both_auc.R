library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
df <- read_csv('cross-cell-both_auc.csv')
# 将 Metric 转换为因子并指定顺序
df$Metric <- factor(df$Metric, levels = c("auROC", "per cell auROC", "per peak auROC", 
                                          "auPRC", "per cell auPRC", "per peak auPRC"))
df$Dataset<-factor(df$Dataset, levels = c("h_pbmc","m_brain","h_brain",'h_gonads','m_palates'))
df$Test<-factor(df$Test,levels=c('cross-cell','cross-both'))

##########------------------------------------------############################
# 添加一个新的列来标识
df <- df %>%
  mutate(MetricType = ifelse(grepl("auPRC", Metric), "auPRC", "auROC"))%>%
  mutate(MetricType = factor(MetricType, levels = c("auROC", "auPRC")))

# 统计每个分面类型的数量
metric_counts <- df %>%
  group_by(MetricType) %>%
  summarise(count = n(), .groups = "drop")

##########------------------------------------------############################
y_pos <- df %>% 
  group_by(Metric) %>% 
  summarise(max_val = max(Value))

annotation_df <- df %>%
  pivot_wider(names_from = Test, values_from = Value) %>%
  mutate(diff =  `cross-cell` - `cross-both`) %>%
  group_by(Metric) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    p_value = t.test(`cross-cell`, `cross-both`, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_value_label=ifelse(p_value < 0.05,
                              paste0("p = ", signif(p_value, digits = 3)
                              ),NA_character_))%>%
  mutate(MetricType = factor(case_when(
    grepl("auROC", Metric) ~ "auROC",
    grepl("auPRC", Metric) ~ "auPRC"),
    levels = c("auROC", "auPRC")))
# 合并标注数据
annotation_df <- left_join(annotation_df, y_pos, by = "Metric")
annotation_df$y_position <- ifelse(annotation_df$p_value < 0.05,
                                   annotation_df$max_val + 0.04, 
                                   NA_real_) 

# 绘制分面图，使用不同的纵坐标尺度并去除没有数据的分面轴
ggplot(df, aes(x = Metric, y = Value)) +
  geom_boxplot(aes(group = interaction(Metric, Test), fill = Test),
               # position = position_dodge(width = 0.8), 
               alpha = 1,outlier.shape = NA) +  # 绘制箱线图
  geom_jitter(
    aes(group = interaction(Metric, Test), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.3,   # 控制抖动幅度
      dodge.width = 0.6     # 控制不同 Scenario 的间距
    ),size = 3,alpha = 1,color = "black", #stroke = 1
    )+
  scale_fill_manual(values = c("cross-cell" = "#FF4500",'cross-both'='#FFB07C')) +  # 自定义箱线图颜色
  scale_shape_manual(
    values = c(
      "h_pbmc" = 1,    # 空心圆
      "m_brain" = 2,    # 空心三角
      "h_brain" = 5,    # 空心菱形
      "h_gonads" = 0,    # 空心方框
      "m_palates" = 6    # 空心倒三角
    ))+
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  # 增加 x 轴 Metric 之间的间距
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +  # 图例顺序
  facet_wrap(~MetricType, scales = "free")   +
  theme_set(theme_bw())+
  # theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16, face = "bold", color = "black", margin = margin(t =2)),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold",margin = margin(b = 10)),
    legend.text = element_text(size = 16),
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3),color = "black", size = 0.5),  # 设置分面标签背景
    strip.text = element_text(size = 16, face = "bold", color = "black"),  # 设置分面标签的文字样式
    strip.text.x = element_text(size = 16, face = "bold",color = "black"),  # 设置分面标签的文本样式
    # 控制主面板边框粗细
    panel.border = element_rect(size = 0.5, color = "black"),  # 边框线条粗细设为 0.5
  ) +
  labs(title = "Cross-both/cell Prediction", x = "", y = "Value")+ 
  geom_text(
    data = annotation_df[annotation_df$p_value<0.05,],
    aes(x = Metric, y = y_position,group = MetricType,  # x=1.5 将注释居中于分面面板
        label = paste0(p_value_label, "\n95%CI = [", 
                       round(lower_ci, 2), ",", 
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4,
    color = "black",hjust = 0.5)
ggsave("cross-both_auc_box2.pdf", plot = last_plot(), 
       width = 10, height = 6, units = "in",dpi=300)
