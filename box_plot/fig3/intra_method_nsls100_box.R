library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
df <- read_csv('intra_method_nsls100.csv')

# 将 Metric 转换为因子并指定顺序
df$Metric <- factor(df$Metric, levels = c("ns(10)", "ns(50)", "ns(100)", 
                                          "ls(10)", "ls(50)", "ls(100)"))
df$Method <- factor(df$Method, levels = c('raw_atac',"Seurat", "LS_Lab", "MultiVI", "XChrom"))
df$Dataset<-factor(df$Dataset, levels = c("h_pbmc","m_brain","h_brain",'h_gonads','m_palates'))

##########------------------------------------------############################
# 添加一个新的列来标识 ns 和 ls
df <- df %>%
  mutate(MetricType = ifelse(grepl("ns", Metric), "neighbor score", "label score"))%>%
  mutate(MetricType = factor(MetricType, levels = c("neighbor score", "label score")))

# 统计每个分面类型的数量
metric_counts <- df %>%
  group_by(MetricType) %>%
  summarise(count = n(), .groups = "drop")

#################################### 分面 ######################################
y_pos <- df %>% 
  group_by(Metric) %>% 
  summarise(max_val = max(Value))

# 获取所有需要比较的方法（排除XChrom自身）
other_methods <- unique(df$Method[df$Method != "XChrom"])

# 提取XChrom的数据作为基准
xchrom_data <- df %>% 
  filter(Method == "XChrom") %>% 
  select(Dataset, Metric, MetricType, Value_XChrom = Value)

# 生成比较数据
comparison_data <- df %>%
  filter(Method != "XChrom") %>% 
  inner_join(xchrom_data, by = c("Dataset", "Metric", "MetricType")) %>%
  mutate(diff = Value_XChrom - Value)

# 计算统计指标
annotation_df <- comparison_data %>%
  group_by(Metric, Method) %>%  # 按Metric和对比方法分组
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    p_value = t.test(Value_XChrom, Value, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_value_label = ifelse(p_value < 0.05,
                           paste0("p = ", signif(p_value, digits = 3)),
                           NA_character_)
  ) %>%
  # 合并MetricType和y轴位置信息
  left_join(distinct(df, Metric, MetricType), by = "Metric") %>% 
  left_join(
    df %>% group_by(Metric) %>% summarise(max_val = max(Value)),
    by = "Metric"
  ) %>% 
  mutate(y_position = ifelse(p_value < 0.05, max_val + 0.1, NA_real_))

# 绘制分面图，使用不同的纵坐标尺度并去除没有数据的分面轴
ggplot(df, aes(x = Metric, y = Value,)) +
  geom_boxplot(aes(group = interaction(Metric, Method), fill = Method), 
               alpha = 1, outlier.shape = NA) +  # 绘制箱线图
  geom_jitter(
    aes(group = interaction(Metric, Method), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.2,   # 控制抖动幅度
      dodge.width = 0.8     # 控制不同 Scenario 的间距
    ),size = 3,alpha = 1,color = "black",#stroke = 0.5
  )+
  scale_fill_manual(values = c('raw_atac' = 'grey', "Seurat" = "skyblue", 'LS_Lab' = "pink",
                               'MultiVI' = '#CC99E6', "XChrom" = "#FFB07C")) +  # 自定义箱线图颜色
  scale_shape_manual(
    values = c(
      "h_pbmc" = 1,    # 空心圆
      "m_brain" = 2,    # 空心三角
      "h_brain" = 5,    # 空心菱形
      "h_gonads" = 0,    # 空心方框
      "m_palates" = 6    # 空心倒三角
    ))+
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  # 增加 x 轴 Metric 之间的间距
  guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1)) +  # 图例顺序
  facet_wrap(~MetricType, scales = "free")+
  theme_set(theme_bw())+
  # theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold", color = "black", margin = margin(t =2)),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold",margin = margin(b = 10)),
    legend.text = element_text(size = 18),
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3),color = "black", size = 0.5),  # 设置分面标签背景
    strip.text = element_text(size = 18, face = "bold", color = "black"),  # 设置分面标签的文字样式
    strip.text.x = element_text(size = 18, face = "bold",color = "black"),  # 设置分面标签的文本样式
  ) +
  labs(title = "Cross-cell Prediction", x = "", y = "Value") +
  geom_text(
    data = annotation_df[annotation_df$p_value<0.05,],
    aes(x = Metric, y = y_position,group = MetricType,  # x=1.5 将注释居中于分面面板
        label = paste0(p_value_label, "\n95%CI = [",
                       round(lower_ci, 2), ",",
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 5,
    color = "black",hjust = 0.5)



ggsave("intra_method_nsls_box.pdf", plot = last_plot(), 
       width = 20, height = 8, units = "in")
