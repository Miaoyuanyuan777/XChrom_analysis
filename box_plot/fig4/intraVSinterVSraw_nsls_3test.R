library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)


df <- read_csv("test_r4_intrainter_nsls.csv")
r=4

# 将 Metric 转换为因子并指定顺序
df$Metric <- factor(df$Metric, levels = c("ns(10)", "ns(50)", "ns(100)", 
                                          "ls(10)", "ls(50)", "ls(100)"))
df$Scenario<-factor(df$Scenario,levels=c('raw_atac','inter','intra'))
df$Dataset<-factor(df$Dataset,levels = c('s1d2','s2d1','s3d10','s4d1','s2d4'))


# 绘制分面图，使用不同的纵坐标尺度并去除没有数据的分面轴
df <- df %>%
  mutate(MetricType = ifelse(grepl("ns", Metric), "neighbor score", "label score"))%>%
  mutate(MetricType = factor(MetricType, levels = c("neighbor score", "label score")))

y_pos <- df %>% 
  group_by(Metric) %>% 
  summarise(max_val = max(Value))


perform_tests <- function(differences, digits_round = r,
                          exact = FALSE, correct = TRUE) {
  d <- differences[is.finite(differences)]
  if (!is.null(digits_round)) d <- round(d, digits_round)  # 固化 ties（关键！）
  
  n <- length(d)
  if (n < 2) return(list(p_value = NA_real_, test_method = "insufficient_n"))
  
  if (sd(d) == 0) {
    if (all(abs(d) < .Machine$double.eps)) {
      return(list(p_value = 1, test_method = "no-diff"))
    } else {
      # 常数非零，用精确符号检验（强烈建议报告这个 p，一眼看方向占比）
      k <- sum(d > 0); n_nonzero <- sum(d != 0)
      p <- binom.test(k, n_nonzero, p = 0.5, alternative = "two.sided")$p.value
      return(list(p_value = p, test_method = "Sign test (binom)"))
    }
  }
  # 非常数向量再做 Shapiro；n>=3 才有意义
  shapiro_p <- if (n >= 3) {
    tryCatch(shapiro.test(d)$p.value, error = function(e) NA_real_)
  } else NA_real_
  
  if (!is.na(shapiro_p) && shapiro_p > 0.05) {
    # 配对 t 等价于对差值做单样本 t：mu=0
    p <- t.test(d, mu = 0)$p.value
    return(list(p_value = p, test_method = "paired t-test"))
  } else {
    # 非正态或样本太少：Wilcoxon，考虑 ties
    p <- suppressWarnings(wilcox.test(d,alternative = "two.sided",mu = 0, exact = FALSE, correct = TRUE)$p.value)
    return(list(p_value = p, test_method = "Wilcoxon signed-rank"))
  }
}


## intra vs inter
annotation_df <- df %>%
  pivot_wider(names_from = Scenario, values_from = Value) %>%
  mutate(diff =  intra - inter) %>%
  group_by(Metric) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    test = list(perform_tests(diff)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = as.numeric(purrr::map_dbl(test, "p_value")),
    test_method = purrr::map_chr(test, "test_method"),
    p_value_label=ifelse(p_value < 0.05,
                              paste0("p = ", signif(p_value, digits = 3)
                              ),NA_character_))%>%
  mutate(
    comparison = "intra vs inter",
    MetricType = factor(case_when(
    grepl("ls", Metric) ~ "label score",
    grepl("ns", Metric) ~ "neighbor score"),
    levels = c("neighbor score", "label score")))

# 合并标注数据
annotation_df <- left_join(annotation_df, y_pos, by = "Metric")
annotation_df$y_position <- ifelse(annotation_df$p_value < 0.05,
                                   annotation_df$max_val + 0.15, 
                                   NA_real_) 



## check inter vs raw_ATAC
annotation_df2 <- df %>%
  pivot_wider(names_from = Scenario, values_from = Value) %>%
  mutate(diff =  inter - raw_atac) %>%
  group_by(Metric) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    test = list(perform_tests(diff)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = as.numeric(purrr::map_dbl(test, "p_value")),
    test_method = purrr::map_chr(test, "test_method"),
    p_value_label=ifelse(p_value < 0.05,
                              paste0("p = ", signif(p_value, digits = 3)
                              ),NA_character_))%>%
  mutate(
    comparison = "inter vs raw_ATAC",
    MetricType = factor(case_when(
    grepl("ls", Metric) ~ "label score",
    grepl("ns", Metric) ~ "neighbor score"),
    levels = c("neighbor score", "label score")))

# 合并标注数据
annotation_df2 <- left_join(annotation_df2, y_pos, by = "Metric")
annotation_df2$y_position <- ifelse(annotation_df2$p_value < 0.05,
                                   annotation_df2$max_val + 0.04, 
                                   NA_real_) 

## check intra vs raw_ATAC
annotation_df3 <- df %>%
  pivot_wider(names_from = Scenario, values_from = Value) %>%
  mutate(diff =  intra - raw_atac) %>%
  group_by(Metric) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    test = list(perform_tests(diff)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = as.numeric(purrr::map_dbl(test, "p_value")),
    test_method = purrr::map_chr(test, "test_method"),
    p_value_label=ifelse(p_value < 0.05,
                              paste0("p = ", signif(p_value, digits = 3)
                              ),NA_character_))%>%
  mutate(
    comparison = "intra vs raw_ATAC",
    MetricType = factor(case_when(
    grepl("ls", Metric) ~ "label score",
    grepl("ns", Metric) ~ "neighbor score"),
    levels = c("neighbor score", "label score")))

# 合并标注数据
annotation_df3 <- left_join(annotation_df3, y_pos, by = "Metric")
annotation_df3$y_position <- ifelse(annotation_df3$p_value < 0.05,
                                    annotation_df3$max_val + 0.04, 
                                    NA_real_) 

combined_annot <- bind_rows(annotation_df, annotation_df2, annotation_df3) %>%
  #filter(p_value < 0.05) %>%
  group_by(Metric) %>%
  ungroup()

out <- combined_annot %>%
  mutate(test = map_chr(test, ~ toJSON(.x, auto_unbox = TRUE)))
write_csv(out, paste0('intrainter_nsls_r','_r',r,"_test.csv"))


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
    data = combined_annot[combined_annot$p_value<0.05,],
    aes(x = Metric, y = y_position,group = MetricType,  # x=1.5 将注释居中于分面面板
        label = paste0(p_value_label, "\n95%CI = [", 
                       round(lower_ci, 2), ",", 
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4,
    color = "black",hjust = 0.5)

ggsave(paste0("intrainter_nsls_r",r,'.pdf'), plot = last_plot(), 
       width = 10, height = 6, units = "in",dpi=300)

