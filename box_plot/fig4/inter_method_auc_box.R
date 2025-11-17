library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)


df <- read_csv('test_r4_inter_method_auc.csv')
r=4

# Convert Metric & Method & Dataset to factor and specify the order
df$Metric <- factor(df$Metric, levels = c("auROC", "per cell auROC", "per peak auROC", 
                                          "auPRC", "per cell auPRC", "per peak auPRC"))
df$Method <- factor(df$Method, levels = c("Seurat", "LS_Lab", "MultiVI", "XChrom"))
df$Dataset<-factor(df$Dataset,levels = c('s1d1-s1d2','s1d1-s2d1','s2d5-s3d10','s3d3-s4d1','s4d8-s2d4'))
df %>%
  group_by(Metric, Method) %>%
  summarise(Dataset_Count = n_distinct(Dataset), .groups = "drop")

#################################### facet ######################################
# Add 2 new column to group auROC(overall,per-cell,per-peak) and auPRC(overall,per-cell,per-peak),respectively

df <- df %>%
  mutate(MetricType = ifelse(grepl("auROC", Metric), "auROC", "auPRC"))%>%
  mutate(MetricType = factor(MetricType, levels = c("auROC", "auPRC")))

# Obtain all the methods that need to be compared (excluding XChrom itself)
other_methods <- unique(df$Method[df$Method != "XChrom"])

# Extract the data of XChrom as the benchmark
xchrom_data <- df %>% 
  filter(Method == "XChrom") %>% 
  select(Dataset, Metric, MetricType, Value_XChrom = Value)

# Generate comparative data
comparison_data <- df %>%
  filter(Method != "XChrom") %>% 
  inner_join(xchrom_data, by = c("Dataset", "Metric", "MetricType")) %>%
  mutate(diff = Value_XChrom - Value)



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




# Calculate statistical results
annotation_df <- comparison_data %>%
  group_by(Metric, Method) %>%  # Group by Metric & Method
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
    p_value_label = ifelse(p_value < 0.05,
                           paste0("p = ", signif(p_value, digits = 3)),
                           NA_character_)
  ) %>%
  # Merge the MetricType and Y-axis position information
  left_join(distinct(df, Metric, MetricType), by = "Metric") %>% 
  left_join(
    df %>% group_by(Metric) %>% summarise(max_val = max(Value)),
    by = "Metric"
  ) %>% 
  mutate(y_position = ifelse(p_value < 0.05, max_val + 0.06, NA_real_))


out <- annotation_df %>%
  mutate(test = map_chr(test, ~ toJSON(.x, auto_unbox = TRUE)))
write_csv(out, paste0('interMethod_aucprc_r',r,"_test.csv"))




ggplot(df, aes(x = Metric, y = Value,)) +
  geom_boxplot(aes(group = interaction(Metric, Method), fill = Method), 
               alpha = 1, outlier.shape = NA) + 
  geom_jitter(
    aes(group = interaction(Metric, Method), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.2,   
      dodge.width = 0.8    
    ),size = 3,alpha = 0.8,color = "black",stroke = 0.8)+
  scale_fill_manual(values = c("Seurat" = "skyblue", 'LS_Lab' = "pink",'MultiVI'='#CC99E6',"XChrom" = "#FFB07C")) + 
  scale_shape_manual(
    values = c(
      "s1d1-s1d2" = 1,   
      "s1d1-s2d1" = 2,   
      "s4d8-s2d4" = 5,   
      "s3d3-s4d1" = 0,    
      "s2d5-s3d10" = 6    
    ))+
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) + 
  guides(shape = guide_legend(order = 1), fill = guide_legend(order = 2)) + 
  facet_wrap(~MetricType, scales = "free")+
  theme_set(theme_bw())+
  # theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black", margin = margin(t =2)),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold",margin = margin(b = 10)),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3),color = "black", size = 0.5),  
    strip.text = element_text(size = 12, face = "bold", color = "black"),  
    strip.text.x = element_text(size = 12, face = "bold",color = "black"), 
  ) +
  labs(title = "Metrics Across Scenario", x = "", y = "Value") +
  geom_text(
    data = annotation_df[annotation_df$p_value<0.05,],
    aes(x = Metric, y = y_position,group = MetricType, 
        label = paste0(p_value_label, "\n95%CI = [",
                       round(lower_ci, 2), ",",
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4,
    color = "black",hjust = 0.5)

ggsave(paste0('inter_method_aucprc_r',r,'.pdf'), plot = last_plot(), 
       width = 20, height = 6, units = "in",dpi=300)

