library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)

df <- read_csv("test_r4_nsls.csv")

df$Model <- factor(df$Model, levels = c('raw_atac',"human_model", "mouse_model"))
df$Species <- factor(df$Species, levels = c('human','mouse','macaque','marmoset'))

metric_ = 'ns(100)'
r = 4

##############----------ns/ls--------#########
## "ns(10)","ns(50)" / "ns(100)"
## "ls(10)","ls(50)" / "ls(100)"

#######################
filtered_df <- df %>% 
  filter(Metric == metric_) %>% 
  group_by(Species, Samples) %>%
  filter(n() != 2) 

y_pos_all <- filtered_df %>%
  group_by(Species) %>%
  summarise(max_val = max(Value)) %>%
  mutate(y_position = max_val + 0.1)



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




######################
# 1: human_model vs mouse_model
annotation_df1 <- filtered_df %>%
  filter(Model %in% c("human_model", "mouse_model")) %>%
  pivot_wider(names_from = Model, values_from = Value) %>%
  mutate(diff = mouse_model - human_model) %>%
  group_by(Species) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE)/sqrt(n()),
    lower_ci = mean_diff - qt(0.975, n()-1)*se_diff,
    upper_ci = mean_diff + qt(0.975, n()-1)*se_diff,
    test = list(perform_tests(diff)),
    .groups = "drop"
  ) %>%
  mutate(comparison = "human vs mouse", x_pos = 2.5,
         p_value = as.numeric(purrr::map_dbl(test, "p_value")),
         test_method = purrr::map_chr(test, "test_method"),
         p_value_label = sprintf("p = %.3f", p_value)
         ) %>% 
  left_join(y_pos_all, by = "Species")

# 2: mouse_model vs raw_atac
annotation_df2 <- filtered_df %>%
  filter(Model %in% c("mouse_model", "raw_atac")) %>%
  pivot_wider(names_from = Model, values_from = Value) %>%
  mutate(diff = mouse_model - raw_atac) %>%
  group_by(Species) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE)/sqrt(n()),
    lower_ci = mean_diff - qt(0.975, n()-1)*se_diff,
    upper_ci = mean_diff + qt(0.975, n()-1)*se_diff,
    test = list(perform_tests(diff)),
    .groups = "drop"
  ) %>%
  mutate(comparison = "mouse vs raw", x_pos = 2.0,
         p_value = as.numeric(purrr::map_dbl(test, "p_value")),
         test_method = purrr::map_chr(test, "test_method"),
         p_value_label = sprintf("p = %.3f", p_value)
         ) %>% 
  left_join(y_pos_all, by = "Species")

# 3: human_model vs raw_atac
annotation_df3 <- filtered_df %>%
  filter(Model %in% c("human_model", "raw_atac")) %>%
  pivot_wider(names_from = Model, values_from = Value) %>%
  mutate(diff = human_model - raw_atac) %>%
  group_by(Species) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE)/sqrt(n()),
    lower_ci = mean_diff - qt(0.975, n()-1)*se_diff,
    upper_ci = mean_diff + qt(0.975, n()-1)*se_diff,
    test = list(perform_tests(diff)),
    .groups = "drop"
  ) %>%
  mutate(comparison = "human vs raw", x_pos = 1.5,
         p_value = as.numeric(purrr::map_dbl(test, "p_value")),
         test_method = purrr::map_chr(test, "test_method"),
         p_value_label = sprintf("p = %.3f", p_value)
         ) %>% 
  left_join(y_pos_all, by = "Species")

#######################
combined_annot <- bind_rows(annotation_df1, annotation_df2, annotation_df3) %>%
  # filter(p_value < 0.05) %>%
  group_by(Species) %>%
  ungroup()

out <- combined_annot %>%
  mutate(test = map_chr(test, ~ toJSON(.x, auto_unbox = TRUE)))
write_csv(out, paste0(gsub("[()]", "", metric_),'_r',r,"_test.csv"))

##############----------PLOT--------#########
ggplot(filtered_df, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(width = 0.6, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  facet_wrap(~Species, scales = "free_x",strip.position = "bottom",nrow = 1) +
  scale_fill_manual(values = c('raw_atac'='grey', "human_model" = "orange", "mouse_model" = "#7B4F94")) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold", color = "black"), 
    strip.background = element_blank(),  
    panel.spacing = unit(0.5, "lines") 
  ) +
  labs(x = NULL, y = metric_) +
  geom_text(
    data = combined_annot,
    aes(x = x_pos, y = y_position, 
        label = paste0(p_value, "\n95%CI = [", 
                       round(lower_ci, 2), ",", 
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4, color = "black", vjust = 0.5
  )

ggsave(paste0(metric_,'.pdf'), plot = last_plot(), 
       width = 9, height = 4, units = "in",dpi=300)
