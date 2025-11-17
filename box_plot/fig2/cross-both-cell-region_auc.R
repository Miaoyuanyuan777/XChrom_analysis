library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)



df <- read_csv('test_r4_cross-both-cell-region_auc.csv')
r=4

# Convert Metric & Method & Test to factor and specify the order
df$Metric <- factor(df$Metric, levels = c("auROC", "per cell auROC", "per peak auROC", 
                                          "auPRC", "per cell auPRC", "per peak auPRC"))
df$Dataset<-factor(df$Dataset, levels = c("h_pbmc","m_brain","h_brain",'h_gonads','m_palates'))
df$Test<-factor(df$Test,levels=c('cross-cell','cross-region','cross-both'))

######################################
# Add 2 new column to group auROC(overall,per-cell,per-peak) and auPRC(overall,per-cell,per-peak),respectively
df <- df %>%
  mutate(MetricType = ifelse(grepl("auPRC", Metric), "auPRC", "auROC"))%>%
  mutate(MetricType = factor(MetricType, levels = c("auROC", "auPRC")))

# Count the quantity of each facet type
metric_counts <- df %>%
  group_by(MetricType) %>%
  summarise(count = n(), .groups = "drop")

######################################
y_pos <- df %>% 
  group_by(Metric) %>% 
  summarise(max_val = max(Value))

# Obtain all the methods that need to be compared (excluding the cross-both itself)
other_methods <- unique(df$Test[df$Test != "cross-both"])

# Extract the data of cross-both as the benchmark
crossboth_data <- df %>% 
  filter(Test == "cross-both") %>% 
  select(Dataset, Metric, MetricType, Value_cross_both = Value)

# Generate comparative data
comparison_data <- df %>%
  filter(Test != "cross-both") %>% 
  inner_join(crossboth_data, by = c("Dataset", "Metric", "MetricType")) %>%
  mutate(diff = Value_cross_both - Value)




perform_tests <- function(differences, digits_round = r,
                          exact = FALSE, correct = TRUE) {
  d <- differences[is.finite(differences)]
  if (!is.null(digits_round)) d <- round(d, digits_round)  
  
  n <- length(d)
  if (n < 2) return(list(p_value = NA_real_, test_method = "insufficient_n"))
  
  if (sd(d) == 0) {
    if (all(abs(d) < .Machine$double.eps)) {
      return(list(p_value = 1, test_method = "no-diff"))
    } else {
      k <- sum(d > 0); n_nonzero <- sum(d != 0)
      p <- binom.test(k, n_nonzero, p = 0.5, alternative = "two.sided")$p.value
      return(list(p_value = p, test_method = "Sign test (binom)"))
    }
  }
  shapiro_p <- if (n >= 3) {
    tryCatch(shapiro.test(d)$p.value, error = function(e) NA_real_)
  } else NA_real_
  
  if (!is.na(shapiro_p) && shapiro_p > 0.05) {
    p <- t.test(d, mu = 0)$p.value
    return(list(p_value = p, test_method = "paired t-test"))
  } else {
    p <- suppressWarnings(wilcox.test(d,alternative = "two.sided",mu = 0, exact = FALSE, correct = TRUE)$p.value)
    return(list(p_value = p, test_method = "Wilcoxon signed-rank"))
  }
}





# Calculate statistical results
annotation_df <- comparison_data %>%
  group_by(Metric, Test) %>%  # Group by Metric & Test
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
  mutate(y_position = ifelse(p_value < 0.05, max_val + 0.04, NA_real_))




out <- annotation_df %>%
  mutate(test = map_chr(test, ~ toJSON(.x, auto_unbox = TRUE)))
write_csv(out, paste0('cross-both-cell-region_aucprc_r',r,"_test.csv"))






ggplot(df, aes(x = Metric, y = Value)) +
  geom_boxplot(aes(group = interaction(Metric, Test), fill = Test),
               alpha = 1, outlier.shape = NA) +
  geom_jitter(
    aes(group = interaction(Metric, Test), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ), size = 3, alpha = 1, color = "black"
  ) +
  scale_fill_manual(values = c("cross-cell" = "#FF7F50", 'cross-both' = '#FFA500', 'cross-region' = '#FFB07C')) +
  scale_shape_manual(
    values = c(
      "h_pbmc" = 1,
      "m_brain" = 2,
      "h_brain" = 5,
      "h_gonads" = 0,
      "m_palates" = 6
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  facet_wrap(~MetricType, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16, face = "bold", color = "black", margin = margin(t = 2)),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
    legend.text = element_text(size = 16),
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3), color = "black", size = 0.5),
    strip.text = element_text(size = 16, face = "bold", color = "black"),
    strip.text.x = element_text(size = 16, face = "bold", color = "black"),
    panel.border = element_rect(size = 0.5, color = "black")
  ) +
  labs(title = "Cross-both/cell/region Prediction", x = "", y = "Value") + 
  geom_text(
    data = annotation_df %>% filter(p_value < 0.05),
    aes(x = Metric, y = y_position, label = paste0(p_value_label, "\n95%CI = [", 
                                                   round(lower_ci, 2), ",", 
                                                   round(upper_ci, 2), "]")),
    inherit.aes = FALSE, size = 4,
    color = "black", hjust = 0.5
  )



ggsave(paste0('cross-both_auc_r',r,'.pdf'), plot = last_plot(), 
       width = 10, height = 6, units = "in",dpi=300)

