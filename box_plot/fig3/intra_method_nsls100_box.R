library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)


df <- read_csv('test_r4_intra_method_nsls100.csv')
r=4


# Convert Metric & Method & Dataset to factor and specify the order
df$Metric <- factor(df$Metric, levels = c("ns(10)", "ns(50)", "ns(100)", 
                                          "ls(10)", "ls(50)", "ls(100)"))
df$Method <- factor(df$Method, levels = c('raw_atac',"Seurat", "LS_Lab", "MultiVI", "XChrom"))
df$Dataset<-factor(df$Dataset, levels = c("h_pbmc","m_brain","h_brain",'h_gonads','m_palates'))

######################################
# Add 2 new column to group NS(k=10,50,100) and LS(k=10,50,100),respectively
df <- df %>%
  mutate(MetricType = ifelse(grepl("ns", Metric), "neighbor score", "label score"))%>%
  mutate(MetricType = factor(MetricType, levels = c("neighbor score", "label score")))

# Count each facet type
metric_counts <- df %>%
  group_by(MetricType) %>%
  summarise(count = n(), .groups = "drop")

#################################### facet ######################################
y_pos <- df %>% 
  group_by(Metric) %>% 
  summarise(max_val = max(Value))

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
  mutate(y_position = ifelse(p_value < 0.05, max_val + 0.1, NA_real_))



out <- annotation_df %>%
  mutate(test = map_chr(test, ~ toJSON(.x, auto_unbox = TRUE)))
write_csv(out, paste0('cross-cell_method_nsls','_r',r,"_test.csv"))




ggplot(df, aes(x = Metric, y = Value,)) +
  geom_boxplot(aes(group = interaction(Metric, Method), fill = Method), 
               alpha = 1, outlier.shape = NA) + 
  geom_jitter(
    aes(group = interaction(Metric, Method), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.2,  
      dodge.width = 0.8    
    ),size = 3,alpha = 1,color = "black",#stroke = 0.5
  )+
  scale_fill_manual(values = c('raw_atac' = 'grey', "Seurat" = "skyblue", 'LS_Lab' = "pink",
                               'MultiVI' = '#CC99E6', "XChrom" = "#FFB07C")) + 
  scale_shape_manual(
    values = c(
      "h_pbmc" = 1,   
      "m_brain" = 2,    
      "h_brain" = 5,    
      "h_gonads" = 0,   
      "m_palates" = 6   
    ))+
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  
  guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1)) + 
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
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3),color = "black", size = 0.5), 
    strip.text = element_text(size = 18, face = "bold", color = "black"),  
    strip.text.x = element_text(size = 18, face = "bold",color = "black"), 
  ) +
  labs(title = "Cross-cell Prediction", x = "", y = "Value") +
  geom_text(
    data = annotation_df[annotation_df$p_value<0.05,],
    aes(x = Metric, y = y_position,group = MetricType,  
        label = paste0(p_value_label, "\n95%CI = [",
                       round(lower_ci, 2), ",",
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 5,
    color = "black",hjust = 0.5)



ggsave(paste0('cross-cell_method_nsls_r',r,'.pdf'), plot = last_plot(), 
       width = 20, height = 8, units = "in")

