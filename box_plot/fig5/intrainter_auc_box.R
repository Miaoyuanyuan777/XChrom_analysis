library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)


df <- read_csv("test_r4_intrainter_auc.csv")
r=4
# Convert Metric & Scenario & Dataset to factor and specify the order
df$Metric <- factor(df$Metric, levels = c("auROC", "per cell auROC", "per peak auROC", 
                                          "auPRC", "per cell auPRC", "per peak auPRC"))
df$Scenario<-factor(df$Scenario,levels=c('inter','intra'))
df$Dataset<-factor(df$Dataset,levels = c('s1d2','s2d1','s3d10','s4d1','s2d4'))

#################################### facet ######################################
# Add 2 new column to group auROC(overall,per-cell,per-peak) and auPRC(overall,per-cell,per-peak),respectively

df <- df %>%
  mutate(MetricType = ifelse(grepl("auROC", Metric), "auROC", "auPRC"))%>%
  mutate(MetricType = factor(MetricType, levels = c("auROC", "auPRC")))

#######################################################
y_pos <- df %>% 
  group_by(Metric) %>% 
  summarise(max_val = max(Value))

perform_tests <- function(differences, digits_round = NULL, exact = FALSE, correct = TRUE) {
  d <- differences[is.finite(differences)]
  if (!is.null(digits_round)) d <- round(d, digits_round)
  
  n <- length(d)
  mean_d <- mean(d)
  sd_d <- sd(d)
  
  res <- list(
    p_value = NA_real_, test_method = "insufficient_n",
    estimate_diff = NA_real_, estimate_type = NA_character_, #  (Mean diff / Pseudo-median)
    effect_size = NA_real_, effect_size_type = NA_character_, #  (Cohen's d / Wilcoxon r)
    lower_ci = NA_real_, upper_ci = NA_real_
  )
  
  if (n < 2) return(res)
  
  if (sd_d == 0) {
    if (all(abs(d) < .Machine$double.eps)) {
      res$p_value <- 1; res$test_method <- "no-diff"
      res$estimate_diff <- 0; res$estimate_type <- "mean_diff"
      res$effect_size <- 0; res$effect_size_type <- "Cohen_d"
      res$lower_ci <- 0; res$upper_ci <- 0
      return(res)
    } else {
      k <- sum(d > 0); n_nonzero <- sum(d != 0)
      bt <- binom.test(k, n_nonzero, p = 0.5, alternative = "two.sided")
      res$p_value <- bt$p.value; res$test_method <- "Sign test (binom)"
      res$estimate_diff <- mean_d; res$estimate_type <- "mean_diff"
      res$effect_size <- k / n_nonzero; res$effect_size_type <- "Proportion_positive"
      res$lower_ci <- mean_d - qt(0.975, df = n-1) * (sd_d / sqrt(n))
      res$upper_ci <- mean_d + qt(0.975, df = n-1) * (sd_d / sqrt(n))
      return(res)
    }
  }
  
  shapiro_p <- if (n >= 3) {
    tryCatch(shapiro.test(d)$p.value, error = function(e) NA_real_)
  } else NA_real_
  
  if (!is.na(shapiro_p) && shapiro_p > 0.05) {
    # t-test
    tt <- t.test(d, mu = 0)
    res$p_value <- tt$p.value
    res$test_method <- "paired t-test"
    
    # Mean Difference
    res$estimate_diff <- tt$estimate[["mean of x"]] 
    res$estimate_type <- "Mean_Difference"
    
    # Cohen's d
    res$effect_size <- mean_d / sd_d  
    res$effect_size_type <- "Cohen_d"
    
    res$lower_ci <- tt$conf.int[1]
    res$upper_ci <- tt$conf.int[2]
  } else {
    # Wilcoxon
    wt <- suppressWarnings(tryCatch(
      wilcox.test(d, alternative = "two.sided", mu = 0, exact = FALSE, correct = TRUE, conf.int = TRUE),
      error = function(e) wilcox.test(d, alternative = "two.sided", mu = 0, exact = FALSE, correct = TRUE)
    ))
    
    res$p_value <- wt$p.value
    res$test_method <- "Wilcoxon signed-rank"
    
    # Pseudo-median
    if ("estimate" %in% names(wt)) {
      res$estimate_diff <- wt$estimate[["(pseudo)median"]]
      res$estimate_type <- "Pseudo_median"
    } else {
      res$estimate_diff <- median(d)
      res$estimate_type <- "Median_Difference" # Fallback
    }
    
    # Wilcoxon r = Z / sqrt(N)
    z_val <- qnorm(wt$p.value / 2)
    res$effect_size <- abs(z_val) / sqrt(n)
    res$effect_size_type <- "Wilcoxon_r"
    
    if ("conf.int" %in% names(wt)) {
      res$lower_ci <- wt$conf.int[1]
      res$upper_ci <- wt$conf.int[2]
    } else {
      res$lower_ci <- mean_d - qt(0.975, df = n-1) * (sd_d / sqrt(n))
      res$upper_ci <- mean_d + qt(0.975, df = n-1) * (sd_d / sqrt(n))
    }
  }
  return(res)
}


df_wide <- df %>% 
  select(Dataset, Metric, MetricType, Scenario, Value) %>%
  pivot_wider(names_from = Scenario, values_from = Value)

methods_list <- levels(df$Scenario)

target_method <- "inter"
other_methods <- setdiff(methods_list, target_method)
# (d <- df_m[[m1]] - df_m[[m2]])  Scenario inter - intra
pairs <- lapply(other_methods, function(m) c(target_method, m)) 

results_list <- list()

for(m in unique(df$Metric)) {
  m_type <- unique(df$MetricType[df$Metric == m])
  df_m <- df_wide %>% filter(Metric == m)
  
  for(p in pairs) {
    m1 <- p[1] 
    m2 <- p[2]
    
    if(m1 %in% colnames(df_m) && m2 %in% colnames(df_m)) {
      d <- df_m[[m1]] - df_m[[m2]]
      
      test_res <- perform_tests(d)
      n_val <- sum(!is.na(d))
      
      results_list[[length(results_list) + 1]] <- tibble(
        Metric = m,
        MetricType = m_type,
        Group1 = m1,
        Group2 = m2,
        N = n_val, 
        estimate_diff = test_res$estimate_diff, 
        estimate_type = test_res$estimate_type,
        lower_ci = test_res$lower_ci,  
        upper_ci = test_res$upper_ci,
        effect_size = test_res$effect_size,
        effect_size_type = test_res$effect_size_type,
        test_method = test_res$test_method,
        p_value = test_res$p_value
      )
    }
  }
}

pairwise_df <- bind_rows(results_list) %>%
  mutate(
    signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )
write_csv(pairwise_df, paste0('interintra_auc','_t_wilcox_test.csv'))
sig_df <- pairwise_df %>% filter(signif != "ns")
metric_indices <- df %>% 
  distinct(MetricType, Metric) %>% 
  group_by(MetricType) %>% 
  arrange(Metric) %>% 
  mutate(metric_x = row_number()) %>% 
  ungroup()

k <- length(methods_list)
dodge_width <- 0.8
offsets <- seq(-dodge_width/2 + dodge_width/(2*k), dodge_width/2 - dodge_width/(2*k), length.out = k)
names(offsets) <- methods_list
metric_max <- df %>% group_by(Metric) %>% summarise(max_val = max(Value, na.rm = TRUE))
plot_sig_df <- sig_df %>%
  left_join(metric_indices, by = c("MetricType", "Metric")) %>%
  left_join(metric_max, by = "Metric") %>%
  mutate(
    x1_offset = offsets[Group1],
    x2_offset = offsets[Group2],
    x_min = metric_x + pmin(x1_offset, x2_offset),
    x_max = metric_x + pmax(x1_offset, x2_offset),
    dist = abs(x1_offset - x2_offset) 
  ) %>%

  arrange(MetricType, Metric, dist) %>%
  group_by(MetricType, Metric) %>%
  mutate(
    step = row_number(),
    y_pos = max_val + 0.05 + 0.04 * (step - 1), 
    y_tip = y_pos - 0.01 
  ) %>%
  ungroup()



p <- ggplot(df, aes(x = Metric, y = Value,)) +
  geom_boxplot(aes(group = interaction(Metric, Scenario), fill = Scenario), 
               alpha = 1, outlier.shape = NA) + 
  geom_jitter(
    aes(group = interaction(Metric, Scenario), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.3,  
      dodge.width = 0.6     
    ),size = 3,alpha = 0.8,color = "black",stroke = 1)+
  scale_fill_manual(values = c("intra" = "#aed09c", "inter" = "#FFD700")) + 
  scale_shape_manual(
    values = c(
      "s1d2" = 1,   
      "s2d1" = 2,   
      "s2d4" = 5,   
      "s4d1" = 0,   
      "s3d10" = 6   
    ))+
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  
  guides(shape = guide_legend(order = 1), fill = guide_legend(order = 2)) + 
  facet_wrap(~MetricType, scales = "free")+
  theme_set(theme_bw())+
  # theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold", color = "black", margin = margin(t =2)),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold",margin = margin(b = 10)),
    legend.text = element_text(size = 18),
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3),color = "black", size = 0.5),  
    strip.text = element_text(size = 18, face = "bold", color = "black"),  
    strip.text.x = element_text(size = 18, face = "bold",color = "black"), 
  ) +
  labs(title = "Metrics Across Scenario", x = "", y = "Value")



if(nrow(plot_sig_df) > 0) {
  p <- p + 
    geom_segment(data = plot_sig_df, aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos), inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = plot_sig_df, aes(x = x_min, xend = x_min, y = y_pos, yend = y_tip), inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = plot_sig_df, aes(x = x_max, xend = x_max, y = y_pos, yend = y_tip), inherit.aes = FALSE, linewidth = 0.5) +
    geom_text(data = plot_sig_df, aes(x = (x_min + x_max)/2, y = y_pos - 0.00, label = signif), inherit.aes = FALSE, size = 6, vjust = 0, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
}

p

ggsave(paste0("intrainter_auc_r",r,'.pdf'), plot = p, 
       width = 10, height = 6, units = "in",dpi = 300)

