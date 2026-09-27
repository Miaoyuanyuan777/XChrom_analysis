library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

df <- read_csv('denoise_atac_nsls-ari-asw.csv')
if("Methods" %in% names(df)) {
  df <- df %>% rename(Method = Methods)
}

r <- 4

df$Metric <- factor(df$Metric, levels = c("ns(10)", "ns(50)", "ns(100)", 
                                          "ls(10)", "ls(50)", "ls(100)",
                                          "ari", "ami", "asw"))
df$Method <- factor(df$Method, levels = c("raw_ATAC", "SCALE", "scCASE", "scBasset", "XChrom"))
df$Dataset <- factor(df$Dataset, levels = c("h_PBMC", "m_brain", "h_brain", "h_gonads", "m_palates"))

df <- df %>%
  mutate(MetricType = case_when(
    grepl("ns", Metric) ~ "neighbor score",
    grepl("ls", Metric) ~ "label score",
    Metric %in% c("ari", "ami", "asw") ~ "clustering metrics",
    TRUE ~ "other"
  )) %>%
  mutate(MetricType = factor(MetricType, levels = c("neighbor score", "label score", "clustering metrics")))



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
  select(Dataset, Metric, MetricType, Method, Value) %>%
  pivot_wider(names_from = Method, values_from = Value)

methods_list <- levels(df$Method)

target_method <- "XChrom"
other_methods <- setdiff(methods_list, target_method)
# (d <- df_m[[m1]] - df_m[[m2]])
pairs <- lapply(other_methods, function(m) c(target_method, m)) 

results_list <- list()

for(m in unique(df$Metric)) {
  m_type <- unique(df$MetricType[df$Metric == m])
  df_m <- df_wide %>% filter(Metric == m)
  
  for(p in pairs) {
    m1 <- p[1] #  "XChrom"
    m2 <- p[2] # "raw_ATAC", "SCALE", "scCASE", "scBasset" 
    
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

write_csv(pairwise_df, paste0('denoise_','XChrom_vs_otherM_wilcox_test.csv'))


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


p <- ggplot(df, aes(x = Metric, y = Value)) +
  geom_boxplot(aes(fill = Method), alpha = 0.7, outlier.shape = NA, position = position_dodge(width = dodge_width)) + 
  geom_jitter(
    aes(shape = Dataset, group = Method),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = dodge_width),
    size = 1.5, alpha = 1, color = "black"
  ) +
  scale_fill_manual(values = c(
    'raw_ATAC' = 'grey',
    'scCASE'   = '#88D498', # "#A8E6CF"
    'SCALE'    = "#48D1CC",# '#984EA1',"#48D1CC","#7FFFD4","#45B39D"
    'scBasset' = '#1A66B3',
    'XChrom'   = '#FFB07C'
  )) + 
  scale_shape_manual(
    values = c(
      "h_PBMC" = 1,    
      "m_brain" = 2,    
      "h_brain" = 5,    
      "h_gonads" = 0,    
      "m_palates" = 6    
    )) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  
  guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1)) +  
  facet_wrap(~MetricType, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold", color = "black", margin = margin(t = 2)),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 10)),
    legend.text = element_text(size = 18),
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3), color = "black", linewidth = 0.5), 
    strip.text = element_text(size = 18, face = "bold", color = "black")
  ) +
  labs(title = "Cross-region Prediction", x = "", y = "Value")

if(nrow(plot_sig_df) > 0) {
  p <- p + 
    geom_segment(data = plot_sig_df, aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos), inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = plot_sig_df, aes(x = x_min, xend = x_min, y = y_pos, yend = y_tip), inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = plot_sig_df, aes(x = x_max, xend = x_max, y = y_pos, yend = y_tip), inherit.aes = FALSE, linewidth = 0.5) +
    geom_text(data = plot_sig_df, aes(x = (x_min + x_max)/2, y = y_pos - 0.01, label = signif), inherit.aes = FALSE, size = 6, vjust = 0, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
}

p

###################################################
########## 
##################################################

methods_group1 <- c("raw_ATAC", "scBasset", "XChrom")
methods_group2 <- c("raw_ATAC","scCASE", "SCALE","XChrom")

color_palette <- c(
  'raw_ATAC' = 'grey',
  'scCASE'   = '#88D498', 
  'SCALE'    = "#48D1CC",
  'scBasset' = '#1A66B3',
  'XChrom'   = '#FFB07C'
)


shape_palette <- c(
  "h_PBMC" = 1,    
  "m_brain" = 2,    
  "h_brain" = 5,    
  "h_gonads" = 0,    
  "m_palates" = 6    
)

plot_subgroup <- function(sub_methods, title_text, pairwise_df, df, metric_indices, metric_max, all_methods) {
  df_sub <- df %>% filter(Method %in% sub_methods)
  sig_df_sub <- pairwise_df %>% 
    filter(signif != "ns" & Group1 %in% sub_methods & Group2 %in% sub_methods)
  
  k_total <- length(all_methods)
  dodge_width <- 0.8
  offsets_all <- seq(-dodge_width/2 + dodge_width/(2*k_total), dodge_width/2 - dodge_width/(2*k_total), length.out = k_total)
  names(offsets_all) <- all_methods
  
  plot_sig_sub <- sig_df_sub %>%
    left_join(metric_indices, by = c("MetricType", "Metric")) %>%
    left_join(metric_max, by = "Metric") %>%
    mutate(
      x1_offset = offsets_all[Group1],
      x2_offset = offsets_all[Group2],
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
  
  p <- ggplot(df_sub, aes(x = Metric, y = Value)) +
    geom_boxplot(aes(fill = Method), alpha = 1, outlier.shape = NA, position = position_dodge(width = dodge_width)) + 
    geom_jitter(
      aes(shape = Dataset, group = Method),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = dodge_width),
      size = 3, alpha = 1, color = "black"
    ) +
    scale_fill_manual(values = color_palette, limits = all_methods) +
    scale_shape_manual(values = shape_palette) +
    scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  
    guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1)) +  
    facet_wrap(~MetricType, scales = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold", color = "black", margin = margin(t = 2)),
      axis.text.y = element_text(size = 18, face = "bold", color = "black"),
      axis.title = element_text(size = 18, color = "black", face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 10)),
      legend.text = element_text(size = 18),
      strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3), color = "black", linewidth = 0.5), 
      strip.text = element_text(size = 18, face = "bold", color = "black")
    ) +
    labs(title = title_text, x = "", y = "Value")
  
  if(nrow(plot_sig_sub) > 0) {
    p <- p + 
      geom_segment(data = plot_sig_sub, aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos), inherit.aes = FALSE, linewidth = 0.5) +
      geom_segment(data = plot_sig_sub, aes(x = x_min, xend = x_min, y = y_pos, yend = y_tip), inherit.aes = FALSE, linewidth = 0.5) +
      geom_segment(data = plot_sig_sub, aes(x = x_max, xend = x_max, y = y_pos, yend = y_tip), inherit.aes = FALSE, linewidth = 0.5) +
      geom_text(data = plot_sig_sub, aes(x = (x_min + x_max)/2, y = y_pos - 0.01, label = signif), inherit.aes = FALSE, size = 6, vjust = 0, fontface = "bold") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
  }
  
  return(p)
}

p1 <- plot_subgroup(
  sub_methods = methods_group1, 
  title_text = "Cross-region Prediction (Group 1)", 
  pairwise_df = pairwise_df, 
  df = df, 
  metric_indices = metric_indices, 
  metric_max = metric_max,
  all_methods = methods_list 
)
print(p1)
ggsave(paste0('crosscell_', r, '_XChrom_vs_scBasset.pdf'), p1, width = 10, height = 6)

p2 <- plot_subgroup(
  sub_methods = methods_group2, 
  title_text = "Cross-region Prediction (Group 2)", 
  pairwise_df = pairwise_df, 
  df = df, 
  metric_indices = metric_indices, 
  metric_max = metric_max,
  all_methods = methods_list
)
print(p2)
ggsave(paste0('crosscell_', r, '_XChrom_vs_scCASE.pdf'), p2, width = 12, height = 6)