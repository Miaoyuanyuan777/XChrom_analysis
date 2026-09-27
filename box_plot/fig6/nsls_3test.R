library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)

df <- read_csv("test_r4_nsls.csv")

df$Model <- factor(df$Model, levels = c('raw_atac',"human_model", "mouse_model"))
df$Species <- factor(df$Species, levels = c('human','mouse','macaque','marmoset'))
r <- 4

perform_tests <- function(differences, digits_round = r, exact = FALSE, correct = TRUE) {
  d <- differences[is.finite(differences)]
  if (!is.null(digits_round)) d <- round(d, digits_round) 
  
  n <- length(d)
  mean_d <- mean(d)
  sd_d <- sd(d)
  
  res <- list(
    p_value = NA_real_, test_method = "insufficient_n",
    estimate_diff = NA_real_, estimate_type = NA_character_, 
    effect_size = NA_real_, effect_size_type = NA_character_, 
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
    tt <- t.test(d, mu = 0)
    res$p_value <- tt$p.value
    res$test_method <- "paired t-test"
    res$estimate_diff <- as.numeric(tt$estimate) 
    res$estimate_type <- "Mean_Difference"
    res$effect_size <- mean_d / sd_d  # Cohen's d
    res$effect_size_type <- "Cohen_d"
    res$lower_ci <- tt$conf.int[1]
    res$upper_ci <- tt$conf.int[2]
  } else {
    wt <- suppressWarnings(tryCatch(
      wilcox.test(d, alternative = "two.sided", mu = 0, exact = FALSE, correct = TRUE, conf.int = TRUE),
      error = function(e) wilcox.test(d, alternative = "two.sided", mu = 0, exact = FALSE, correct = TRUE)
    ))
    
    res$p_value <- wt$p.value
    res$test_method <- "Wilcoxon signed-rank"
    
    if ("estimate" %in% names(wt)) {
      res$estimate_diff <- as.numeric(wt$estimate) 
      res$estimate_type <- "Pseudo_median"
    } else {
      res$estimate_diff <- median(d)
      res$estimate_type <- "Median_Difference"
    }
    
    # Wilcoxon r = Z / sqrt(N)
    z_val <- qnorm(wt$p.value / 2)
    res$effect_size <- abs(z_val) / sqrt(n)
    res$effect_size_type <- "Wilcoxon_r"
    
    if ("conf.int" %in% names(wt)) {
      res$lower_ci <- wt$conf.int[1]
      res$upper_ci <- wt$conf.int[2]
      if (is.na(res$lower_ci) || is.nan(res$lower_ci)) {
        res$lower_ci <- mean_d - qt(0.975, df = n-1) * (sd_d / sqrt(n))
        res$upper_ci <- mean_d + qt(0.975, df = n-1) * (sd_d / sqrt(n))
      }
    } else {
      res$lower_ci <- mean_d - qt(0.975, df = n-1) * (sd_d / sqrt(n))
      res$upper_ci <- mean_d + qt(0.975, df = n-1) * (sd_d / sqrt(n))
    }
  }
  return(res)
}

compare_pairs <- list(
  c(Group2 = "raw_atac",    Group1 = "human_model"),
  c(Group2 = "raw_atac",    Group1 = "mouse_model"),
  c(Group2 = "mouse_model", Group1 = "human_model") 
)

test_data <- map_dfr(compare_pairs, function(pair) {
  g2 <- pair["Group2"]
  g1 <- pair["Group1"]
  df_g2 <- df %>% filter(Model == g2 & !is.na(Value)) %>% select(Metric, Species, Samples, Value_Group2 = Value)
  df_g1 <- df %>% filter(Model == g1 & !is.na(Value)) %>% select(Metric, Species, Samples, Value_Group1 = Value)
  inner_join(df_g2, df_g1, by = c("Metric", "Species", "Samples")) %>%
    mutate(
      Group2 = g2,
      Group1 = g1,
      diff = Value_Group2 - Value_Group1 
    )
})

# ==========================================
# ==========================================
annotation_all <- test_data %>%
  group_by(Metric, Species, Group1, Group2) %>%
  summarise(
    N        = n(),
    test_res = list(perform_tests(diff, digits_round = r)), 
    .groups  = "drop"
  ) %>%
  mutate(
    estimate_diff    = map_dbl(test_res, "estimate_diff"),
    estimate_type    = map_chr(test_res, "estimate_type"),
    lower_ci         = map_dbl(test_res, "lower_ci"),
    upper_ci         = map_dbl(test_res, "upper_ci"),
    effect_size      = map_dbl(test_res, "effect_size"),
    effect_size_type = map_chr(test_res, "effect_size_type"),
    test_method      = map_chr(test_res, "test_method"),
    p_value          = map_dbl(test_res, "p_value"),
    signif           = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(-test_res)

write_csv(annotation_all, paste0("all_metrics_test_results_nsls", ".csv"))

# ==========================================
# ==========================================
all_metrics <- unique(df$Metric)
methods_list <- c('raw_atac', 'human_model', 'mouse_model')

for (metric_ in all_metrics) {
  
  message(paste0("Generating plot for: ", metric_))
  
  df_sub <- df %>% filter(Metric == metric_, !is.na(Value))
  anno_sub <- annotation_all %>% filter(Metric == metric_)
  
  species_max <- df_sub %>% 
    group_by(Species) %>% 
    summarise(max_val = max(Value, na.rm = TRUE), .groups = "drop")
  
  anno_sig <- anno_sub %>% filter(signif != "ns")
  if (nrow(anno_sig) > 0) {
    anno_sig <- anno_sig %>%
      left_join(species_max, by = "Species") %>%
      mutate(
        x1 = as.numeric(factor(Group1, levels = methods_list)),
        x2 = as.numeric(factor(Group2, levels = methods_list)),
        x_min = pmin(x1, x2),
        x_max = pmax(x1, x2),
        dist = abs(x1 - x2)
      ) %>%
      arrange(Species, dist, x_min) %>% 
      group_by(Species) %>%
      mutate(
        step   = row_number(),
        y_line = max_val + 0.02 + 0.04 * (step - 1),  
        y_tip  = y_line - 0.015,                      
        y_star = y_line + 0.002                       
      ) %>%
      ungroup()
  }
  
  p <- ggplot(df_sub, aes(x = Model, y = Value, fill = Model)) +
    geom_boxplot(width = 0.6, alpha = 0.9, outlier.shape = NA) + 
    geom_jitter(width = 0.2, size = 2, alpha = 0.9) + 
    facet_wrap(~ Species, scales = "free_x", strip.position = "bottom", nrow = 1) + 
    scale_fill_manual(values = c("raw_atac"    = "grey",
                                 "human_model" = "orange", 
                                 "mouse_model" = "#7B4F94")) +
    theme_bw() +
    theme(
      axis.text.x      = element_blank(), 
      axis.ticks.x     = element_blank(), 
      axis.text.y      = element_text(size = 13, color = "black"),
      legend.text      = element_text(size = 16),
      axis.title       = element_text(size = 16, color = "black", face = "bold"),
      legend.title     = element_blank(),
      plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
      strip.text       = element_text(size = 16, face = "bold", color = "black"), 
      strip.background = element_blank(), 
      panel.spacing    = unit(0.5, "lines"),
      legend.position  = "right" 
    ) +
    labs(title = "", x = "", y = metric_)
  
  if (nrow(anno_sig) > 0) {
    p <- p +
      geom_segment(data = anno_sig, aes(x = x_min, xend = x_max, y = y_line, yend = y_line), 
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = anno_sig, aes(x = x_min, xend = x_min, y = y_line, yend = y_tip), 
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = anno_sig, aes(x = x_max, xend = x_max, y = y_line, yend = y_tip), 
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = anno_sig, aes(x = (x_min + x_max)/2, y = y_star, label = signif), 
                inherit.aes = FALSE, size = 6, fontface = "bold", color = "black", vjust = 0) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  }
  
  print(p)
  clean_metric_name <- gsub("[ /]", "_", metric_) 
  ggsave(paste0(clean_metric_name, '_r', r, '.pdf'), plot = p, 
         width = 9, height = 4, units = "in", dpi = 300)
}
