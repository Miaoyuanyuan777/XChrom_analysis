library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

# ============================================================
# 1. read data
# ============================================================

df <- read_csv("cross-region_auc.csv", show_col_types = FALSE) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c(
        "auROC", "per cell auROC", "per peak auROC",
        "auPRC", "per cell auPRC", "per peak auPRC"
      )
    ),
    Method = factor(Method, levels = c("scBasset", "XChrom")),
    Dataset = factor(
      Dataset,
      levels = c(
        "h_pbmc", "m_brain", "h_brain",
        "h_gonads", "m_palates"
      )
    ),
    MetricType = factor(
      if_else(grepl("auPRC", as.character(Metric)), "auPRC", "auROC"),
      levels = c("auROC", "auPRC")
    )
  )

if (anyNA(df[c("Metric", "Method", "Dataset", "Value")])) {
  stop("Metric、Method、Dataset、Value has NA values.")
}

if (any(duplicated(df[c("Dataset", "Metric", "Method")]))) {
  stop("Dataset–Metric–Method repeats found. Each combination should be unique.")
}

# ============================================================
# 2. paired statistical tests function
# ============================================================

perform_tests <- function(differences) {
  d <- differences[is.finite(differences)]
  n <- length(d)
  
  res <- list(
    n_nonzero = sum(d != 0),
    estimate_diff = NA_real_,
    estimate_type = NA_character_,
    lower_ci = NA_real_,
    upper_ci = NA_real_,
    effect_size = NA_real_,
    effect_size_type = NA_character_,
    test_method = "insufficient_n",
    p_value = NA_real_
  )
  
  if (n < 2L) return(res)
  
  mean_d <- mean(d)
  sd_d <- sd(d)
  
  if (all(d == 0)) {
    res$estimate_diff <- 0
    res$estimate_type <- "Mean_Difference"
    res$test_method <- "No observed difference"
    res$p_value <- 1
    return(res)
  }
  
  if (sd_d == 0) {
    bt <- binom.test(sum(d > 0), n, p = 0.5)
    
    res$estimate_diff <- mean_d
    res$estimate_type <- "Mean_Difference"
    res$test_method <- "Exact sign test"
    res$p_value <- bt$p.value
    return(res)
  }
  
  shapiro_p <- if (n >= 3L) {
    tryCatch(
      shapiro.test(d)$p.value,
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  
  # ----------------------------------------------------------
  # paired t-test
  # ----------------------------------------------------------
  
  if (is.finite(shapiro_p) && shapiro_p > 0.05) {
    tt <- tryCatch(
      t.test(d, mu = 0),
      error = function(e) NULL
    )
    
    if (!is.null(tt)) {
      res$estimate_diff <- unname(as.numeric(tt$estimate))
      res$estimate_type <- "Mean_Difference"
      res$lower_ci <- unname(tt$conf.int[1])
      res$upper_ci <- unname(tt$conf.int[2])
      res$effect_size <- mean_d / sd_d
      res$effect_size_type <- "Cohen_dz"
      res$test_method <- "Paired t-test"
      res$p_value <- tt$p.value
      
      return(res)
    }
  }
  
  # ----------------------------------------------------------
  # Wilcoxon signed-rank test
  # ----------------------------------------------------------
  
  wt <- suppressWarnings(
    tryCatch(
      wilcox.test(
        d,
        mu = 0,
        alternative = "two.sided",
        exact = FALSE,
        correct = TRUE,
        conf.int = TRUE
      ),
      error = function(e) {
        wilcox.test(
          d,
          mu = 0,
          alternative = "two.sided",
          exact = FALSE,
          correct = TRUE,
          conf.int = FALSE
        )
      }
    )
  )
  
  res$test_method <- "Wilcoxon signed-rank (normal approximation)"
  res$p_value <- wt$p.value
  res$estimate_type <- "Pseudo_median"
  
  if (
    !is.null(wt$estimate) &&
    length(wt$estimate) == 1L &&
    is.finite(wt$estimate)
  ) {
    res$estimate_diff <- unname(as.numeric(wt$estimate))
  }
  
  if (
    !is.null(wt$conf.int) &&
    length(wt$conf.int) == 2L &&
    all(is.finite(wt$conf.int))
  ) {
    res$lower_ci <- unname(wt$conf.int[1])
    res$upper_ci <- unname(wt$conf.int[2])
  }
  
  n_nonzero <- res$n_nonzero
  
  if (
    n_nonzero > 0L &&
    is.finite(wt$p.value) &&
    wt$p.value > 0
  ) {
    res$effect_size <- abs(qnorm(wt$p.value / 2)) / sqrt(n_nonzero)
    res$effect_size_type <- "Approx_Wilcoxon_r"
  }
  
  return(res)
}

# ============================================================
# 3. statistical tests for XChrom vs scBasset
# ============================================================

df_wide <- df %>%
  select(Dataset, Metric, MetricType, Method, Value) %>%
  pivot_wider(names_from = Method, values_from = Value)

if (
  !all(c("XChrom", "scBasset") %in% names(df_wide)) ||
  anyNA(df_wide[c("XChrom", "scBasset")])
) {
  stop("XChrom and scBasset not paired properly. Check the input data.")
}

pairwise_df <- df_wide %>%
  group_by(Metric, MetricType) %>%
  group_modify(~ {
    d <- .x$XChrom - .x$scBasset
    test_res <- perform_tests(d)
    
    tibble(
      Group1 = "XChrom",
      Group2 = "scBasset",
      N = sum(is.finite(d)),
      N_nonzero = test_res$n_nonzero,
      estimate_diff = test_res$estimate_diff,
      estimate_type = test_res$estimate_type,
      lower_ci = test_res$lower_ci,
      upper_ci = test_res$upper_ci,
      effect_size = test_res$effect_size,
      effect_size_type = test_res$effect_size_type,
      test_method = test_res$test_method,
      p_value = test_res$p_value
    )
  }) %>%
  ungroup() %>%
  mutate(
    signif = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write_csv(
  pairwise_df,
  "cross-region_XChrom_vs_scBasset_statistical_tests.csv"
)

print(pairwise_df, width = Inf)

# ============================================================
# 4. significant results for plotting
# ============================================================

sig_df <- pairwise_df %>%
  filter(!is.na(p_value), p_value < 0.05)

metric_indices <- df %>%
  distinct(MetricType, Metric) %>%
  arrange(MetricType, Metric) %>%
  group_by(MetricType) %>%
  mutate(metric_x = row_number()) %>%
  ungroup()

dodge_width <- 0.8
methods_list <- levels(df$Method)
k <- length(methods_list)

offsets <- seq(
  -dodge_width / 2 + dodge_width / (2 * k),
  dodge_width / 2 - dodge_width / (2 * k),
  length.out = k
)
names(offsets) <- methods_list

metric_max <- df %>%
  group_by(Metric) %>%
  summarise(max_val = max(Value), .groups = "drop")

plot_sig_df <- sig_df %>%
  left_join(metric_indices, by = c("MetricType", "Metric")) %>%
  left_join(metric_max, by = "Metric") %>%
  mutate(
    x1_offset = offsets[Group1],
    x2_offset = offsets[Group2],
    x_min = metric_x + pmin(x1_offset, x2_offset),
    x_max = metric_x + pmax(x1_offset, x2_offset)
  ) %>%
  group_by(MetricType, Metric) %>%
  mutate(
    y_pos = max_val + 0.05 + 0.04 * (row_number() - 1),
    y_tip = y_pos - 0.01
  ) %>%
  ungroup()

# ============================================================
# 5. plotting
# ============================================================

p <- ggplot(df, aes(x = Metric, y = Value)) +
  geom_boxplot(
    aes(group = interaction(Metric, Method), fill = Method),
    position = position_dodge(width = dodge_width),
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(group = interaction(Metric, Method), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = dodge_width
    ),
    size = 3,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "scBasset" = "#1A66B3",
      "XChrom" = "#FFB07C"
    )
  ) +
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
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  facet_wrap(~MetricType, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, size = 16,
      face = "bold", color = "black",
      margin = margin(t = 2)
    ),
    axis.text.y = element_text(
      size = 16, face = "bold", color = "black"
    ),
    axis.title = element_text(
      size = 16, color = "black", face = "bold"
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.title = element_text(
      hjust = 0.5, size = 20, face = "bold",
      margin = margin(b = 10)
    ),
    strip.background = element_rect(
      fill = rgb(0.5, 0.5, 0.5, 0.3),
      color = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(
      size = 16, face = "bold", color = "black"
    ),
    panel.border = element_rect(
      linewidth = 0.5, color = "black"
    )
  ) +
  labs(
    title = "Cross-region Prediction",
    x = "",
    y = "Value"
  )

if (nrow(plot_sig_df) > 0) {
  p <- p +
    geom_segment(
      data = plot_sig_df,
      aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos),
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = plot_sig_df,
      aes(x = x_min, xend = x_min, y = y_pos, yend = y_tip),
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = plot_sig_df,
      aes(x = x_max, xend = x_max, y = y_pos, yend = y_tip),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = plot_sig_df,
      aes(
        x = (x_min + x_max) / 2,
        y = y_pos + 0.005,
        label = signif
      ),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold"
    )
}

print(p)

ggsave(
  "cross-region_auc_box2.pdf",
  plot = p,
  width = 10,
  height = 6,
  units = "in"
)
