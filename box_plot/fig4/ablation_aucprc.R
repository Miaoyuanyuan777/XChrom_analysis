library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

df <- read_csv('ablation_rna_seq_data.csv')

df_auc <- df %>%
  filter(grepl("auROC|auPRC", Metric))

df_auc$Metric <- factor(df_auc$Metric, levels = c("auROC", "per cell auROC", "per peak auROC", 
                                                  "auPRC", "per cell auPRC", "per peak auPRC"))

df_auc$Test <- factor(df_auc$Test, levels = c("cross_cell", "cross_region", "cross_both"))

studies_levels <- c("RAW_atac", 
                    "RNA_shuffle", "RNA_random_easier","RNA_random","RNA_random_complex",
                    "RNA_embed","RNA_embed_dense",
                    "X_pca_z_layernorm", "X_pca_nodense", "X_pca_complex",
                    "SEQ_embed","SEQ_sei","SEQ_easier","RNA_scVI", "RAW_XChrom")
df_auc$Studies <- factor(df_auc$Studies, levels = studies_levels)

df_auc <- df_auc %>% filter(!is.na(Metric) & !is.na(Studies) & !is.na(Test))

df_auc <- df_auc %>%
  mutate(MetricType = ifelse(grepl("auPRC", Metric), "auPRC", "auROC")) %>%
  mutate(MetricType = factor(MetricType, levels = c("auROC", "auPRC")))

ggplot(df_auc, aes(x = Metric, y = Value, fill = Studies)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, size = 0.3) +
  scale_fill_manual(
    values = c(
      "RAW_XChrom"  = "#FF7F00", 
      # "RAW_atac"    = "#999999", 
      "X_pca_nodense"     =  "#C6DBEF",
      # "X_pca_z_layernorm" = "#99D8C9",
      "X_pca_complex" = "#6BAED6",
      "RNA_scVI"    = "#0066cc", 
      # "RNA_shuffle" = "#6BAED6", 
      # "RNA_random"  = "#C6DBEF",
      # "RNA_random_easier" = "#0066cc",
      # "RNA_random_complex" = "#33ccff",
      "RNA_embed" = "#33ccff",
      # "RNA_embed_dense" = "#31A354",
      "SEQ_embed"  = "#CAB2D6",
      "SEQ_sei" = "#8E62B5",# "#6A3D9A"
      "SEQ_easier" = "#5B4B8A"
    )
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  
  facet_grid(Test ~ MetricType, scales = "free") +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black", margin = margin(t = 2)),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
    legend.text = element_text(size = 14),
    legend.position = "bottom", 
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3), color = "black", size = 0.5), 
    strip.text = element_text(size = 14, face = "bold", color = "black"),  
    strip.text.x = element_text(size = 14, face = "bold", color = "black"),  
    panel.border = element_rect(size = 0.5, color = "black")
  ) +
  labs(title = "XChrom Ablation Studies Performance", x = "", y = "Value")

ggsave("ablation_auc_cross3task.pdf", plot = last_plot(), 
       width = 12, height = 12, units = "in", dpi = 300)


#########################
#########################
# plot
####################

df <- read_csv('ablation_rna_seq_data.csv')

df <- df %>% 
  mutate(across(c(Studies, Test, Metric), str_trim))

df_auc <- df %>%
  filter(grepl("auROC|auPRC", Metric))

df_auc$Metric <- factor(df_auc$Metric, levels = c("auROC", "per cell auROC", "per peak auROC", 
                                                  "auPRC", "per cell auPRC", "per peak auPRC"))

df_auc$Test <- factor(df_auc$Test, levels = c("cross_cell", "cross_region", "cross_both"))

df_auc <- df_auc %>% filter(!is.na(Metric) & !is.na(Studies) & !is.na(Test))

df_auc <- df_auc %>%
  mutate(MetricType = ifelse(grepl("auPRC", Metric), "auPRC", "auROC")) %>%
  mutate(MetricType = factor(MetricType, levels = c("auROC", "auPRC")))

study_colors <- c(
  "RAW_XChrom"  = "#FF7F00", 
  # "RAW_atac"    = "#999999", 
  "X_pca_nodense"     =  "#C6DBEF",
  # "X_pca_z_layernorm" = "#99D8C9",
  "X_pca_complex" = "#6BAED6",
  "RNA_scVI"    = "#0066cc", 
  # "RNA_shuffle" = "#6BAED6", 
  # "RNA_random"  = "#C6DBEF",
  # "RNA_random_easier" = "#0066cc",
  # "RNA_random_complex" = "#33ccff",
  "RNA_embed" = "#33ccff",
  # "RNA_embed_dense" = "#31A354",
  "SEQ_embed"  = "#CAB2D6",
  "SEQ_sei" = "#8E62B5",# "#6A3D9A"
  "SEQ_easier" = "#5B4B8A"
)

plot_ablation_group <- function(data, target_studies, output_name) {
  
  sub_data <- data %>% filter(Studies %in% target_studies)
  
  sub_data$Studies <- factor(sub_data$Studies, levels = target_studies)
  
  p <- ggplot(sub_data, aes(x = Metric, y = Value, fill = Studies)) +
    geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, size = 0.3) +
    scale_fill_manual(values = study_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
    facet_grid(Test ~ MetricType, scales = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold", color = "black", margin = margin(t = 2)),
      axis.text.y = element_text(size = 14, face = "bold", color = "black"),
      axis.title = element_text(size = 16, color = "black", face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 10)),
      legend.text = element_text(size = 14),
      legend.position = "bottom", 
      strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3), color = "black", size = 0.5), 
      strip.text = element_text(size = 14, face = "bold", color = "black"),  
      strip.text.x = element_text(size = 14, face = "bold", color = "black"),  
      panel.border = element_rect(size = 0.5, color = "black")
    ) +
    labs(title = paste("XChrom Ablation Studies -", output_name), x = "", y = "Value")
  print(p)
}

# SEQ_random，RAW_XChrom，RAW_atac
plot_ablation_group(df_auc, 
                    c("RAW_XChrom",
                      "SEQ_embed",
                      "SEQ_easier",
                      "SEQ_sei"), 
                    "group5_seq")

# ggsave("ablation_dna_aucprc.pdf", width = 12, height = 9, units = "in", dpi = 300)


# X_pca_nodense,X_pca_complex,RNA_embed,RNA_scVI,RNA_shuffle
plot_ablation_group(df_auc, 
                    c("RAW_XChrom",
                      "X_pca_nodense", 
                      "X_pca_complex", 
                      "RNA_embed",
                      "RNA_scVI"
                      ), 
                    "group1_pca")

#######################################
########  auROC,auPRC 
#######################################

library(dplyr)
library(cowplot)
df <- read_csv('ablation_rna_seq_data.csv')
df <- df %>% mutate(across(c(Studies, Test, Metric), str_trim))
df_auc <- df %>% filter(grepl("auROC|auPRC", Metric))
df_auc$Metric <- factor(df_auc$Metric, levels = c("auROC", "per cell auROC", "per peak auROC", 
                                                  "auPRC", "per cell auPRC", "per peak auPRC"))
df_auc$Test <- factor(df_auc$Test, levels = c("cross_cell", "cross_region", "cross_both"))
df_auc <- df_auc %>% filter(!is.na(Metric) & !is.na(Studies) & !is.na(Test))
df_auc <- df_auc %>%
  mutate(MetricType = ifelse(grepl("auPRC", Metric), "auPRC", "auROC")) %>%
  mutate(MetricType = factor(MetricType, levels = c("auROC", "auPRC")))

study_colors <- c(
  "RAW_XChrom"  = "#FF7F00", 
  # "RAW_atac"    = "#999999", 
  "X_pca_nodense"     =  "#C6DBEF",
  # "X_pca_z_layernorm" = "#99D8C9",
  "X_pca_complex" = "#6BAED6",
  "RNA_scVI"    = "#0066cc", 
  # "RNA_shuffle" = "#6BAED6", 
  # "RNA_random"  = "#C6DBEF",
  # "RNA_random_easier" = "#0066cc",
  # "RNA_random_complex" = "#33ccff",
  "RNA_embed" = "#33ccff",
  # "RNA_embed_dense" = "#31A354",
  "SEQ_embed"  = "#CAB2D6",
  "SEQ_sei" = "#8E62B5",# "#6A3D9A"
  "SEQ_easier" = "#5B4B8A"
)

library(ggplot2)
library(dplyr)
library(patchwork)

my_get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  leg_idx <- which(sapply(g$grobs, function(x) grepl("guide-box", x$name)))
  if (length(leg_idx) > 0) {
    return(g$grobs[[leg_idx[1]]])
  }
  return(NULL)
}

plot_ablation_group <- function(data, target_studies, output_name) {
  
  sub_data <- data %>% filter(Studies %in% target_studies)
  sub_data$Studies <- factor(sub_data$Studies, levels = target_studies)
  
  my_theme <- theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold", color = "black"),
      axis.text.y = element_text(size = 14, face = "bold", color = "black"),
      axis.title = element_text(size = 16, color = "black", face = "bold"),
      legend.text = element_text(size = 14),
      
      strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3), color = "black", linewidth = 0.5), 
      strip.text = element_text(size = 14, face = "bold", color = "black"),  
      strip.text.x = element_text(size = 14, face = "bold", color = "black"),  
      panel.border = element_rect(linewidth = 0.5, color = "black")
    )
  
  # 1. auROC
  p_roc <- ggplot(sub_data %>% filter(MetricType == "auROC"), aes(x = Metric, y = Value, fill = Studies)) +
    geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, linewidth = 0.3) +
    scale_fill_manual(values = study_colors, name = NULL) + 
    facet_grid(Test ~ MetricType, scales = "free_x") +
    coord_cartesian(ylim = c(0.45, 0.85)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    my_theme +
    labs(x = NULL, y = "auROC Value")
  
  # 2. auPRC
  p_prc <- ggplot(sub_data %>% filter(MetricType == "auPRC"), aes(x = Metric, y = Value, fill = Studies)) +
    geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, linewidth = 0.3) +
    scale_fill_manual(values = study_colors, name = NULL) + 
    facet_grid(Test ~ MetricType, scales = "free_x") +
    coord_cartesian(ylim = c(0.2, 0.6)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    my_theme +
    labs(x = NULL, y = "auPRC Value")
  

  shared_legend <- my_get_legend(
    p_roc + theme(legend.position = "bottom")
  )
  p_roc <- p_roc + theme(legend.position = "none")
  p_prc <- p_prc + theme(legend.position = "none")
  plot_row <- plot_grid(p_roc, p_prc, ncol = 2, align = "h", axis = "tb")
  title_gg <- ggdraw() + 
    draw_label(paste("XChrom Ablation Studies -", output_name), 
               fontface = 'bold', size = 20, hjust = 0.5)
  
  p_combined <- plot_grid(
    title_gg,
    plot_row,
    shared_legend,
    ncol = 1,
    rel_heights = c(0.08, 1, 0.1) 
  )
  
  print(p_combined)
}
# SEQ_random，RAW_XChrom，RAW_atac
plot_ablation_group(df_auc,
                    c("SEQ_embed","SEQ_easier","SEQ_sei","RAW_XChrom"),
                    "group5_seq")

ggsave("ablation_dna_aucprc.pdf", width = 12, height = 9, units = "in", dpi = 300)

plot_ablation_group(df_auc, 
                    c("RNA_scVI",
                      "RNA_embed",
                      "X_pca_nodense", 
                      "X_pca_complex", 
                      "RAW_XChrom"
                      ), 
                    "group1_pca")
ggsave("ablation_rna_aucprc.pdf", width = 12, height = 9, units = "in", dpi = 300)