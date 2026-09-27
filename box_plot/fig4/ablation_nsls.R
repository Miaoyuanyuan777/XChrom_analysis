library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

df <- read_csv('ablation_rna_seq_data.csv')
df <- df %>% 
  mutate(across(c(Studies, Test, Metric), str_trim))
df_nsls <- df %>%
  filter(grepl("^ls\\(|^ns\\(", Metric))
df_nsls$Metric <- factor(df_nsls$Metric, levels = c("ns(10)", "ns(50)", "ns(100)",
                                                    "ls(10)", "ls(50)", "ls(100)"))
df_nsls$Test <- factor(df_nsls$Test, levels = c("test_cells", "denoise"))
studies_levels <- c("RAW_atac", 
                    "RNA_shuffle", "RNA_random_easier","RNA_random","RNA_random_complex",
                    "RNA_embed","RNA_embed_dense",
                    "X_pca_z_layernorm", "X_pca_nodense", "X_pca_complex",
                    "SEQ_embed","SEQ_sei","SEQ_easier","RNA_scVI", "RAW_XChrom")
df_nsls$Studies <- factor(df_nsls$Studies, levels = studies_levels)

df_nsls <- df_nsls %>% filter(!is.na(Metric) & !is.na(Studies) & !is.na(Test))
df_nsls <- df_nsls %>%
  mutate(MetricType = ifelse(grepl("^ls", Metric), "Label Score (ls)", "Neighborhood Score (ns)")) %>%
  mutate(MetricType = factor(MetricType, levels = c("Neighborhood Score (ns)","Label Score (ls)")))

ggplot(df_nsls, aes(x = Metric, y = Value, fill = Studies)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, size = 0.3) +
  
  scale_fill_manual(
    values = c(
      "RAW_XChrom"  = "#FF7F00", 
      # "RAW_atac"    = "#999999", 
      "X_pca_nodense"     =  "#C6DBEF",
      # "X_pca_z_layernorm" = "#99D8C9",
      "X_pca_complex" = "#6BAED6",
      "RNA_scVI"    = "#99D8C9", 
      # "RNA_shuffle" = "#6BAED6", 
      # "RNA_random"  = "#C6DBEF", 
      # "RNA_random_easier" = "#0066cc",
      # "RNA_random_complex" = "#33ccff",
      "RNA_embed" = "#1C9099",
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
  labs(title = "XChrom Ablation Studies - NS/LS Performance", x = "", y = "Value")

# 保存图片
# ggsave("ablation_nsls_barplot.pdf", plot = last_plot(), 
#        width = 12, height = 10, units = "in", dpi = 300)


#########################
#########################
# plot
####################
df <- read_csv('ablation_rna_seq_data.csv')
df <- df %>% 
  mutate(across(c(Studies, Test, Metric), str_trim))
df_nsls <- df %>%
  filter(grepl("^ls\\(|^ns\\(", Metric))
df_nsls$Metric <- factor(df_nsls$Metric, levels = c("ns(10)", "ns(50)", "ns(100)",
                                                    "ls(10)", "ls(50)", "ls(100)"))
df_nsls$Test <- factor(df_nsls$Test, levels = c("test_cells", "trainval_cells"))
df_nsls <- df_nsls %>% filter(!is.na(Metric) & !is.na(Studies) & !is.na(Test))
df_nsls <- df_nsls %>%
  mutate(MetricType = ifelse(grepl("^ls", Metric), "Label Score (ls)", "Neighbor Score (ns)")) %>%
  mutate(MetricType = factor(MetricType, levels = c("Neighbor Score (ns)","Label Score (ls)")))
study_colors <- c(
  "RAW_XChrom"  = "#FF7F00", 
  "RAW_atac"    = "#999999",
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
    geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, linewidth = 0.3) +
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

plot_ablation_group(df_nsls,
                    c( "RAW_atac","SEQ_embed","SEQ_easier","SEQ_sei","RAW_XChrom"),
                    "group5_seq")

plot_ablation_group(df_nsls, 
                    c("RAW_atac",
                      "RNA_scVI",
                      "RNA_embed",
                      "X_pca_nodense", 
                      "X_pca_complex", 
                      "RAW_XChrom"
                      # "RNA_shuffle"
                    ), 
                    "group1_pca")

#######################################
########  ns,ls
#######################################
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
  if (nrow(sub_data) == 0) {
    stop(paste("\n❌ error", 
               paste(target_studies, collapse = ", "), 
               "\nNo matching studies found in the data!"))
  }
  
  sub_data$Studies <- factor(sub_data$Studies, levels = target_studies)
  
  ns_data <- sub_data %>% filter(grepl("ns", MetricType, ignore.case = TRUE))
  ls_data <- sub_data %>% filter(grepl("ls", MetricType, ignore.case = TRUE))
  if (nrow(ns_data) == 0) stop("❌ error: No rows found containing 'ns' in MetricType!")
  if (nrow(ls_data) == 0) stop("❌ error: No rows found containing 'ls' in MetricType!")
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
  p_ns <- ggplot(ns_data, aes(x = Metric, y = Value, fill = Studies)) +
    geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, linewidth = 0.3) +
    scale_fill_manual(values = study_colors, name = NULL) + 
    facet_grid(Test ~ MetricType, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    my_theme +
    labs(x = NULL, y = "NS Value")
  p_ls <- ggplot(ls_data, aes(x = Metric, y = Value, fill = Studies)) +
    geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, linewidth = 0.3) +
    scale_fill_manual(values = study_colors, name = NULL) + 
    facet_grid(Test ~ MetricType, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    my_theme +
    labs(x = NULL, y = "LS Value")
  

  shared_legend <- my_get_legend(p_ns + theme(legend.position = "bottom"))
  p_ns <- p_ns + theme(legend.position = "none")
  p_ls <- p_ls + theme(legend.position = "none")
  plot_row <- plot_grid(p_ns, p_ls, ncol = 2, align = "h", axis = "tb")
  title_gg <- ggdraw() + draw_label(paste("XChrom Ablation Studies -", output_name), fontface = 'bold', size = 20, hjust = 0.5)
  
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
plot_ablation_group(df_nsls,
                    c( "RAW_atac","SEQ_embed","SEQ_easier","SEQ_sei","RAW_XChrom"),
                    "group5_seq")
ggsave("ablation_dna_nsls.pdf", width = 12, height = 9, units = "in", dpi = 300)
plot_ablation_group(df_nsls, 
                    c("RAW_atac",
                      "RNA_scVI",
                      "RNA_embed",
                      "X_pca_nodense", 
                      "X_pca_complex", 
                      "RAW_XChrom"
                      # "RNA_shuffle"
                    ), 
                    "group1_pca")
ggsave("ablation_rna_nsls.pdf", width = 12, height = 9, units = "in", dpi = 300)