library(tidyverse)
library(ggplot2)
library(ggh4x)
library(scales) 

df <- read_csv("disp_metrics.csv")
colnames(df) <- c("Ratio", "Context", "Metric", "Value")

df_clean <- df %>%
  mutate(
    Ratio_Clean = case_when(
      grepl("raw", Ratio, ignore.case = TRUE) ~ "raw",
      grepl("raw", Context, ignore.case = TRUE) ~ "raw",
      TRUE ~ as.character(Ratio)
    ),
    Context_Clean = case_when(
      grepl("b_cells", Context, ignore.case = TRUE) ~ "test_b_cells",
      grepl("other_cells", Context, ignore.case = TRUE) ~ "test_other_cells",
      grepl("test_cells", Context, ignore.case = TRUE) ~ "test_cells",
      TRUE ~ Context
    )
  )

# # raw -> 0 -> 0.01 -> 0.02 -> 0.1
other_ratios <- setdiff(unique(df_clean$Ratio_Clean), "raw")
other_ratios_sorted <- other_ratios[order(as.numeric(other_ratios))]
final_levels <- c("raw", other_ratios_sorted)

df_clean$Ratio_Clean <- factor(df_clean$Ratio_Clean, levels = final_levels)


############
# auroc,auprc, test_b_cells, test_other_cells
############

df_p1 <- df_clean %>%
  filter(grepl("auROC|auPRC", Metric)) %>%
  filter(grepl("test_b_cells|test_other_cells", Context)) %>%
  mutate(
    MetricType = ifelse(grepl("auROC", Metric), "auROC", "auPRC"),
    MetricType = factor(MetricType, levels = c("auROC", "auPRC")),
    Scope = case_when(
      grepl("per cell", Metric) ~ "per cell",
      grepl("per peak", Metric) ~ "per peak",
      TRUE ~ "overall"
    ),
    Scope = factor(Scope, levels = c("overall", "per cell", "per peak"))
  )

ggplot(df_p1, aes(x = Ratio_Clean, y = Value, color = Context_Clean, group = Context_Clean)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_grid(MetricType ~ Scope, scales = "free_y") +
  scale_color_brewer(palette = "Set1", name = "Cell Context") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "#EFEFEF"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "auROC and auPRC Trends across Ratios",
    x = "Ratio",
    y = "Score Value"
  )

ggsave("B_marker_peaks_trends_aucprc.pdf", width = 8, height = 6, dpi = 300)

#########################################
##### k=10
##########################################

df_p5 <- df_clean %>%
  filter(Metric %in% c("ns(10)", "ls(10)")) %>%
  filter(Context %in% c("test_b_cells", "test_other_cells", "raw_test_b_cells", "raw_test_other_cells")) %>%
  mutate(
    Metric_Category = ifelse(grepl("^ns", Metric), "Neighbor Score (NS)", "Label Score (LS)"),
    Cell_Group = ifelse(grepl("b_cells", Context), "test_b_cells", "test_other_cells")
  )

df_lines <- df_p5 %>% 
  filter(!grepl("raw", Context))

df_raw <- df_p5 %>%
  filter(grepl("raw", Context)) %>%
  group_by(Metric_Category, Cell_Group) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

df_raw_expanded <- df_raw %>%
  expand_grid(Ratio = unique(df_lines$Ratio))




ggplot() +
  geom_line(data = df_lines, 
            aes(x = as.factor(Ratio), y = Value, color = Cell_Group, group = Cell_Group, linetype = "Test Data"), 
            linewidth = 1) + 
  geom_point(data = df_lines, 
             aes(x = as.factor(Ratio), y = Value, color = Cell_Group, group = Cell_Group), 
             size = 2.5) + 
  
  geom_line(data = df_raw_expanded, 
            aes(x = as.factor(Ratio), y = Value, color = Cell_Group, group = Cell_Group, linetype = "Baseline (raw_ATAC)"), 
            linewidth = 0.8) +
  
  facet_wrap(~ Metric_Category, scales = "free_y") +
  
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("Test Data" = "solid", "Baseline (raw_ATAC)" = "dashed")) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "#EFEFEF"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Performance of NS(10) and LS(10) over Ratio",
    x = "Ratio",
    y = "Score Value",
    color = "Cell Context",  
    linetype = "Data Type"
  )

ggsave("B_marker_peaks_trends_nsls10.pdf", width = 8, height = 5, dpi = 300)

#####################################
##### Jaccard index，B specific peaks
#######################################
# ==========================================
# B_marker_peaks 
# ==========================================

df_b_marker <- df_clean %>%
  filter(grepl("B_marker_peaks", Context, ignore.case = TRUE)) 


p_b_marker <- ggplot(df_b_marker, aes(x = Ratio_Clean, y = Value, color = Metric, group = Metric)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 0, hjust = 1, face = "bold"), 
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Metric Trends across Ratios (Context: B_marker_peaks)",
    x = "Ratio",
    y = "Value",
    color = "Metric"
  )

print(p_b_marker)
ggsave("B_marker_peaks_trends_overlap.pdf",width = 4, height = 5, dpi = 300)