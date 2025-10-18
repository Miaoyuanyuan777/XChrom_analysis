library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

df <- read_csv("nsls.csv")

df$Model <- factor(df$Model, levels = c('raw_atac',"human_model", "mouse_model"))
df$Species <- factor(df$Species, levels = c('human','mouse','macaque','marmoset'))

#######################
filtered_df <- df %>% 
  filter(Metric == "ls(50)") %>% 
  group_by(Species, Samples) %>%
  filter(n() != 2)

y_pos_all <- filtered_df %>%
  group_by(Species) %>%
  summarise(max_val = max(Value)) %>%
  mutate(y_position = max_val + 0.08)

#######################
#1: human_model vs mouse_model
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
    p_value = wilcox.test(human_model, mouse_model, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(comparison = "human vs mouse", x_pos = 2.5) %>% 
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
    p_value = wilcox.test(raw_atac, mouse_model, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(comparison = "mouse vs raw", x_pos = 2.0) %>% 
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
    p_value = wilcox.test(raw_atac, human_model, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(comparison = "human vs raw", x_pos = 1.5) %>% 
  left_join(y_pos_all, by = "Species")

#######################
combined_annot <- bind_rows(annotation_df1, annotation_df2, annotation_df3) %>%
  filter(p_value < 0.05) %>%
  group_by(Species) %>%
  ungroup()

##############----------plot--------#########
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
  labs(x = NULL, y = "label score(k=10)") +
  geom_text(
    data = combined_annot,
    aes(x = x_pos, y = y_position, 
        label = paste0(p_value, "\n95%CI = [", 
                       round(lower_ci, 2), ",", 
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4, color = "black", vjust = 0.5
  )

ggsave("ls10.pdf", plot = last_plot(), 
       width = 9, height = 4, units = "in",dpi=300)
