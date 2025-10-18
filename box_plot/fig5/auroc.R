library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
df <- read_csv("auroc.csv")

# make model factor
df$Model <- factor(df$Model, levels = c("human_model", "mouse_model"))
df$Species<-factor(df$Species,levels=c('human','mouse','macaque','marmoset'))

##############----------auROC--------#########
filtered_df <- df %>% filter(Metric == "per cell auROC")  
## per cell auROC
## per peak auROC


#########----------filt human & mouse which are not paired-------##########
filtered_df <- filtered_df %>% 
  group_by(Species, Samples) %>%
  filter(n() == 2)
################################################################################
y_pos <- filtered_df %>% 
  group_by(Species) %>% 
  summarise(max_val = max(Value)) %>% 
  mutate(y_position = max_val + 0.02) 

annotation_df <- filtered_df %>%
  pivot_wider(names_from = Model, values_from = Value) %>%
  mutate(diff =  mouse_model - human_model) %>%
  group_by(Species) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    p_value = wilcox.test(human_model, mouse_model, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_value_label = sprintf("p = %.3f", p_value)
  )
annotation_df <- left_join(annotation_df, y_pos, by = "Species")

ggplot(filtered_df, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(width = 0.6, alpha = 0.9, outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) + 
  facet_wrap(~ Species, scales = "free_x", strip.position = "bottom",nrow=1) + 
  scale_fill_manual(values = c("human_model" = "orange", 
                               "mouse_model" = "#7B4F94")) +
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
  labs(title = "", x = "", y = "per cell auROC") +
  geom_text(
    data = annotation_df[annotation_df$p_value<0.05,],
    aes(x = 1.5, y = y_position,  
        label = paste0(p_value_label, "\n95%CI = [", 
                       round(lower_ci, 2), ",", 
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4,
    color = "black",hjust = 0.5  
  )

ggsave("per-cell-auROC.pdf", plot = last_plot(), 
       width = 9, height = 4, units = "in",dpi=300)
