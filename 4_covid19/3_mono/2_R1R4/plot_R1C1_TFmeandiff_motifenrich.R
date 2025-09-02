library(scales) 
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(tidyverse)
library(cowplot)
library(scales)

r1 <- read_csv("./1_differential_analysis/DAresults_R1vsOther.csv")
r1$pval_adj<-r1$padj
motif_r1 <- read_csv("./2_motif_enrichment/C1UP_diff005p_enrichmotif.csv")


# # Screen significant results
r1_sig <- r1 %>%
  filter(pval_adj < 0.05, abs(mean_difference) > 0.03)

merge_data <- function(rna_data, motif_data, cluster){
  inner_join(rna_data %>% 
               mutate(cluster = cluster),motif_data %>% 
               rename(tf = TF_prefix),by = "tf") %>%
    mutate(sig_stars = cut(mlog10padj,
                           breaks = c(-Inf, 1.3, Inf),
                           labels = c(' ',"*")))}

df_r1 <- merge_data(r1_sig, motif_r1, "R1")

## ##################################
r1_diffsort <- df_r1 %>%
  group_by(tf) %>%
  arrange(desc(mean_difference)) %>% 
  pull(tf)

r11 <- df_r1 %>% 
  mutate(tf = factor(tf, levels = r1_diffsort)) %>% 
  mutate(mlog10_pval_adj = -log10(pval_adj + 1e-300))

# Create graphics--------------------------------------------------------
# (1) The bar chart of the mean difference in TF activity above
p_mean_diff <- ggplot(r11, aes(x = tf)) +
  geom_col(aes(y = mean_difference, fill = mlog10_pval_adj),
           width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid",
             color = "grey", linewidth = 0.5) +
  scale_fill_gradient(
    low = "grey90", high = "red",
    limits = c(0, ceiling(max(r11$mlog10_pval_adj))),
    name = "-log10(Adjusted-Pval)",
    guide = guide_colorbar( 
      direction = "horizontal", 
      title.position = "top",  
      barwidth = unit(4, "cm")  
    )) +
  labs(y = "Activity mean difference\n(R1 vs Other)", title = "COVID-19:R1(ncMono)") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )
p_mean_diff

# (2) The bubble chart below shows the significance of motif enrichment
r11 <- r11 %>%
  mutate(fill_value = if_else(mlog10padj >= 1.3, mlog10padj, 0))

p_motif <- ggplot(r11, aes(x = tf)) +
  geom_point(aes(y = 0, fill = fill_value),
             size = 7, 
             shape = 21,
             alpha = 0.8) +
  scale_fill_gradientn(
    colours = c("grey80", "blue", "black"),
    values = scales::rescale(c(0, 1.3, ceiling(max(r11$mlog10padj)))),
    breaks = c(0,1.3,ceiling(max(r11$mlog10padj))),
    limits = c(0, ceiling(max(r11$mlog10padj))),
    name = "-log10(Adjusted-Pval)",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(4, "cm"))
  ) +
  labs(x = "Transcription factor (TF)") +
  scale_y_continuous(
    limits = c(-0.001, 0.001),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1.2,
      margin = margin(t = -6)
    ),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(-2,0,0,0), "cm"),
    aspect.ratio = 0.1,
    legend.position = "right"
  )
p_motif

# (3) Combine two graphics --------------------------------------------------------
aligned_plots <- align_plots(p_mean_diff, p_motif, align = "v")
final_plot <- plot_grid(
  aligned_plots[[1]], 
  aligned_plots[[2]],
  ncol = 1, 
  rel_heights = c(0.2, 0.4) 
)

final_plot

ggsave("./R1C1_bubble-bar.pdf", final_plot, 
       width = 16, height = 8, dpi = 300)