library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
df <- read_csv('intra_method_auc.csv')

# Convert Metric & Method & Dataset to factor and specify the order
df$Metric <- factor(df$Metric, levels = c("auROC", "per cell auROC", "per peak auROC", 
                                          "auPRC", "per cell auPRC", "per peak auPRC"))
df$Method <- factor(df$Method, levels = c("Seurat", "LS_Lab", "MultiVI", "XChrom"))
df$Dataset<-factor(df$Dataset, levels = c("h_pbmc","m_brain","h_brain",'h_gonads','m_palates'))

# Add 2 new column to group auROC(overall,per-cell,per-peak) and auPRC(overall,per-cell,per-peak),respectively
df <- df %>%
  mutate(MetricType = ifelse(grepl("auPRC", Metric), "auPRC", "auROC"))%>%
  mutate(MetricType = factor(MetricType, levels = c("auROC", "auPRC")))

# Obtain all the methods that need to be compared (excluding XChrom itself)
other_methods <- unique(df$Method[df$Method != "XChrom"])

# Extract the data of XChrom as the benchmark
xchrom_data <- df %>% 
  filter(Method == "XChrom") %>% 
  select(Dataset, Metric, MetricType, Value_XChrom = Value)

# Generate comparative data
comparison_data <- df %>%
  filter(Method != "XChrom") %>% 
  inner_join(xchrom_data, by = c("Dataset", "Metric", "MetricType")) %>%
  mutate(diff = Value_XChrom - Value)

# Calculate statistical results
annotation_df <- comparison_data %>%
  group_by(Metric, Method) %>%  # # Group by Metric & Test
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_diff - qt(0.975, df = n()-1) * se_diff,
    upper_ci = mean_diff + qt(0.975, df = n()-1) * se_diff,
    p_value = t.test(Value_XChrom, Value, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_value_label = ifelse(p_value < 0.05,
                           paste0("p = ", signif(p_value, digits = 3)),
                           NA_character_)
  ) %>%
  # Merge the MetricType and Y-axis position information
  left_join(distinct(df, Metric, MetricType), by = "Metric") %>% 
  left_join(
    df %>% group_by(Metric) %>% summarise(max_val = max(Value)),
    by = "Metric"
  ) %>% 
  mutate(y_position = ifelse(p_value < 0.05, max_val + 0.15, NA_real_))

ggplot(df, aes(x = Metric, y = Value,)) +
  geom_boxplot(aes(group = interaction(Metric, Method), fill = Method), 
               alpha = 1, outlier.shape = NA) +  
  geom_jitter(
    aes(group = interaction(Metric, Method), shape = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.2,   
      dodge.width = 0.8    
    ),size = 3,alpha = 1,color = "black",#stroke = 0.8
    )+
  scale_fill_manual(values = c("Seurat" = "skyblue", 'LS_Lab' = "pink",'MultiVI'='#CC99E6',"XChrom" = "#FFB07C")) + 
  scale_shape_manual(
    values = c(
      "h_pbmc" = 1,   
      "m_brain" = 2,    
      "h_brain" = 5,    
      "h_gonads" = 0,   
      "m_palates" = 6    
    ))+
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +  
  guides(shape = guide_legend(order = 1), fill = guide_legend(order = 2)) +  
  facet_wrap(~MetricType, scales = "free")+
  theme_set(theme_bw())+
  # theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16, face = "bold", color = "black", margin = margin(t =2)),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold",margin = margin(b = 10)),
    legend.text = element_text(size = 16),
    strip.background = element_rect(fill = rgb(0.5, 0.5, 0.5, 0.3),color = "black", size = 0.5),  
    strip.text = element_text(size = 16, face = "bold", color = "black"),  
    strip.text.x = element_text(size = 16, face = "bold",color = "black"), 
  ) +
  labs(title = "Cross-cell Prediction", x = "", y = "Value") +
  geom_text(
    data = annotation_df[annotation_df$p_value<0.05,],
    aes(x = Metric, y = y_position,group = MetricType,  
        label = paste0(p_value_label, "\n95%CI = [",
                       round(lower_ci, 2), ",",
                       round(upper_ci, 2), "]")),
    inherit.aes = FALSE,size = 4,
    color = "black",hjust = 0.5)


ggsave("intra_method_auc_box.pdf", plot = last_plot(), 
       width = 20, height = 6, units = "in",dpi=300)
