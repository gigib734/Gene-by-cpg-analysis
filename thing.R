methylation_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DMPs blood PCDH genes - Meta-analysis results DMPs blood PCDH genes (1).csv")
transcription_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DEGs blood PCDH genes - Meta-analysis results DEGs blood PCDH genes.csv.csv")
library(dplyr)
library(tidyverse)
library(readr)
library(ggplot2)
A1_filtered_meth <- methylation_data[grepl("(^|;)PCDHGA1(;|$)", methylation_data$Annotated_genes), ]
A2_filtered_meth <- methylation_data[grepl("(^|;)PCDHGA2(;|$)", methylation_data$Annotated_genes), ]
A3_filtered_meth <- methylation_data[grepl("(^|;)PCDHGA3(;|$)", methylation_data$Annotated_genes), ]
A4_filtered_meth <- methylation_data[grepl("(^|;)PCDHGA4(;|$)", methylation_data$Annotated_genes), ]
B1_filtered_meth <- methylation_data[grepl("(^|;)PCDHGB1(;|$)", methylation_data$Annotated_genes), ]

# Using base R
A1234B1_dif_cpgs <- methylation_data[grepl("PCDHGA1", methylation_data$annotated_genes) + 
                          grepl("PCDHGA2", methylation_data$annotated_genes) + 
                          grepl("PCDHGA3", methylation_data$annotated_genes) + 
                          grepl("PCDHGA4", methylation_data$annotated_genes) + 
                          grepl("PCDHGB1", methylation_data$annotated_genes) == 1, ]

methylation_data$Effect_size <- as.numeric(methylation_data$Effect_size)

library(dplyr)
library(ggplot2)

# Add a gene category column to each table
A1_filtered_meth$Gene <- "PCDHGA1"
A2_filtered_meth$Gene <- "PCDHGA2"
A3_filtered_meth$Gene <- "PCDHGA3"
A4_filtered_meth$Gene <- "PCDHGA4"
B1_filtered_meth$Gene <- "PCDHGB1"

# Combine all tables
combined_table <- bind_rows(A1_filtered_meth, A2_filtered_meth, A3_filtered_meth, 
                            A4_filtered_meth, B1_filtered_meth)

# Add y-values based on gene category
combined_table <- combined_table %>%
  mutate(y_value = case_when(
    Gene == "PCDHGA1" ~ 0.0005,
    Gene == "PCDHGA2" ~ 0.0002,
    Gene == "PCDHGA3" ~ 0.0007,
    Gene == "PCDHGA4" ~ 0.0022,
    Gene == "PCDHGB1" ~ 0.0043
  ))

# Add a column to identify shared CPGs (appearing in all 5 genes)
combined_table <- combined_table %>%
  group_by(CpG) %>%  # Replace 'CpG_column' with your actual CPG column name
  mutate(shared = n_distinct(Gene) == 5) %>%
  ungroup() %>%
  mutate(color_group = ifelse(shared, Gene, "Unique"),
         label_text = ifelse(shared, "", Gene))  # Only label if unique

# Add a column to identify shared CPGs (appearing in all 5 genes)
combined_table <- combined_table %>%
  group_by(CpG) %>%  # Replace 'CpG_column' with your actual CPG column name
  mutate(shared = n_distinct(Gene) == 5) %>%
  ungroup() %>%
  mutate(color_group = ifelse(shared, Gene, "Unique"),
         label_text = ifelse(shared, "", CpG))  # Label with CPG name if unique

# Plot with conditional coloring and labeling
ggplot(combined_table, aes(x = Effect_size, y = y_value, color = color_group)) +
  geom_point(size = 3) +
  geom_text(aes(label = label_text), vjust = -0.5, size = 3) +
  labs(x = "Effect Size",
       y = "Transcription Effect",
       title = "CPGs by Effect Size",
       color = "Gene") +
  theme_minimal() +
  scale_color_manual(values = c("PCDHGA1" = "red", 
                                "PCDHGA2" = "blue", 
                                "PCDHGA3" = "green", 
                                "PCDHGA4" = "purple", 
                                "PCDHGB1" = "orange",
                                "Unique" = "pink"))

ggsave("combined_plot.png", width = 20, height = 6, units = "in")
