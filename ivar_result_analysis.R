
library(tidyverse)

# ---------------------------------------------------------------
# Import All Annotated iVar TSV Files
# ---------------------------------------------------------------
# Get list of all annotated TSV files
tsv_files <- list.files("ivar_results/annotated", pattern = "*_annotated.tsv$", full.names = TRUE)

# Read all files into one master data frame
master_df <- tsv_files %>%
  map_dfr(function(file) {
    # Extract sample ID from file name
    sample_name <- gsub("_annotated.tsv", "", basename(file))
    
    # Read TSV and add Sample column
    read_tsv(file, show_col_types = FALSE) %>%
      mutate(Sample = sample_name)
  })

# ---------------------------------------------------------------
# Filter Passing Variants & Assign Cohorts
# ---------------------------------------------------------------
clean_df <- master_df %>%
  filter(PASS == TRUE) %>% 
  mutate(Cohort = case_when(
    Sample %in% c("DRR067444", "DRR067448", "DRR067449") ~ "Mild",
    Sample %in% c("DRR067462", "DRR067463", "DRR067465") ~ "Severe",
    TRUE ~ "Unknown"
  ))

# View the first few rows of your clean data
head(clean_df)

# ---------------------------------------------------------------
# Variant Count Summary Table
# ---------------------------------------------------------------
gene_summary <- clean_df %>%
  group_by(Cohort, FEATURE) %>%
  summarise(
    Total_Variants = n(),
    Average_Allele_Freq = mean(ALT_FREQ),
    .groups = "drop"
  )
gene_summary_data <- clean_df %>%
  group_by(Cohort, FEATURE) %>%
  summarise(
    Total_Variants = n(),
    Average_Allele_Freq = round(mean(ALT_FREQ), 4),
    Min_Allele_Freq = round(min(ALT_FREQ), 4),
    Max_Allele_Freq = round(max(ALT_FREQ), 4),
    .groups = "drop"
  )

# Export directly to CSV
write_csv(gene_summary_data, "ivar_results/dengue_gene_summary.csv")

message("Successfully saved summary table to ivar_results/dengue_gene_summary.csv")
print("=== VARIANT COUNT BY GENE & COHORT ===")
print(gene_summary)

# ---------------------------------------------------------------
# Bar Chart (Variants per Dengue Protein)
# ---------------------------------------------------------------
plot1 <- ggplot(clean_df, aes(x = FEATURE, fill = Cohort)) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  labs(
    title = "Intra-Host Variant Burden Across Dengue 2 Proteins",
    x = "Mature Dengue Protein",
    y = "Number of Passing Variants",
    fill = "Cohort"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot to folder
ggsave("ivar_results/dengue_gene_variants.png", plot = plot1, width = 7, height = 4.5)

# ---------------------------------------------------------------
# Allele Frequency Across the Dengue Genome
# ---------------------------------------------------------------
plot2 <- ggplot(clean_df, aes(x = POS, y = ALT_FREQ, color = Cohort)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "Intra-Host Variant Allele Frequency Along Dengue Genome",
    x = "Genomic Position (bp)",
    y = "Alternate Allele Frequency (AF)",
    color = "Cohort"
  )

# Save plot to folder
ggsave("ivar_results/dengue_genome_af_scatter.png", plot = plot2, width = 8, height = 4.5)