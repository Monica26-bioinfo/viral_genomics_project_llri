## Bioinformatics Project LLRI 
## Intra-Host Genomic Variation Analysis of Dengue Virus Serotype 2: Mild vs. Severe Infection
## This is the Downstream analysis after generating the VCF files of the samples.

library(tidyverse)
library(vcfR)

# 1. Define Cohorts (Sample type)
sample_info <- tibble(
  sample_id = c("DRR067444", "DRR067448", "DRR067449", "DRR067462", "DRR067463", "DRR067465"),
  cohort    = c("Mild", "Mild", "Mild", "Severe", "Severe", "Severe")
)

# 2. Function to load and parse VCF files with accurate continuous AF
vcf_files <- list.files("variants", pattern = "\\.vcf$", full.names = TRUE)

parse_vcf <- function(vcf_path) {
  sample_name <- gsub("\\.vcf$", "", basename(vcf_path))
  
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  if (nrow(vcf@fix) == 0) return(NULL) # skip if empty
  
  fix_df <- getFIX(vcf) %>% as_tibble()
  
  # Extract AO (Alternate Count) and DP (Depth) matrices safely
  ao_matrix <- extract.gt(vcf, element = "AO", as.numeric = TRUE)
  dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  
  ao_val <- if (!is.null(ao_matrix)) as.numeric(ao_matrix[, 1]) else NA_real_
  dp_val <- if (!is.null(dp_matrix)) as.numeric(dp_matrix[, 1]) else NA_real_
  
  fix_df %>%
    transmute(
      sample_id = sample_name,
      POS       = as.numeric(POS),
      REF, ALT,
      QUAL      = as.numeric(QUAL),
      AO        = ao_val,
      DP        = dp_val,
      # Compute exact continuous AF (AO / DP); fallback to 1.0 if missing
      AF        = ifelse(!is.na(AO) & !is.na(DP) & DP > 0, AO / DP, 1.0)
    )
}

# Load raw variants
all_vars_raw <- map_dfr(vcf_files, parse_vcf) %>%
  left_join(sample_info, by = "sample_id")

cat("\n--- RAW UNFILTERED VARIANT COUNT PER SAMPLE ---\n")
print(table(all_vars_raw$sample_id))

# Apply quality filtering (QUAL >= 20)
all_vars <- all_vars_raw %>%
  mutate(DP = ifelse(is.na(DP), 10, DP)) %>%
  filter(QUAL >= 20)

# 3. Calculate Group Averages (Overall Mutation Burden)
cohort_avg <- all_vars %>%
  group_by(cohort) %>%
  summarise(
    total_samples           = n_distinct(sample_id),
    total_variants          = n(),
    avg_variants_per_sample = round(n() / total_samples, 1),
    .groups                 = "drop"
  )

cat("\n--- COHORT MUTATION AVERAGES ---\n")
print(cohort_avg)

# 4. Calculate Transitions vs Transversions per cohort
titv_analysis <- all_vars %>%
  mutate(
    mut_type = paste0(REF, "->", ALT),
    is_transition = case_when(
      mut_type %in% c("A->G", "G->A", "C->T", "T->C") ~ "Transition",
      mut_type %in% c("A->C", "C->A", "A->T", "T->A", "G->C", "C->G", "G->T", "T->G") ~ "Transversion",
      TRUE ~ "Indel/Other"
    )
  ) %>%
  group_by(cohort, is_transition) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = is_transition, values_from = count, values_fill = 0) %>%
  mutate(Ti_Tv_Ratio = round(Transition / Transversion, 2))

cat("\n--- TRANSITION / TRANSVERSION RATIO ---\n")
print(titv_analysis)

# 5. Bin variants by Allele Frequency (Major vs Minor)
af_spectrum <- all_vars %>%
  mutate(
    variant_class = ifelse(AF >= 0.50, "Major Variant (AF >= 0.50)", "Minor Variant (AF < 0.50)")
  ) %>%
  group_by(cohort, variant_class) %>%
  summarise(
    total_count    = n(),
    avg_per_sample = round(n() / 3, 1),
    .groups        = "drop"
  ) %>%
  complete(cohort, variant_class, fill = list(total_count = 0, avg_per_sample = 0))

cat("\n--- ALLELE FREQUENCY SPECTRUM ---\n")
print(af_spectrum)

# 6. Count variants per sample individually & Statistical Test
sample_counts <- all_vars %>%
  count(sample_id, cohort, name = "variant_count")

stat_test <- t.test(variant_count ~ cohort, data = sample_counts)

cat("\n--- STATISTICAL TEST RESULTS ---\n")
cat("P-value:", stat_test$p.value, "\n")
if (stat_test$p.value < 0.05) {
  cat("Conclusion: The difference in variant counts between Mild and Severe is STATISTICALLY SIGNIFICANT (p < 0.05).\n")
} else {
  cat("Conclusion: The difference is NOT statistically significant (p >= 0.05). Larger sample size needed.\n")
}

# 7. Generate Visualization
variant_plot <- ggplot(sample_counts, aes(x = cohort, y = variant_count, fill = cohort)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(size = 3, position = position_jitter(width = 0.1, seed = 42)) +
  labs(
    title = "Intra-Host Variant Burden: Mild vs Severe Dengue",
    x     = "Clinical Cohort", 
    y     = "Variant Count per Sample"
  ) +
  theme_minimal()

print(variant_plot)