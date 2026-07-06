# ================================================================
# Language: R / tidyverse
#
# Required packages:
# install.packages(c("dplyr", "tidyr", "stringr", "readr"))
# ================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)

dir.create("results", showWarnings = FALSE)

# ------------------------------------------------
# Load practice datasets
# ------------------------------------------------

metadata <- read_csv("data/sample_metadata.csv", show_col_types = FALSE)
qc <- read_csv("data/rna_qc.csv", show_col_types = FALSE)
counts <- read_csv("data/gene_counts.csv", show_col_types = FALSE)
sample_ids <- read_csv("data/sample_ids.csv", show_col_types = FALSE)
fastq_files <- read_csv("data/fastq_files.csv", show_col_types = FALSE)
variants <- read_csv("data/variants.csv", show_col_types = FALSE)
bed <- read_csv("data/genomic_intervals.csv", show_col_types = FALSE)
cells <- read_csv("data/single_cell_qc.csv", show_col_types = FALSE)
single_cell_counts <- read_csv("data/single_cell_counts_long.csv", show_col_types = FALSE)

# ================================================================
# 1. Sample metadata cleaning
# Goal:
# - Standardize condition values
# - Lowercase tissue values
# - Remove rows missing patient_id
# - Identify duplicate patients
# ================================================================

metadata_clean <- metadata %>%
  mutate(
    condition = str_to_lower(condition),
    tissue = str_to_lower(tissue)
  ) %>%
  filter(!is.na(patient_id), patient_id != "")

duplicate_patients <- metadata_clean %>%
  group_by(patient_id) %>%
  filter(n() > 1) %>%
  ungroup()

write_csv(metadata_clean, "results/01_metadata_clean.csv")
write_csv(duplicate_patients, "results/01_duplicate_patients.csv")

# ================================================================
# 2. Merge metadata with sequencing QC metrics
# Goal:
# - Full join metadata and QC
# - Identify samples present in one file but not the other
# ================================================================

merged_qc <- metadata %>%
  mutate(in_metadata = TRUE) %>%
  full_join(
    qc %>% mutate(in_qc = TRUE),
    by = "sample_id"
  ) %>%
  mutate(
    merge_status = case_when(
      in_metadata == TRUE & in_qc == TRUE ~ "both",
      in_metadata == TRUE & is.na(in_qc) ~ "metadata_only",
      is.na(in_metadata) & in_qc == TRUE ~ "qc_only",
      TRUE ~ "unknown"
    )
  )

merge_summary <- merged_qc %>%
  count(merge_status)

write_csv(merged_qc, "results/02_metadata_qc_full_join.csv")
write_csv(merge_summary, "results/02_merge_summary.csv")

# ================================================================
# 3. Calculate mapping percentage
# Goal:
# - Add mapping_rate
# - Flag samples with mapping rate below 80%
# ================================================================

qc_with_mapping <- qc %>%
  mutate(
    mapping_rate = mapped_reads / total_reads * 100,
    low_mapping = mapping_rate < 80
  )

write_csv(qc_with_mapping, "results/03_qc_with_mapping_rate.csv")

# ================================================================
# 4. Reshape gene expression data from wide to long
# Goal:
# - Convert gene-by-sample count matrix into long format
# ================================================================

counts_long <- counts %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression"
  )

write_csv(counts_long, "results/04_counts_long.csv")

# ================================================================
# 5. Summarize expression by condition
# Goal:
# - Join long expression table to metadata
# - Calculate mean expression per gene by condition
# ================================================================

expression_summary <- counts_long %>%
  left_join(
    metadata_clean %>% select(sample_id, condition),
    by = "sample_id"
  ) %>%
  filter(!is.na(condition), condition != "") %>%
  group_by(gene, condition) %>%
  summarize(
    mean_expression = mean(expression, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(expression_summary, "results/05_expression_summary_by_condition.csv")

# ================================================================
# 6. Find genes with largest fold change
# Goal:
# - Convert condition means to wide format
# - Calculate log2 fold change for disease vs control
# ================================================================

expression_wide <- expression_summary %>%
  pivot_wider(
    names_from = condition,
    values_from = mean_expression,
    values_fill = 0
  )

fold_change <- expression_wide %>%
  mutate(
    log2_fold_change = log2((disease + 1) / (control + 1)),
    abs_log2_fold_change = abs(log2_fold_change)
  ) %>%
  arrange(desc(abs_log2_fold_change))

top_fold_change <- fold_change %>%
  slice_head(n = 10)

write_csv(fold_change, "results/06_fold_change_all_genes.csv")
write_csv(top_fold_change, "results/06_top_10_fold_change_genes.csv")

# ================================================================
# 7. Filter genes by expression threshold
# Goal:
# - Keep genes expressed in at least 3 samples
# ================================================================

counts_filtered <- counts %>%
  filter(
    rowSums(across(-gene, ~ .x > 0)) >= 3
  )

write_csv(counts_filtered, "results/07_counts_filtered_min_3_samples.csv")

# ================================================================
# 8. Parse sample IDs
# Example IDs:
# P001_Tumor_RNA_rep1
# Goal:
# - Extract patient_id, tissue, assay, replicate
# ================================================================

sample_ids_parsed <- sample_ids %>%
  separate(
    sample_id,
    into = c("patient_id", "tissue", "assay", "replicate"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    tissue = str_to_lower(tissue),
    replicate = str_remove(replicate, "rep"),
    replicate = as.integer(replicate)
  )

write_csv(sample_ids_parsed, "results/08_sample_ids_parsed.csv")

# ================================================================
# 9. Work with FASTQ filenames
# Example filenames:
# S1_L001_R1_001.fastq.gz
# Goal:
# - Extract sample_id, lane, read
# - Identify missing R1/R2 pairs
# ================================================================

fastq_parsed <- fastq_files %>%
  mutate(
    base = str_remove(filename, "\\.fastq\\.gz$")
  ) %>%
  separate(
    base,
    into = c("sample_id", "lane", "read", "chunk"),
    sep = "_"
  ) %>%
  select(sample_id, lane, read, filename)

missing_pairs <- fastq_parsed %>%
  group_by(sample_id, lane) %>%
  summarize(
    n_reads = n_distinct(read),
    reads_present = paste(sort(unique(read)), collapse = ","),
    .groups = "drop"
  ) %>%
  filter(n_reads < 2)

write_csv(fastq_parsed, "results/09_fastq_parsed.csv")
write_csv(missing_pairs, "results/09_missing_fastq_pairs.csv")

# ================================================================
# 10. Detect duplicate or inconsistent sample IDs
# Goal:
# - Find patient IDs with more than one condition
# ================================================================

inconsistent_patients <- metadata_clean %>%
  group_by(patient_id) %>%
  summarize(
    n_conditions = n_distinct(condition),
    conditions = paste(unique(condition), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_conditions > 1)

problem_patient_rows <- metadata_clean %>%
  semi_join(inconsistent_patients, by = "patient_id")

write_csv(inconsistent_patients, "results/10_inconsistent_patients.csv")
write_csv(problem_patient_rows, "results/10_problem_patient_rows.csv")

# ================================================================
# 11. SQL-style wrangling in R
# Goal:
# - Inner join metadata to QC and calculate mapping rate
# ================================================================

sample_qc <- metadata_clean %>%
  inner_join(qc, by = "sample_id") %>%
  mutate(
    mapping_rate = 100 * mapped_reads / total_reads
  )

write_csv(sample_qc, "results/11_sample_qc_inner_join.csv")

# Equivalent SQL:
# SELECT
#     s.sample_id,
#     s.patient_id,
#     s.condition,
#     q.total_reads,
#     q.mapped_reads,
#     100.0 * q.mapped_reads / q.total_reads AS mapping_rate
# FROM samples s
# JOIN qc_metrics q
#     ON s.sample_id = q.sample_id;

# ================================================================
# 12. Variant annotation wrangling
# Goal:
# - Filter for high-impact variants in cancer genes
# ================================================================

cancer_genes <- c("TP53", "EGFR", "KRAS", "BRAF", "PIK3CA")

filtered_variants <- variants %>%
  filter(
    impact == "HIGH",
    gene %in% cancer_genes
  )

variant_counts <- filtered_variants %>%
  count(gene, sort = TRUE)

write_csv(filtered_variants, "results/12_filtered_high_impact_variants.csv")
write_csv(variant_counts, "results/12_high_impact_variant_counts.csv")

# ================================================================
# 13. Genomic interval wrangling
# Goal:
# - Calculate interval length
# - Filter intervals longer than 100 bp
# - Summarize total covered bases by chromosome
# ================================================================

bed_with_length <- bed %>%
  mutate(
    length = end - start
  )

bed_long <- bed_with_length %>%
  filter(length > 100)

chrom_summary <- bed_with_length %>%
  group_by(chrom) %>%
  summarize(
    total_bp = sum(length),
    n_intervals = n(),
    .groups = "drop"
  )

write_csv(bed_with_length, "results/13_bed_with_length.csv")
write_csv(bed_long, "results/13_bed_intervals_longer_than_100bp.csv")
write_csv(chrom_summary, "results/13_chromosome_interval_summary.csv")

# ================================================================
# 14. Single-cell QC wrangling
# Goal:
# - Keep cells with:
#   n_genes >= 200
#   pct_counts_mt < 10
#   n_genes <= 5000
# - Summarize QC by sample
# ================================================================

filtered_cells <- cells %>%
  filter(
    n_genes >= 200,
    pct_counts_mt < 10,
    n_genes <= 5000
  )

cell_counts_after_filtering <- filtered_cells %>%
  count(sample_id, name = "n_cells_after_filtering")

single_cell_qc_summary <- cells %>%
  group_by(sample_id) %>%
  summarize(
    n_cells = n(),
    median_genes = median(n_genes, na.rm = TRUE),
    median_counts = median(total_counts, na.rm = TRUE),
    median_pct_mt = median(pct_counts_mt, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(filtered_cells, "results/14_filtered_single_cells.csv")
write_csv(cell_counts_after_filtering, "results/14_cell_counts_after_filtering.csv")
write_csv(single_cell_qc_summary, "results/14_single_cell_qc_summary.csv")

# ================================================================
# 15. Pseudobulk aggregation
# Goal:
# - Aggregate single-cell counts by sample_id, cell_type, and gene
# - Convert pseudobulk data to wide format
# ================================================================

pseudobulk <- single_cell_counts %>%
  group_by(sample_id, cell_type, gene) %>%
  summarize(
    count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

pseudobulk_wide <- pseudobulk %>%
  unite("sample_celltype", sample_id, cell_type, sep = "_") %>%
  pivot_wider(
    names_from = sample_celltype,
    values_from = count,
    values_fill = 0
  )

write_csv(pseudobulk, "results/15_pseudobulk_long.csv")
write_csv(pseudobulk_wide, "results/15_pseudobulk_wide.csv")

# ================================================================
# 16. Full mini mock interview pipeline
# Goal:
# Given sample metadata, RNA QC, and gene counts:
# 1. Standardize sample IDs
# 2. Remove samples missing condition labels
# 3. Merge metadata with QC metrics
# 4. Calculate mapping rate
# 5. Remove samples with mapping rate below 80%
# 6. Filter genes expressed in fewer than 3 samples
# 7. Convert counts wide to long
# 8. Calculate mean expression by condition
# 9. Return top 10 genes by absolute log2 fold change
# ================================================================

metadata_pipeline <- metadata %>%
  mutate(
    sample_id = str_trim(sample_id),
    sample_id = str_to_upper(sample_id),
    condition = str_to_lower(condition)
  ) %>%
  filter(
    !is.na(condition),
    condition != "",
    !is.na(patient_id),
    patient_id != ""
  )

qc_pipeline <- qc %>%
  mutate(
    sample_id = str_trim(sample_id),
    sample_id = str_to_upper(sample_id)
  )

counts_pipeline <- counts %>%
  rename_with(
    ~ str_to_upper(str_trim(.x)),
    .cols = -gene
  )

sample_qc_pipeline <- metadata_pipeline %>%
  inner_join(qc_pipeline, by = "sample_id") %>%
  mutate(
    mapping_rate = 100 * mapped_reads / total_reads
  )

passing_samples <- sample_qc_pipeline %>%
  filter(mapping_rate >= 80) %>%
  pull(sample_id)

counts_pipeline_filtered <- counts_pipeline %>%
  select(gene, all_of(passing_samples)) %>%
  filter(
    rowSums(across(-gene, ~ .x > 0)) >= 3
  )

counts_pipeline_long <- counts_pipeline_filtered %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "count"
  )

expression_by_condition_pipeline <- counts_pipeline_long %>%
  left_join(
    sample_qc_pipeline %>% select(sample_id, condition),
    by = "sample_id"
  ) %>%
  group_by(gene, condition) %>%
  summarize(
    mean_expression = mean(count, na.rm = TRUE),
    .groups = "drop"
  )

top_genes_pipeline <- expression_by_condition_pipeline %>%
  pivot_wider(
    names_from = condition,
    values_from = mean_expression,
    values_fill = 0
  ) %>%
  mutate(
    log2_fold_change = log2((disease + 1) / (control + 1)),
    abs_log2_fold_change = abs(log2_fold_change)
  ) %>%
  arrange(desc(abs_log2_fold_change)) %>%
  slice_head(n = 10)

write_csv(sample_qc_pipeline, "results/16_pipeline_sample_qc.csv")
write_csv(counts_pipeline_filtered, "results/16_pipeline_filtered_counts.csv")
write_csv(top_genes_pipeline, "results/16_pipeline_top_10_genes.csv")

# ================================================================
# Final console messages
# ================================================================

cat("\nBioinformatics wrangling practice script completed.\n")
cat("Input data are in the data/ folder.\n")
cat("Output files were written to the results/ folder.\n\n")
cat("Key outputs to inspect:\n")
cat("- results/02_merge_summary.csv\n")
cat("- results/03_qc_with_mapping_rate.csv\n")
cat("- results/06_top_10_fold_change_genes.csv\n")
cat("- results/09_missing_fastq_pairs.csv\n")
cat("- results/14_filtered_single_cells.csv\n")
cat("- results/15_pseudobulk_wide.csv\n")
cat("- results/16_pipeline_top_10_genes.csv\n")
