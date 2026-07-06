# Set your working directory
setwd("C:/Users/ssp5682/Desktop/mohanty_lab/with_Hunter_BSH_Final_subset_output_combined_tables")

library(data.table)
library(tidyverse)
library(pheatmap)

input_FASST_red <- fread("C:/Users/ssp5682/Desktop/mohanty_lab/with_Hunter_BSH_Final_subset_output_combined_tables/Bile_Amidates_heatmap_labels.csv")

fasst_batch_output <- fread("C:/Users/ssp5682/Desktop/mohanty_lab/with_Hunter_BSH_Final_subset_output_combined_tables/matches.tsv")

fasst_batch_output_filtered <- fasst_batch_output %>%
  dplyr::filter(Cosine > 0.7) |> 
  dplyr::filter(`Matching Peaks` >= 2)

fasst_batch <- left_join(fasst_batch_output_filtered, input_FASST_red, by = c("Scan" = "Feature_ID"))

molecules_interest_batch_mode <- fasst_batch |> 
  dplyr::rename(Compound = Bile_amidates) |>
  dplyr::mutate(USI_duplicate = USI, .after = USI)

processed_redu_metadata <- "all_sampleinformation.tsv"

if (!file.exists(file.path(getwd(), processed_redu_metadata))) {
  redu_url <- "https://redu.gnps2.org/dump"
  options(timeout = 600)
  download.file(redu_url, file.path(getwd(), processed_redu_metadata), mode = "wb")
  redu_metadata <- data.table::fread(processed_redu_metadata)
} else {
  redu_metadata <- data.table::fread(processed_redu_metadata)
}

MassiveID_filename <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":"))
  combined <- paste(parts[2], parts[length(parts)], sep = ":")
  return(combined)
}

molecules_interest_batch_mode$USI <- vapply(molecules_interest_batch_mode$USI, MassiveID_filename, FUN.VALUE = character(1))

ReDU_USI <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":"))
  combined <- paste(parts[2], parts[length(parts)], sep = ":")
  return(combined)
}

redu_metadata$USI <- vapply(redu_metadata$USI, ReDU_USI, FUN.VALUE = character(1))

ReDU_MASST <- left_join(molecules_interest_batch_mode, redu_metadata, by = "USI", relationship = "many-to-many")

ReDU_MASST_standardize <- ReDU_MASST |> 
  dplyr::mutate(
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, 'skin of trunk|skin of pes|head or neck skin|axilla skin|skin of manus|arm skin|skin of leg', 'skin'),
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, 'blood plasma|blood serum', 'blood'),
    HealthStatus = str_replace(HealthStatus, 'Chronic Illness', 'chronic illness'),
    HealthStatus = str_replace(HealthStatus, 'Healthy', 'healthy')
  )

df_humans <- ReDU_MASST_standardize |>  
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

list_rattus_mus <- c('10088|Mus', '10090|Mus musculus', '10105|Mus minutoides', '10114|Rattus', '10116|Rattus norvegicus')

df_rodents <- ReDU_MASST_standardize |>  
  dplyr::filter(NCBITaxonomy %in% list_rattus_mus)

# Diagnostic
analyze_counts <- function(df, column_interest) {
  df_body_parts <- df |> distinct(across(all_of(column_interest)))
  df_BodyPartName_counts <- df |> count(across(all_of(column_interest)), name = "Counts_fastMASST")
  compounds <- df |> 
    group_by(across(all_of(column_interest))) |> 
    summarise(Compounds = n_distinct(Compound), CompoundsList = toString(unique(Compound))) |> 
    ungroup()
  combined <- df_body_parts |> 
    left_join(df_BodyPartName_counts, by = column_interest) |> 
    left_join(compounds, by = column_interest)
  return(combined)
}

body_counts_humans <- analyze_counts(df_humans, "UBERONBodyPartName")
head(body_counts_humans)
body_counts_rodents <- analyze_counts(df_rodents, "UBERONBodyPartName")
head(body_counts_rodents)

# Anatomical orders
anatomical_order <- c("blood", "milk", "small intestine", "mucosa", "urine", "feces")
rodents_anatomical_order <- c("upper digestive tract", "liver", "gallbladder", "bile", "duodenum", "jejunum", "ileum", "small intestine", "caecum", "large intestine", "colon", "feces")

# Colors
colors_version <- c("#FFFFFF", "#5B8DB8", "#C07B74" )
color_gradient <- colorRampPalette(colors_version)
gradient_colors <- color_gradient(500)

# ReDU count tables
df_redu_humans <- redu_metadata |> dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")
df_redu_rodents <- redu_metadata |> dplyr::filter(NCBITaxonomy %in% list_rattus_mus)

human_ReDU_LifeStage <- df_redu_humans |> dplyr::count(LifeStage) |> dplyr::rename(LifeStage_counts = n)
human_ReDU_LifeStage$LifeStage_counts <- as.numeric(human_ReDU_LifeStage$LifeStage_counts)

human_ReDU_DOIDCommonName <- df_redu_humans |> dplyr::count(DOIDCommonName) |> dplyr::rename(DOIDCommonName_counts = n)
human_ReDU_DOIDCommonName$DOIDCommonName_counts <- as.numeric(human_ReDU_DOIDCommonName$DOIDCommonName_counts)

human_ReDU_HealthStatus <- df_redu_humans |> dplyr::count(HealthStatus) |> dplyr::rename(HealthStatus_counts = n)
human_ReDU_HealthStatus$HealthStatus_counts <- as.numeric(human_ReDU_HealthStatus$HealthStatus_counts)

human_ReDU_BiologicalSex <- df_redu_humans |> dplyr::count(BiologicalSex) |> dplyr::rename(BiologicalSex_counts = n)
human_ReDU_BiologicalSex$BiologicalSex_counts <- as.numeric(human_ReDU_BiologicalSex$BiologicalSex_counts)

human_ReDU_UBERONBodyPartName <- df_redu_humans |> dplyr::count(UBERONBodyPartName) |> dplyr::rename(UBERONBodyPartName_counts = n)
rodent_ReDU_UBERONBodyPartName <- df_redu_rodents |> dplyr::count(UBERONBodyPartName) |> dplyr::rename(UBERONBodyPartName_counts = n)
rodent_ReDU_BiologicalSex <- df_redu_rodents |> dplyr::count(BiologicalSex) |> dplyr::rename(BiologicalSex_counts = n)

# ---- DISEASE HEATMAP - FULLY NORMALIZED ----
grouped_df_humans <- df_humans |>
  group_by(Compound, DOIDCommonName) |>  
  summarise(Count = n()) |> 
  ungroup()

grouped_df_humans_pivot_table <- grouped_df_humans |>
  pivot_wider(names_from = Compound, values_from = Count, values_fill = list(Count = 0))

merged_DOID_humans <- left_join(grouped_df_humans_pivot_table, human_ReDU_DOIDCommonName, by = "DOIDCommonName")

merged_DOID_humans <- merged_DOID_humans |> 
  dplyr::filter(DOIDCommonName != "missing value") |> 
  dplyr::select(DOIDCommonName, where(~ !all(. == 0))) |> 
  dplyr::filter(rowSums(across(where(is.numeric))) != 0)

merged_DOID_humans$DOIDCommonName <- gsub("Crohn's disease", "crohn's disease", merged_DOID_humans$DOIDCommonName)

columns_to_normalize <- setdiff(names(merged_DOID_humans), c("DOIDCommonName", "DOIDCommonName_counts"))

normalized_merged_DOID_humans <- merged_DOID_humans |>
  dplyr::filter(!DOIDCommonName == "missing value") |> 
  dplyr::mutate(across(all_of(columns_to_normalize), ~ .x / .data$DOIDCommonName_counts)) |> 
  dplyr::select(-DOIDCommonName_counts)

sums <- colSums(dplyr::select(normalized_merged_DOID_humans, where(is.numeric)), na.rm = TRUE)
sums_df <- as.data.frame(t(sums))
sums_df$DOIDCommonName <- 'Sum'
sums_df <- sums_df[, names(normalized_merged_DOID_humans)]
merged_sum_humans_DOID <- bind_rows(normalized_merged_DOID_humans, sums_df)
merged_sum_humans_DOID <- merged_sum_humans_DOID |> 
  dplyr::filter(!is.na(DOIDCommonName)) |>
  dplyr::mutate(across(where(is.numeric), ~replace_na(.x, 0)))

merged_sum_humans_DOID_percentage <- merged_sum_humans_DOID |> 
  dplyr::mutate(across(all_of(columns_to_normalize), ~ .x / .x[n()] * 100))

merged_sum_humans_DOID_percentage_plot <- merged_sum_humans_DOID_percentage |>
  dplyr::filter(DOIDCommonName != "Sum") |> 
  dplyr::arrange(DOIDCommonName) |> 
  tibble::column_to_rownames("DOIDCommonName") |>
  t()

Diseases_humans <- pheatmap(merged_sum_humans_DOID_percentage_plot,
                            color = gradient_colors,
                            cluster_rows = TRUE,
                            cluster_cols = FALSE,
                            angle_col = 90,
                            main = "Health phenotype association",
                            fontsize = 10,
                            cellwidth = 15,
                            cellheight = 15,
                            treeheight_row = 100,
                            fontsize_row = 12,
                            fontsize_col = 12,
                            legend_fontsize = 10,
                            border_color = "grey90")
Diseases_humans
ggsave("Diseases_humans.pdf", plot = Diseases_humans, width = 10, height = 10, dpi = 900)

# Disease barplot - files per disease
disease_file_counts <- df_humans |>
  dplyr::filter(DOIDCommonName != "missing value") |>
  dplyr::group_by(DOIDCommonName) |>
  dplyr::summarise(FileCount = n_distinct(USI)) |>
  dplyr::ungroup() |>
  dplyr::mutate(DOIDCommonName = gsub("Crohn's disease", "crohn's disease", DOIDCommonName))

heatmap_col_order <- colnames(merged_sum_humans_DOID_percentage_plot)

Disease_barplot <- ggplot(disease_file_counts |>
                            dplyr::filter(DOIDCommonName %in% heatmap_col_order) |>
                            dplyr::mutate(DOIDCommonName = factor(DOIDCommonName, levels = heatmap_col_order)),
                          aes(x = DOIDCommonName, y = FileCount)) +
  geom_bar(stat = "identity", fill =  "#5B8DB8") +
  scale_y_continuous(breaks = seq(0, 70, by = 10),
                     limits = c(0, 70),
                     expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank(),
        plot.margin = margin(t = 50, r = 100, b = 50, l = 50, unit = "pt")) +
  labs(x = "", y = "Number of files")
Disease_barplot
ggsave("Disease_barplot.pdf", plot = Disease_barplot, width = 8.27, height = 11.69, dpi = 900)

# Compound barplot aligned to disease heatmap row order
compound_file_counts <- df_humans |>
  dplyr::filter(DOIDCommonName %in% heatmap_col_order) |>
  dplyr::group_by(Compound) |>
  dplyr::summarise(FileCount = n_distinct(USI)) |>
  dplyr::ungroup()

heatmap_row_order <- rownames(merged_sum_humans_DOID_percentage_plot)[Diseases_humans$tree_row$order]

Compound_barplot <- ggplot(compound_file_counts |>
                             dplyr::filter(Compound %in% heatmap_row_order) |>
                             dplyr::mutate(Compound = factor(Compound, levels = rev(heatmap_row_order))),
                           aes(x = FileCount, y = Compound)) +
  geom_bar(stat = "identity", fill = "#5B8DB8") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank(),
        plot.margin = margin(t = 50, r = 50, b = 50, l = 50, unit = "pt")) +
  labs(x = "Number of files", y = "")
Compound_barplot
ggsave("Compound_barplot.pdf", plot = Compound_barplot, width = 8.27, height = 11.69, dpi = 900)

# ---- FULLY NORMALIZED ORGAN HEATMAP - HUMANS ----
grouped_organ_humans <- df_humans |>
  dplyr::filter(UBERONBodyPartName != "missing value") |>
  group_by(Compound, UBERONBodyPartName) |>
  summarise(Count = n(), .groups = "drop")

grouped_organ_humans_pivot <- grouped_organ_humans |>
  pivot_wider(names_from = Compound, values_from = Count, values_fill = 0)

merged_organ_humans <- left_join(grouped_organ_humans_pivot, human_ReDU_UBERONBodyPartName,
                                 by = "UBERONBodyPartName") |>
  dplyr::filter(UBERONBodyPartName != "missing value") |>
  dplyr::select(UBERONBodyPartName, where(~ !all(. == 0)))

columns_to_normalize_organ <- setdiff(names(merged_organ_humans), c("UBERONBodyPartName", "UBERONBodyPartName_counts"))

normalized_organ_humans <- merged_organ_humans |>
  dplyr::mutate(across(all_of(columns_to_normalize_organ), ~ .x / .data$UBERONBodyPartName_counts)) |>
  dplyr::select(-UBERONBodyPartName_counts)

sums_organ_humans <- colSums(dplyr::select(normalized_organ_humans, where(is.numeric)), na.rm = TRUE)
sums_df_organ_humans <- as.data.frame(t(sums_organ_humans))
sums_df_organ_humans$UBERONBodyPartName <- 'Sum'
sums_df_organ_humans <- sums_df_organ_humans[, names(normalized_organ_humans)]
merged_sum_organ_humans <- bind_rows(normalized_organ_humans, sums_df_organ_humans)
merged_sum_organ_humans <- merged_sum_organ_humans |>
  dplyr::filter(!is.na(UBERONBodyPartName)) |>
  dplyr::mutate(across(where(is.numeric), ~replace_na(.x, 0)))

merged_sum_organ_humans_percentage <- merged_sum_organ_humans |>
  dplyr::mutate(across(all_of(columns_to_normalize_organ), ~ .x / .x[n()] * 100))

organ_humans_plot <- merged_sum_organ_humans_percentage |>
  dplyr::filter(UBERONBodyPartName != "Sum") |>
  dplyr::arrange(UBERONBodyPartName) |>
  tibble::column_to_rownames("UBERONBodyPartName") |>
  t()

organ_humans_plot <- as.data.frame(organ_humans_plot) |>
  dplyr::mutate(across(everything(), as.numeric))

organ_humans_plot <- organ_humans_plot[, anatomical_order]

Organ_humans_normalized <- pheatmap(organ_humans_plot,
                                    color = gradient_colors,
                                    cluster_rows = TRUE,
                                    cluster_cols = FALSE,
                                    angle_col = 90,
                                    main = "Organ distribution in humans (normalized)",
                                    fontsize = 10,
                                    cellwidth = 15,
                                    cellheight = 15,
                                    treeheight_row = 100,
                                    fontsize_row = 12,
                                    fontsize_col = 12,
                                    legend_fontsize = 10,
                                    border_color = "grey90")
ggsave("Organ_humans_normalized.pdf", plot = Organ_humans_normalized, width = 10, height = 10, dpi = 900)

# ---- FULLY NORMALIZED ORGAN HEATMAP - RODENTS ----
grouped_organ_rodents <- df_rodents |>
  dplyr::filter(UBERONBodyPartName != "missing value") |>
  group_by(Compound, UBERONBodyPartName) |>
  summarise(Count = n(), .groups = "drop")

grouped_organ_rodents_pivot <- grouped_organ_rodents |>
  pivot_wider(names_from = Compound, values_from = Count, values_fill = 0)

merged_organ_rodents <- left_join(grouped_organ_rodents_pivot, rodent_ReDU_UBERONBodyPartName,
                                  by = "UBERONBodyPartName") |>
  dplyr::filter(UBERONBodyPartName != "missing value") |>
  dplyr::select(UBERONBodyPartName, where(~ !all(. == 0)))

columns_to_normalize_organ_rodents <- setdiff(names(merged_organ_rodents), c("UBERONBodyPartName", "UBERONBodyPartName_counts"))

normalized_organ_rodents <- merged_organ_rodents |>
  dplyr::mutate(across(all_of(columns_to_normalize_organ_rodents), ~ .x / .data$UBERONBodyPartName_counts)) |>
  dplyr::select(-UBERONBodyPartName_counts)

sums_organ_rodents <- colSums(dplyr::select(normalized_organ_rodents, where(is.numeric)), na.rm = TRUE)
sums_df_organ_rodents <- as.data.frame(t(sums_organ_rodents))
sums_df_organ_rodents$UBERONBodyPartName <- 'Sum'
sums_df_organ_rodents <- sums_df_organ_rodents[, names(normalized_organ_rodents)]
merged_sum_organ_rodents <- bind_rows(normalized_organ_rodents, sums_df_organ_rodents)
merged_sum_organ_rodents <- merged_sum_organ_rodents |>
  dplyr::filter(!is.na(UBERONBodyPartName)) |>
  dplyr::mutate(across(where(is.numeric), ~replace_na(.x, 0)))

merged_sum_organ_rodents_percentage <- merged_sum_organ_rodents |>
  dplyr::mutate(across(all_of(columns_to_normalize_organ_rodents), ~ .x / .x[n()] * 100))

organ_rodents_plot <- merged_sum_organ_rodents_percentage |>
  dplyr::filter(UBERONBodyPartName != "Sum") |>
  dplyr::arrange(UBERONBodyPartName) |>
  tibble::column_to_rownames("UBERONBodyPartName") |>
  t()

organ_rodents_plot <- as.data.frame(organ_rodents_plot) |>
  dplyr::mutate(across(everything(), as.numeric))

organ_rodents_plot <- organ_rodents_plot[, rodents_anatomical_order]

Organ_rodents_normalized <- pheatmap(organ_rodents_plot,
                                     color = gradient_colors,
                                     cluster_rows = TRUE,
                                     cluster_cols = FALSE,
                                     angle_col = 90,
                                     main = "Organ distribution in rodents (normalized)",
                                     fontsize = 10,
                                     cellwidth = 15,
                                     cellheight = 15,
                                     treeheight_row = 100,
                                     fontsize_row = 12,
                                     fontsize_col = 12,
                                     legend_fontsize = 10,
                                     border_color = "grey90")
ggsave("Organ_rodents_normalized.pdf", plot = Organ_rodents_normalized, width = 10, height = 10, dpi = 900)

# ---- FULLY NORMALIZED BIOLOGICAL SEX HEATMAP - HUMANS ----
grouped_sex_humans <- df_humans |>
  dplyr::filter(BiologicalSex != "missing value") |>
  group_by(Compound, BiologicalSex) |>
  summarise(Count = n(), .groups = "drop")

grouped_sex_humans_pivot <- grouped_sex_humans |>
  pivot_wider(names_from = Compound, values_from = Count, values_fill = 0)

merged_sex_humans <- left_join(grouped_sex_humans_pivot, human_ReDU_BiologicalSex,
                               by = "BiologicalSex") |>
  dplyr::filter(BiologicalSex != "missing value") |>
  dplyr::select(BiologicalSex, where(~ !all(. == 0)))

columns_to_normalize_sex <- setdiff(names(merged_sex_humans), c("BiologicalSex", "BiologicalSex_counts"))

normalized_sex_humans <- merged_sex_humans |>
  dplyr::mutate(across(all_of(columns_to_normalize_sex), ~ .x / .data$BiologicalSex_counts)) |>
  dplyr::select(-BiologicalSex_counts)

sums_sex_humans <- colSums(dplyr::select(normalized_sex_humans, where(is.numeric)), na.rm = TRUE)
sums_df_sex_humans <- as.data.frame(t(sums_sex_humans))
sums_df_sex_humans$BiologicalSex <- 'Sum'
sums_df_sex_humans <- sums_df_sex_humans[, names(normalized_sex_humans)]
merged_sum_sex_humans <- bind_rows(normalized_sex_humans, sums_df_sex_humans)
merged_sum_sex_humans <- merged_sum_sex_humans |>
  dplyr::filter(!is.na(BiologicalSex)) |>
  dplyr::mutate(across(where(is.numeric), ~replace_na(.x, 0)))

merged_sum_sex_humans_percentage <- merged_sum_sex_humans |>
  dplyr::mutate(across(all_of(columns_to_normalize_sex), ~ .x / .x[n()] * 100))

sex_humans_plot <- merged_sum_sex_humans_percentage |>
  dplyr::filter(BiologicalSex != "Sum") |>
  dplyr::arrange(BiologicalSex) |>
  tibble::column_to_rownames("BiologicalSex") |>
  t()

sex_humans_plot <- as.data.frame(sex_humans_plot) |>
  dplyr::mutate(across(everything(), as.numeric))

BiologicalSex_humans_normalized <- pheatmap(sex_humans_plot,
                                            color = gradient_colors,
                                            cluster_rows = TRUE,
                                            cluster_cols = FALSE,
                                            angle_col = 90,
                                            main = "Biological sex distribution in humans (normalized)",
                                            fontsize = 10,
                                            cellwidth = 15,
                                            cellheight = 15,
                                            treeheight_row = 100,
                                            fontsize_row = 12,
                                            fontsize_col = 12,
                                            legend_fontsize = 10,
                                            border_color = "grey90")
ggsave("BiologicalSex_humans_normalized.pdf", plot = BiologicalSex_humans_normalized, width = 10, height = 10, dpi = 900)

# ---- FULLY NORMALIZED BIOLOGICAL SEX HEATMAP - RODENTS ----
grouped_sex_rodents <- df_rodents |>
  dplyr::filter(BiologicalSex != "missing value") |>
  group_by(Compound, BiologicalSex) |>
  summarise(Count = n(), .groups = "drop")

grouped_sex_rodents_pivot <- grouped_sex_rodents |>
  pivot_wider(names_from = Compound, values_from = Count, values_fill = 0)

merged_sex_rodents <- left_join(grouped_sex_rodents_pivot, rodent_ReDU_BiologicalSex,
                                by = "BiologicalSex") |>
  dplyr::filter(BiologicalSex != "missing value") |>
  dplyr::select(BiologicalSex, where(~ !all(. == 0)))

columns_to_normalize_sex_rodents <- setdiff(names(merged_sex_rodents), c("BiologicalSex", "BiologicalSex_counts"))

normalized_sex_rodents <- merged_sex_rodents |>
  dplyr::mutate(across(all_of(columns_to_normalize_sex_rodents), ~ .x / .data$BiologicalSex_counts)) |>
  dplyr::select(-BiologicalSex_counts)

sums_sex_rodents <- colSums(dplyr::select(normalized_sex_rodents, where(is.numeric)), na.rm = TRUE)
sums_df_sex_rodents <- as.data.frame(t(sums_sex_rodents))
sums_df_sex_rodents$BiologicalSex <- 'Sum'
sums_df_sex_rodents <- sums_df_sex_rodents[, names(normalized_sex_rodents)]
merged_sum_sex_rodents <- bind_rows(normalized_sex_rodents, sums_df_sex_rodents)
merged_sum_sex_rodents <- merged_sum_sex_rodents |>
  dplyr::filter(!is.na(BiologicalSex)) |>
  dplyr::mutate(across(where(is.numeric), ~replace_na(.x, 0)))

merged_sum_sex_rodents_percentage <- merged_sum_sex_rodents |>
  dplyr::mutate(across(all_of(columns_to_normalize_sex_rodents), ~ .x / .x[n()] * 100))

sex_rodents_plot <- merged_sum_sex_rodents_percentage |>
  dplyr::filter(BiologicalSex != "Sum") |>
  dplyr::arrange(BiologicalSex) |>
  tibble::column_to_rownames("BiologicalSex") |>
  t()

sex_rodents_plot <- as.data.frame(sex_rodents_plot) |>
  dplyr::mutate(across(everything(), as.numeric))

BiologicalSex_rodents_normalized <- pheatmap(sex_rodents_plot,
                                             color = gradient_colors,
                                             cluster_rows = TRUE,
                                             cluster_cols = FALSE,
                                             angle_col = 90,
                                             main = "Biological sex distribution in rodents (normalized)",
                                             fontsize = 10,
                                             cellwidth = 15,
                                             cellheight = 15,
                                             treeheight_row = 100,
                                             fontsize_row = 12,
                                             fontsize_col = 12,
                                             legend_fontsize = 10,
                                             border_color = "grey90")
ggsave("BiologicalSex_rodents_normalized.pdf", plot = BiologicalSex_rodents_normalized, width = 10, height = 10, dpi = 900)

# ---- FULLY NORMALIZED LIFESTAGE HEATMAP - HUMANS ----
lifestage_order <- c("Infancy (<2 yrs)", 
                     "Early Childhood (2 yrs < x <=8 yrs)",
                     "Adolescence (8 yrs < x <= 18 yrs)",
                     "Early Adulthood (18 yrs < x <= 45 yrs)",
                     "Middle Adulthood (45 yrs < x <= 65 yrs)",
                     "Later Adulthood (>65 yrs)")

grouped_lifestage_humans <- df_humans |>
  group_by(Compound, LifeStage) |>
  summarise(Count = n()) |>
  ungroup()

grouped_lifestage_humans_pivot <- grouped_lifestage_humans |>
  pivot_wider(names_from = Compound, values_from = Count, values_fill = list(Count = 0))

merged_LifeStage_humans <- left_join(grouped_lifestage_humans_pivot, human_ReDU_LifeStage, by = "LifeStage")

merged_LifeStage_humans <- merged_LifeStage_humans |>
  dplyr::filter(LifeStage != "missing value") |>
  dplyr::select(LifeStage, where(~ !all(. == 0))) |>
  dplyr::filter(rowSums(across(where(is.numeric))) != 0)

columns_to_normalize_LS <- setdiff(names(merged_LifeStage_humans), c("LifeStage", "LifeStage_counts"))

normalized_merged_LifeStage_humans <- merged_LifeStage_humans |>
  dplyr::mutate(across(all_of(columns_to_normalize_LS), ~ .x / .data$LifeStage_counts)) |>
  dplyr::select(-LifeStage_counts)

sums_LS <- colSums(dplyr::select(normalized_merged_LifeStage_humans, where(is.numeric)), na.rm = TRUE)
sums_df_LS <- as.data.frame(t(sums_LS))
sums_df_LS$LifeStage <- 'Sum'
sums_df_LS <- sums_df_LS[, names(normalized_merged_LifeStage_humans)]
merged_sum_humans_LS <- bind_rows(normalized_merged_LifeStage_humans, sums_df_LS)
merged_sum_humans_LS <- merged_sum_humans_LS |>
  dplyr::filter(!is.na(LifeStage)) |>
  dplyr::mutate(across(where(is.numeric), ~replace_na(.x, 0)))

merged_sum_humans_LS_percentage <- merged_sum_humans_LS |>
  dplyr::mutate(across(all_of(columns_to_normalize_LS), ~ .x / .x[n()] * 100))

LifeStage_humans_plot <- merged_sum_humans_LS_percentage |>
  dplyr::filter(LifeStage != "Sum") |>
  dplyr::mutate(LifeStage = factor(LifeStage, levels = lifestage_order)) |>
  dplyr::arrange(LifeStage) |>
  tibble::column_to_rownames("LifeStage") |>
  t()

LifeStage_humans_normalized <- pheatmap(LifeStage_humans_plot,
                                        color = gradient_colors,
                                        cluster_rows = TRUE,
                                        cluster_cols = FALSE,
                                        angle_col = 90,
                                        main = "Life stage distribution in humans (normalized)",
                                        fontsize = 10,
                                        cellwidth = 15,
                                        cellheight = 15,
                                        treeheight_row = 100,
                                        fontsize_row = 12,
                                        fontsize_col = 12,
                                        legend_fontsize = 10,
                                        border_color = "grey90")
ggsave("LifeStage_humans_normalized.pdf", plot = LifeStage_humans_normalized, width = 10, height = 10, dpi = 900)

# Compound barplot aligned to lifestage heatmap row order
compound_file_counts_lifestage <- df_humans |>
  dplyr::filter(LifeStage %in% lifestage_order) |>
  dplyr::group_by(Compound) |>
  dplyr::summarise(FileCount = n_distinct(USI)) |>
  dplyr::ungroup()

heatmap_row_order_lifestage <- rownames(LifeStage_humans_plot)[LifeStage_humans_normalized$tree_row$order]

Compound_barplot_lifestage <- ggplot(compound_file_counts_lifestage |>
                                       dplyr::filter(Compound %in% heatmap_row_order_lifestage) |>
                                       dplyr::mutate(Compound = factor(Compound, levels = rev(heatmap_row_order_lifestage))),
                                     aes(x = FileCount, y = Compound)) +
  geom_bar(stat = "identity", fill = "#5B8DB8") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid = element_blank(),
        plot.margin = margin(t = 50, r = 50, b = 50, l = 50, unit = "pt")) +
  labs(x = "Number of files", y = "")
Compound_barplot_lifestage
ggsave("Compound_barplot_lifestage.pdf", plot = Compound_barplot_lifestage, width = 8.27, height = 11.69, dpi = 900)
