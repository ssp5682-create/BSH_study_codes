Figure_4: Repository Scale Phenotypic Associations
Script for generating the heatmaps in Figure 4 is provided (`bile_amidates_analysis_github.R`) along with the required input tables. The input tables are the compound/scan lookup table (`Bile_Amidates_heatmap_labels.csv`) and the fastMASST search output (`matches.tsv`), obtained by searching the query MS/MS spectra against public data using fastMASST/FASST on GNPS2. Sample metadata is retrieved directly from the ReDU repository (https://redu.gnps2.org/dump) at runtime.
The following heatmaps are generated, each with associated barplots:
Organ distribution (humans and rodents)
Health phenotype association (humans)
Biological sex distribution (humans and rodents)
Life stage distribution (humans)
Final figure panels were assembled in Adobe Illustrator: font sizes were adjusted manually after export, and exact file-count labels on the barplots were added by hand using the values in the `disease_file_counts`, `compound_file_counts`, and `compound_file_counts_lifestage` tables produced by the script.
