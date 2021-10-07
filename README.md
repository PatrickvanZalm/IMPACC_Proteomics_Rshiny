# IMPACC_Proteomics_Rshiny

Rshiny code for looking at the IMPACC data (https://www.impaccstudy.org/).

App includes the following key components:

## Load data

Will load in the clinical data, DDA, DIA and targeted proteomics data. User can select some filters, normalisation and cut-off filter.

User will be given some simple statistics such as total number of proteins, and comparison of cut-off filters between the three acquisition methods.

## Clinical data

Will give user info on the clinical data in specific. Includes standard tabular patients inclusion, clinical grouping, age distribution, sex, enrollment sites and ethnicity.

## Clustering

User can do some dimensionality reduction including PCA, t-SNE and UMAP. User can also generate a PVCA plot to assess effects of clinical parameters on distribution in a PCA.

## Longitudinal

Will run ANOVA, within-subject ANOVA or non-parametric ANOVA between the visit numbers. Will give some info on sample size, as well as plotting the trajectories of proteins across the visits. User can run enrichment analysis too.

## Volcanoplotter

This project includes many clinical variables. We developed a simple script that will run a t-test between any clinical variable of interest and plot it in a volcanoplot. Again, user can run some enrichment analysis on significant proteins.

# Heatmap

User can generate a heatmap of specified proteins, as well as specified clinical factors for annotation.

