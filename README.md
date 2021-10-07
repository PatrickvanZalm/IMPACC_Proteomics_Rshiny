# IMPACC_Proteomics_Rshiny
Rshiny code for looking at the IMPACC data (https://www.impaccstudy.org/)

This app does the following:
1. Load in clinical data
2. Load in proteomics data (DDA, DIA, targeted)

Once data is loaded in user will first have a quick overview on some simple statistics related to the three different data files.

Next, user can have a look at clinical data including grouping for disease, disease progression, age, sex, enrollment site and ethnicity.

User can use different clustering methods, including PCA, t-SNE, UMAP and PVCA.

User can do some initial longitudinal analysis, using ANOVA or within subject ANOVA. More in depth results can be seen for specifiek proteins.

To assess differences between specific groups, user can specify two groups and run t-tests, following by a volcanoplot. Results can be enriched using STRING-db.

Heatmap can be generated. It will automatically take the significant proteins from the longitdinal analysis, but user can specify whatever proteins he/she might be interested in. 
