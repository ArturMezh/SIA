# 
Codes accompanying our manuscript "The ratio of cytotoxic lymphocytes to M2-like macrophages is prognostic in immunogenic tumors and predicts immunotherapy response"


The two scripts “CRC_data preparation.R” and “VALIDATION_data preparation” are used to pre-process the data for colorectal cancer cohort and for validation cohorts respectively.
These two scripts were used to read the output files from the inForm-processed Vectra Polaris system (Akoya) digitalized multispectral images, apply the thresholds to define marker positivity, using (thresholds established by GeneVia Technologies (Tampere, Finland), align the TMA position to patients (‘cases’), compute the surface area, compute the density of immune cell subsets and match the datasets to clinical meta-data.
As result, two output databases are generated: Database_working_CRC.txt and Database_working_Validation_cohort.txt. These databases are used for further analysis.



Script “CRC_survival_analysis.R” performs:
1)	Survival analysis of immune cell infiltration in colon cancer cohort stages I-III (Figure 1B)
2)	Survival analysis in colon cancer cohort stages I-III of SIA (OS and PHS, Figure 2A)
3)	Survival analysis in colon cancer cohort stages I-III of IS (OS and PHS, Figure 2B)
4)	performs iAUC computations for clinical parameters, SIA and IS and makes related illustrations (Figure 2C)
5)	performs multivariable analysis for clinical parameters, SIA and IS, for table (Table 1) and related illustrations (Figure 2D)
6)	Survival analysis in colon cancer cohort stage II of SIA (OS, Figure S1A)
7)	Survival analysis in CRC stage IV of SIA (OS, Figure S1B)


Script “VALIDATION_cohorts_survival_analysis.R” performs:
1)	Survival analysis in validation cohorts (Figure 3A)
2)	Computes Immmunoscore-like metric, and performs iAUC computations (Figure 3B)


Scripts for single cell RNA data processing.

Scripts to prepare single cell RNA data from different cohorts: 
“scRNA CRC cohort data preparation.R”
“scRNA Lung cohort data preparation.R”
“scRNA Melanoma cohort data preparation.R”
“scRNA normal tissue Atlas data preparation.R”
“scRNA Melanoma ICI cohort data preparation.R”
“scRNA RCC ICI cohort preparation.R”

For single cell RNA data analysis: 

Script “scRNA CRC analysis.R” analyses the data from scCRC cohort (E-MTAB-8410, Lee et al., 2020) and
1)	plots the tSNE maps for Figure S3
2)	plots heatmaps of genes expressed in macrophage subgroups (Figure 4A, and Figure S2)
3)	presents differentially expressed genes (Table S3)

Script “scRNA Lung analysis.R” analyses the data from Lung cancer cohort (E-MTAB-6653, Lambrechts et al., 2018) and
1)	plots the tSNE maps for Figure S4
2)	plots heatmaps of genes expressed in macrophage subgroups (Figure S2)
3)	makes violine plots of genes expressed in macrophage subgroups (Figure 4B)
4)	presents differentially expressed genes (Table S4)

Script “scRNA Melanoma analysis.R” analyses the data from uveal melanoma (Durante et al., 2020 and
1)	plots heatmaps of genes expressed in macrophage subgroups (Figure S2)
2)	makes violine plots of genes expressed in macrophage subgroups (Figure 4C)
3)	presents differentially expressed genes (Table S5)


Script “scRNA normal tissue Atlas.R” analyses the data from 15 different non-malignant organs of the same individual (GSE159929, He et al., 2020) and
1)	analyses percentages of different cell subsets (further transformed outside R environment for Table S6)
2)	makes violine plots of genes expressed in macrophage subgroups (Figure 4D)

Script “scRNA Melanoma ICI analysis.R” analyses the data from melanoma patients (GSE120575, Sade-Feldman et al., 2018) and
1)	makes bar plots illustrating SIA in different response groups and plots AUC curves (Figure 5B)

Script “scRNA RCC ICI analysis.R” analyses the data from renal cell carcinoma patients treated with ICI (Single Cell Portal: dbGaP: phs002065.v1.p1, Bi et al., 2021) and
1)	makes illustration demonstrating SIA in different response groups (Figure 5C)





