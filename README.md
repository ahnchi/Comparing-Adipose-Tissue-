# Comparing-Adipose-Tissue
Repository for "Years of endurance exercise training favorably remodels abdominal subcutaneous adipose tissue in adults with overweight/obesity"

#This repository contains global proteomics and phosphoproteomics data for following article:
"Years of endurance exercise training favorably remodels abdominal subcutaneous adipose tissue in adults with overweight/obesity"

#Global proteomics
  - Any proteins that has more than one NA were filtered out.
  - Overrepresentation pathway analysis was completed through Database for Annotation, Visualization and Integrated Discovery (DAVID).
  - Rank based pathway analysis was completed through STRING databases. 

#Phosphoproteomics
  - Phosphopeptides were pre-processed for median-based normalization and imputation.
  - Imputation was done using NAguideR (http://www.omicsolution.org/wukong/NAguideR/). 'Impseq' method was selected based on classic and proteomic criteria suggested by Wang, Shisheng, et al. "NAguideR: performing and prioritizing missing value imputations for consistent bottom-up proteomic analyses." Nucleic acids research 48.14 (2020): e83-e83.
  - Kinase-enrichment analysis was done through KEA3 (https://maayanlab.cloud/kea3/), using ptmsigdb. Predicted kinase data are: 'kinase_EX.csv' and 'kinase_SED.csv'. 
