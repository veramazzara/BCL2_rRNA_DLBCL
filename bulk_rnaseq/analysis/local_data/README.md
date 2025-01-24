**Data**

This folder, local_data, includes data necessary for running some scripts of the analysis folder as well as data generated as output of these scripts.

Specifically:

- **total_counts_bcl2.txt**      : the txt file includes raw counts (see bcl2_rnaseq_deseq.R)
- **rlog_bcl2_global.rds**       : this .rds object contains the deseq object after the application of the rlog transformation in the DESeq2 workflow (see bcle_ranseq_deseq.R)
- **3d_pca_bcl2.rds**            : this rds object includes 3d pca cooridnates used to create a 3D PCA plot (see 3D_pca_bcl2.ipynb)
- **DE_analysis_filtered_*.csv** : these csv files contain the differential expressed genes for each comparison of interest (see Supplementary File 2, bcl2_rnaseq_deseq.R)


