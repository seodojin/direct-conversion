# Molecular characterization of directly reprogrammed neurons using single cell RNA sequencing

Direct reprogramming from fibroblast to iNeuron, using shPTBP1

This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq data for our paper:

## Scripts

**pre_processing.R**

-   perform QC and select cells, normalize the data, and scale the data.

-   perform (non)linear dimensional reduction and cluster the cells.

**automatically_assign_cell_types.R**

-   Load and cluster the data, assign cell types, and automatically detect the tissue type of the dataset.

**DE_analysis.R**

-   Perform DE analysis within the same cell type across conditions. - neural clusters vs. other clusters

**Trajectory_analysis.R**

-   Run principal component analysis, calculate pseudotime, and perform trajectory analysis with Slingshot.

## Data availability

Please direct any requests for data to the corresponding author.
