# Molecular characterization of directly reprogrammed neurons using single cell RNA sequencing

direct reprogramming from fibroblast to iNeuron, using shPTBP1

This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq data for our paper:

## Scripts

**pre_processing.R**

:   QC and selecting cells

:   Normalizing the data

:   Scaling the data

:   Perform (non) linear dimensional reduction

:   Cluster the cells

**automatically_assign_cell_types.R**

:   Load and cluster the data

:   Cell type assignment

:   Automatically detect a tissue type of the dataset

**DE_analysis.R**

:   Perform DE analysis within the same cell type across conditions - neural clusters vs. other clusters

**Trajectory_analysis.R**

:   Run principal component analysis

:   Pseudotime Calculation

:   Trajectory Analysis with Slingshot

## Data

[access data](https://drive.google.com/drive/folders/11PFSiti3EtbPt2UwwIpIlMXDQNfXhRNq)
