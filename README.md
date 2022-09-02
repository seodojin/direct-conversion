# direct-conversion

direct reprogramming from fibroblast to iNeuron, using shPTBP1

## Scripts

store metadata.R

:   single cell RNAseq raw data ??? 불러??? QC, filter cells & features, normalization, scaling, clustering (PCA, UMAP), UMAP clustering ??? 결과??? rds object ??? 변????????? ??????

:   input files (from google drive)

    -   1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv

    -   1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv

    output files (to google drive)

    -   20220831.rds

automatically assign cell types.R

:   scRNA-seq ???????????? umap clustering 결과??? 불러????????? cluster annotation

    annotation 결과??? plot ??? ??????

    cell type marker genes ????????????

    -   cells \<- brain (tissue type) \<- ScTypeDB_full \<- PanglaoDB, CellMarker database

    -   myofibroblast \<- all tissue \<- ScTypeDB_full \<- PanglaoDB, CellMarker database

    -   fibroblast \<- connective tissue \<- PanglaoDB

:   input files (google drive)

    -   20220831.rds

    -   ScTypeDB_full.xlsx

        <div>

        source: <https://github.com/IanevskiAleksandr/sc-type/>

        </div>

:   output files

    -   annotation 20220831.rds

    -   20220826 legend add umap split version.png

DE analysis.R

:   neural clusters vs. other clusters ??? DE 분석??? GO 분석

:   input files

    -   annotation 20220826.rds

    output files

    -   glutamatergic.csv

    -   gabaergic.csv

    -   neuro.csv

    -   20220805 volcano plot.png

20220509 sctype function-1.R 20220509 sctype function-2.R

:   use ScType to automatically assign cell types : load two additional ScType functions.

## Manuscript

[manuscript draft](https://docs.google.com/document/d/1l4pwT3x1fijsgsGIXOH8NUQPtDORrtGNQYo2YgqYKiE/edit?usp=sharing) (google drive)

## Data

[access data](https://drive.google.com/drive/folders/11PFSiti3EtbPt2UwwIpIlMXDQNfXhRNq)
