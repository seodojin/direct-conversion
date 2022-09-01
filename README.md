# direct-conversion

direct reprogramming from fibroblast to iNeuron, using shPTBP1

## Scripts

store metadata.R

:   single cell RNAseq raw data λ₯? λΆλ¬? QC, filter cells & features, normalization, scaling, clustering (PCA, UMAP), UMAP clustering ? κ²°κ³Όλ₯? rds object λ‘? λ³??΄? ???₯

:   input files (from google drive)

    -   1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv

    -   1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv

    output files (to google drive)

    -   20220831.rds

automatically assign cell types.R

:   scRNA-seq ?°?΄?°? umap clustering κ²°κ³Όλ₯? λΆλ¬?€?¬? cluster annotation

    annotation κ²°κ³Ό?? plot ? ???₯

    cell type marker genes ?°?΄?°?

    -   cells \<- brain (tissue type) \<- ScTypeDB_full \<- PanglaoDB,  CellMarker database

    -   myofibroblast \<- all tissue \<- ScTypeDB_full \<- PanglaoDB,  CellMarker database

    -   fibroblast \<- all tissue \<- PanglaoDB

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

:   neural clusters vs. other clusters ? DE λΆμκ³? GO λΆμ

:   input files

    -   annotation 20220826.rds

    output files

    -   glutamatergic.csv

    -   gabaergic.csv

    -   neuro.csv

    -   20220805 volcano plot.png

20220509 sctype function-1.R

:   λ­ν? ?€?¬λ¦½ν¬?Έμ§? ?΄?€ ?¨?κ° ? ???΄ ??μ§, ?΄?€ ?€?¬λ¦½ν¬?? ?¬?©??μ§ automatically assign cell types.R

20220509 sctype function-2.R

## Manuscript

[manuscript draft](https://docs.google.com/document/d/1l4pwT3x1fijsgsGIXOH8NUQPtDORrtGNQYo2YgqYKiE/edit?usp=sharing) (google drive)

## Data

[access data](https://drive.google.com/drive/folders/11PFSiti3EtbPt2UwwIpIlMXDQNfXhRNq)

## Resources

[markdown tutorial](https://www.markdowntutorial.com/)
