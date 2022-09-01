# direct-conversion

direct reprogramming from fibroblast to iNeuron, using shPTBP1

## Scripts

store metadata.R

:   single cell RNAseq raw data Î•? Î∂àÎü¨?Ñú QC, filter cells & features, normalization, scaling, clustering (PCA, UMAP), UMAP clustering ?ùò Í≤∞Í≥ºÎ•? rds object Î°? Î≥Ä?ôò?ï¥?Ñú ??Ä?û•

:   input files (from google drive)

    -   1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv

    -   1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv

    output files (to google drive)

    -   20220831.rds

automatically assign cell types.R

:   scRNA-seq ?ç∞?ù¥?Ñ∞?ùò umap clustering Í≤∞Í≥ºÎ•? Î∂àÎü¨?ì§?ó¨?Ñú cluster annotation

    annotation Í≤∞Í≥º??Ä plot ?ùÑ ??Ä?û•

    cell type marker genes ?ç∞?ù¥?Ñ∞?Öã

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

:   neural clusters vs. other clusters ?ùò DE Î∂ÑÏÑùÍ≥? GO Î∂ÑÏÑù

:   input files

    -   annotation 20220826.rds

    output files

    -   glutamatergic.csv

    -   gabaergic.csv

    -   neuro.csv

    -   20220805 volcano plot.png

20220509 sctype function-1.R

:   Î≠òÌïò?äî ?ä§?Å¨Î¶ΩÌÅ¨?ù∏ÏßÄ? ?ñ¥?ñ§ ?ï®?àòÍ∞Ä ?†ï?ùò?êò?ñ¥ ?ûà?äîÏßÄ, ?ñ¥?ñ§ ?ä§?Å¨Î¶ΩÌÅ¨?óê?Ñú ?Ç¨?ö©?êò?äîÏßÄ automatically assign cell types.R

20220509 sctype function-2.R

## Manuscript

[manuscript draft](https://docs.google.com/document/d/1l4pwT3x1fijsgsGIXOH8NUQPtDORrtGNQYo2YgqYKiE/edit?usp=sharing) (google drive)

## Data

[access data](https://drive.google.com/drive/folders/11PFSiti3EtbPt2UwwIpIlMXDQNfXhRNq)

## Resources

[markdown tutorial](https://www.markdowntutorial.com/)
