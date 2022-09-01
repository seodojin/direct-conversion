# direct-conversion

direct reprogramming from fibroblast to iNeuron, using shPTBP1

## Scripts

store metadata.R

:   single cell RNAseq raw data 를 불러서 QC, filter cells & features, normalization, scaling, clustering (PCA, UMAP), UMAP clustering 의 결과를 rds object 로 변환해서 저장

:   input files (from google drive)

    -   1129_neuron_SampleTag01_hs_empty_RSEC_MolsPerCell.csv

    -   1129_neuron_SampleTag02_hs_sh_RSEC_MolsPerCell.csv

    output files (to google drive)

    -   20220831.rds

automatically assign cell types.R

:   scRNA-seq 데이터의 umap clustering 결과를 불러들여서 cluster annotation

    annotation 결과와 plot 을 저장

    cell type marker genes 데이터셋

    -   cells \<- brain (tissue type) \<- ScTypeDB_full \<- PanglaoDB, HPA

    -   myofibroblast \<- all tissue \<- ScTypeDB_full \<- PanglaoDB, HPA

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

:   neural clusters vs. other clusters 의 DE 분석과 GO 분석

:   input files

    -   annotation 20220826.rds

    output files

    -   glutamatergic.csv

    -   gabaergic.csv

    -   neuro.csv

    -   20220805 volcano plot.png

20220509 sctype function-1.R

:   뭘하는 스크립크인지? 어떤 함수가 정의되어 있는지, 어떤 스크립크에서 사용되는지 automatically assign cell types.R

20220509 sctype function-2.R

## Manuscript

[manuscript draft](https://docs.google.com/document/d/1l4pwT3x1fijsgsGIXOH8NUQPtDORrtGNQYo2YgqYKiE/edit?usp=sharing) (google drive)

## Data

[access data](https://drive.google.com/drive/folders/11PFSiti3EtbPt2UwwIpIlMXDQNfXhRNq)

## Resources

[markdown tutorial](https://www.markdowntutorial.com/)
