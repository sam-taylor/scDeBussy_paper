# CellAlignDTW_paper
Code to reproduce the figures in the paper


## B cell development



## COVID19 



## Richter transformation

The source data of CLL to Richter transformation came from [Nadeu et al 2022 Nat Med](https://www.nature.com/articles/s41591-022-01927-8) and the seurat objects were downloaded from [Zenodo](https://zenodo.org/records/6631966) into the `data/` folder.

```{bash}
cd data/
wget https://zenodo.org/record/6631966/files/Nadeu2022_scRNAseq.zip
unzip Nadeu2022_scRNAseq.zip
cd Nadeu2022_scRNAseq
```

The scripts that are used to extract and preprocess data are located in folder `richter_transformation`, as described in the following table. 

| Script           | Description                                               |
|------------------|-----------------------------------------------------------|
| 1_extract_data.R | Extract count matrix and annotations from the source data |
| 2_preprocess.py  | Preprocess the expression matrix into h5ad file           |
| 3_palantir.py    | Run palantir pseudotime on h5ad files                     |


## Downstream analysis


