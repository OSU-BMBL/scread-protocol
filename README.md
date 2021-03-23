# scREAD protocol

This is the repository for the *STAR protocols* manuscript: **Use of scREAD to Explore and Analyze Single-cell and Single-nucleus RNA-Seq data for Alzheimer’s Disease**. 

The protocol is baed on **scREAD** (single-cell RNA-Seq database for Alzheimer's Disease). It is a first-of-its-kind database to provide comprehensive analysis results of all the existing single-cell RNA-Seq and single-nucleus RNA-Seq data of Alzheimer's Disease in the public domain. The database is freely available at: [https://bmbls.bmi.osumc.edu/scread/](https://bmbls.bmi.osumc.edu/scread/)

The original scREAD paper was published in *iScience*: [scREAD: A Single-Cell RNA-Seq Database for Alzheimer's Disease](https://www.sciencedirect.com/science/article/pii/S2589004220309664)

If you have any questions or feedback regarding this notebook, please contact Cankun Wang <cankun.wang@osumc.edu>.

## How to use the protocols?

- Calculating overlapping DEGs from the same cell type across datasets (Optional section 6 in the manuscript)
  - Run /overlapping_genes/overlap.rmd locally
  - Use Google Colab version: [https://colab.research.google.com/drive/1lInXa6jD4yc7RGJc0EWDfy5NNoXT1qye?usp=sharing](https://colab.research.google.com/drive/1lInXa6jD4yc7RGJc0EWDfy5NNoXT1qye?usp=sharing) 
- Optional section 7: Running scREAD backend analysis workflow locally (Optional section 7 in the manuscript)
  - open [workflow readme](https://github.com/OSU-BMBL/scread-protocol/tree/master/workflow)

## Directory structure

- overlapping_genes
  - scread_db.rdata (91MB): scREAD dataset information, differential gene expression analysis results. 
  - overlap_functions.R: functions to obtain overlapping genes.
  - overlap.rmd: R markdown version to calculate overlapping genes.
- workflow
  - custom_marker.csv: A manually created marker gene list file used for identified cell types.
  - functions.R: Visualization functions used in R.
  - build_control_atlas.R: build control cells atlas Seurat object from count matrix file.
  - transfer_cell_type.R: filter out control-like cells in disease dataset
  - run_analysis.R: run analysis workflow, and export tables in scREAD database format.
  - example_control.csv. The example control dataset.
  - example_disease.csv. The example disease dataset.

## Authors

- [Cankun Wang](https://github.com/Wang-Cankun)
- [Yujia Xiang](https://github.com/Candlelight-XYJ)

