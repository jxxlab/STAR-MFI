# rabbit_placenta
Spatial transcriptomics of rabbit placental development reveals hemochorial placenta characteristics via cross-species comparisons
## Overview

This repository contains analysis scripts for single-cell, single-nucleus, and spatial transcriptomics of rabbit placental development, including cross-species comparisons.

# System requirements

## Software dependencies
The analyses were performed using the following software environments:

## Operating system: 
-Operating system: CentOS Linux 7 (Core)
-Linux kernel: 3.10.0-1160.108.1.e17.x86_64
-R version: 4.1.3
-Python version: 3.10.6

## Dependencies
- Seurat v4.3.0
- Scrublet v0.2.1
- Harmony v0.1.1
- ClusterProfiler package v4.2.2
- org.Hs.eg.db v3.14.0
- Monocle3 v1.3.1
- SeuratWrappers v0.3.1
- ClusterGVis v0.1.2
- SCENIC, v0.12.1
- homologene package v1.4.68.19.3.27
- stereopy v1.5.0

# Installation guide
  see requirements.txt
## Install from BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::version()    
install.packages("remotes")
options(timeout = 60000000000)
x=c("Seurat","devtools",'R.utils',"pagoda2","colorRamps","tidyverse",'IRkernel','scattermore','argparser',    "msigdbr",'harmony', 'glmpca','qs2',"jackstraw",'arrow',
    'flexmix', "factoextra", arulesViz","Hmisc","BDgraph","qgraph","psych",'Rfast2',"vembedr","magick",'irr')
x[!x %in% rownames(installed.packages())]
install.packages(x[!x %in% rownames(installed.packages())])

## Install from devtools
options(timeout = 60000000000)
options(download.file.method = "wget")
options(download.file.method = "curl")
devtools::install_github(c('satijalab/seurat-wrappers','satijalab/seurat-data','mojaveazure/seurat-disk',"mojaveazure/loomR",
                           'cole-trapnell-lab/monocle3',"aertslab/SCopeLoomR","mkearney/rmd2jupyter",  #"velocyto-team/velocyto.R",
                           "cailab-tamu/scTenifoldKnk","PaulingLiu/scibet","satijalab/azimuth","sqjin/CellChat","xzhoulab/SPARK",
                           "dmcable/spacexr",))
                           
## Install from install.packages
install.packages("ipkg")
ipkg::install_github("bnprks/BPCells/r")
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
x=c("BPCells", "presto", "glmGamPoi","ClusterGVis")
x[!x %in% rownames(installed.packages())]
install.packages(x[!x %in% rownames(installed.packages())])

install.packages("TTR")
install.packages("https://cran.r-project.org/src/contrib/Archive/smoother/smoother_1.3.tar.gz", repos = NULL, type = "source")
BiocManager::version()
x=c( "monocle", "destiny", "scRNAseq","celaref","SingleR", "clusterProfiler", 
     "org.Mm.eg.db", "org.Hs.eg.db", "AnnotationHub", "ComplexHeatmap", 'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma',
     'S4Vectors', 'SingleCellExperiment',  'SummarizedExperiment', 'batchelor', 'lme4',  'HDF5Array', "Cardinal", 'ggrastr',
     "limma", "KEGGgraph", "siggenes","BiocParallel", "Rserve","crmn","MSnbase", "multtest","RBGL","fgsea","httr")
x[!x %in% rownames(installed.packages())]
BiocManager::install(x[!x %in% rownames(installed.packages())])

## Install from PyPi
conda create -n scenic-velo-cpdb python=3.10.6 ipykernel=6.15.1 tornado=6.1 typing_extensions=4.4.0 r-base=4.4.3 rpy2
conda activate scenic-velo-cpdb
pip install -i https://pypi.mirrors.ustc.edu.cn/simple  scvelo==0.2.5  scanpy==1.9.6  cellphonedb==5.0.0  MulticoreTSNE==0.1    matplotlib==3.6.0 plotnine==0.10.1   ktplotspy==0.2.1  anndata2ri   # xarray  adjustText  scvi-tools
pip install  -i https://pypi.mirrors.ustc.edu.cn/simple  -r pyscenic121s-req.txt
pip install pyscenic

# Instructions for use
Single-cell and Single-nucleus RNA-seq Analysis:

scripts/
├── linux code.txt
├── scRNA-seq sample main seurat pipeline.R
├── snRNA-seq sample main seurat pipeline.R
├── maternal or fetal status.txt
├── scrublet-python.txt
├── integration.R
├── gene-modification.R
├── sc and snRNA-seq-uamp-projection.R
├── sc and snRNA-seq-cca-merge.R
├── SCENIC.py
├── pseudotime.R
├── cellphoneDB-python.txt
├── hura-cross-species comparative.R
└── mora-cross-species comparative.R

# License
This project is covered under the GPL-3.0 License.









