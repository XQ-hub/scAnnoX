# scAnnoX

# scAnnoX - An R Package Integrating Multiple Public Tools for Single-Cell Annotation
scAnnoX is an R package that integrates 10 different cell identity algorithms for single-cell sequencing data into a unified framework, and allows comparison between different algorithms according to experimental results. Among them, 10 algorithms include SingleR, Seurat, sciBet, scmap, CHETAH, scSorter, sc.type, cellID, scCATCH, SCINA. It is designed to assist researchers to analyze scRNA-seq data more efficiently, so that they can make more informed decisions among the complex selection of single cell identification algorithms, and more easily test, evaluate, and compare multiple algorithms in an integrated environment. The scAnnoX package is available at https://github.com/XQ-hub/scAnnoX.

# Install

```R
devtools::install_github('XQ-hub/scAnnoX')
```

# Usage

```R
library(scAnnoX)
library(SingleCellExperiment)

# A preprocessed dataset of type Seurat.
test.obj <- readRDS('../scAnnoX/data/GSE81608.test.obj')
ref.obj <- readRDS('../scAnnoX/data/GSE81608.ref.obj')

# Using SingleR as an example.
pred.obj <- autoAnnoTools(
    test.obj,
    ref.obj = ref.obj,
    ref.ctype = 'ActureAnno',
    marker.lst = marker.lst,
    method = 'SingleR',
    select.marker = 'Seurat',
    top.k = 30,
    strategy = 'refernce-based'
) 

# Using SCINA as an example.
Idents(ref.obj) <- ref.obj$CellType
marker.lst <- findMarkerToolsForSc(ref.obj, to.list = TRUE, top.k = 30)
pred.obj <- autoAnnoTools(
    test.obj,
    ref.obj = ref.obj,
    ref.ctype = 'CellType',
    marker.lst = marker.lst,
    method = 'SCINA',
    select.marker = 'Seurat',
    top.k = 30,
    strategy = 'marker-based'
) 
```

# Contributors

scAnnoX was developed by Xiaoqian Huang. Please contact Xiaoqian Huang for any questions or suggestions.

# scAnnoX
# scAnnoX
