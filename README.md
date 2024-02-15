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

# SessionInfo
```r
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22000)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10       compiler_4.2.2    later_1.3.1       urlchecker_1.0.1  prettyunits_1.1.1
 [6] profvis_0.3.8     remotes_2.4.2.1   tools_4.2.2       digest_0.6.31     pkgbuild_1.4.2   
[11] pkgload_1.3.2.1   memoise_2.0.1     lifecycle_1.0.3   rlang_1.1.1       shiny_1.7.4.1    
[16] cli_3.6.0         rstudioapi_0.15.0 curl_5.0.0        fastmap_1.1.1     stringr_1.5.0    
[21] desc_1.4.2        fs_1.6.2          htmlwidgets_1.6.2 vctrs_0.6.3       devtools_2.4.5   
[26] rprojroot_2.0.3   glue_1.6.2        R6_2.5.1          processx_3.8.2    sessioninfo_1.2.2
[31] callr_3.7.3       purrr_1.0.1       magrittr_2.0.3    ps_1.7.5          promises_1.2.0.1 
[36] ellipsis_0.3.2    htmltools_0.5.5   usethis_2.2.2     mime_0.12         xtable_1.8-4     
[41] httpuv_1.6.11     stringi_1.7.12    miniUI_0.1.1.1    cachem_1.0.7      crayon_1.5.2 
```