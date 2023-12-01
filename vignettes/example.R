library(SingleCellExperiment)
library(devtools)
devtools::load_all('../sc_annoTools')

# 预处理过的Seurat类型的数据集
test.obj <- readRDS('../annoTools/data/GSE81608.test.obj')
ref.obj <- readRDS('../annoTools/data/GSE81608.ref.obj')

# 以 SingleR为例
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

# 以 SCINA为例
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