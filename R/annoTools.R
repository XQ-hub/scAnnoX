#------------------------------------------------------------------------------------------
# Referece-based auto-annotation tools.

#' singleRAnno

#' Automated single-cell annotation using SingleR (Aran et al., 2019).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param ref.obj Reference Seurat object.
#' @param ref.ctype Please specify information about the cell type in the reference object, otherwise: defaults to NULL.
#' i.e. Idents are extracted set to the cell type.
#' @param ... More arguments can be assessed using the SingleR function in the SingleR package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export singleRAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' obj.seu <- singleRAnno(obj.seu, ref.obj)
singleRAnno <- function(obj.seu, ref.obj, ref.ctype = NULL, ...) {
    test.expr <- GetAssayData(obj.seu, slot = "data") %>% as.matrix()
    ref.expr <- GetAssayData(ref.obj, slot = "data") %>% as.matrix()

    if (!is.null(ref.ctype)) {
        label.ref <- ref.obj@meta.data[, ref.ctype]
    } else {
        label.ref <- Idents(ref.obj)
    }
    res.pred <- SingleR::SingleR(test = test.expr, ref = ref.expr, labels = label.ref, ...)
    obj.seu$SingleR <- res.pred$labels
    return(obj.seu)
}

#' seuratAnno

#' Automated single-cell annotation using Seurat (Stuart et al., 2019).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param ref.obj Reference Seurat object.
#' @param ref.ctype Please specify information about the cell type in the reference object, otherwise: defaults to NULL.
#' i.e. Idents are extracted set to the cell type.
#' @param npcs Number of PCs used for prediction. Default: 30.
#' @param ... More arguments can be assessed using the TransferData function in the Seurat package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export seuratAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' obj.seu <- seuratAnno(obj.seu, ref.obj)
seuratAnno <- function(obj.seu, ref.obj, ref.ctype = NULL, npcs = 30, ...) {
    anchors <- FindTransferAnchors(reference = ref.obj, query = obj.seu, dims = 1:npcs)
    if (!is.null(ref.ctype)) {
        label.ref <- ref.obj@meta.data[, ref.ctype]
    } else {
        label.ref <- Idents(ref.obj)
    }
    res.pred <- TransferData(anchorset = anchors, refdata = label.ref, dims = 1:npcs, ...)
    obj.seu$SeuratAnno <- res.pred[, 1]
    return(obj.seu)
}

#' sciBetAnno

#' Automated single-cell annotation using sciBet (Stuart et al., 2019).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param ref.obj Reference Seurat object.
#' @param ref.ctype Please specify information about the cell type in the reference object, otherwise: defaults to NULL.
#' i.e. Idents are extracted set to the cell type.
#' @param ... More arguments can be assessed using the SciBet_R function in the scibetR package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export sciBetAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' obj.seu <- sciBetAnno(obj.seu, ref.obj)
sciBetAnno <- function(obj.seu, ref.obj, ref.ctype = NULL, ...) {
    test.df <- GetAssayData(obj.seu, slot = "data") %>%
        as.data.frame() %>%
        t()
    ref.df <- GetAssayData(ref.obj, slot = "data") %>%
        as.data.frame() %>%
        t()
    if (!is.null(ref.ctype)) {
        label.ref <- ref.obj@meta.data[, ref.ctype]
    } else {
        label.ref <- Idents(ref.obj)
    }
    ref.df <- cbind.data.frame(ref.df, label = label.ref %>% as.vector())
    res.pred <- scibetR::SciBet_R(ref.df, test.df, ...)
    obj.seu$sciBet <- res.pred
    return(obj.seu)
}

#' scmapAnno

# ‘ Automated single-cell annotation using sc-map (Kiselev et al., 2018).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param ref.obj Reference Seurat object.
#' @param ref.ctype Please specify information about the cell type in the reference object, otherwise: defaults to NULL.
#' i.e. Idents are extracted set to the cell type.
#' @param slot.data Extract expression in Seurat slot. Default: "data".
#' @param ... More arguments can be assessed using the scmapCell function in the scmap package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export scmapAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' obj.seu <- scmapAnno(obj.seu, ref.obj)
scmapAnno <- function(obj.seu, ref.obj, ref.ctype = NULL, slot.data = c("data", "count"), ...) {
    system.file(package = "scAnnoX") %>%
        file.path(., "scmap-master") %>%
        load_all(.)

    ref.expr <- GetAssayData(ref.obj, slot = slot.data) %>% as.matrix()
    if (!is.null(ref.ctype)) {
        label.ref <- ref.obj@meta.data[, ref.ctype]
    } else {
        label.ref <- Idents(ref.obj)
    }
    sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(ref.expr)), colData = label.ref)

    if (match.arg(slot.data) == "count") {
        logcounts(sce) <- log2(normcounts(sce) + 1)
    } else {
        logcounts(sce) <- normcounts(sce)
    }
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sce[!duplicated(rownames(sce)), ] %>%
        selectFeatures(.) %>%
        indexCell(.)

    test.expr <- GetAssayData(obj.seu, slot = slot.data) %>% as.matrix()
    sce.test <- SingleCellExperiment(assays = list(normcounts = as.matrix(test.expr)))
    if (match.arg(slot.data) == "count") {
        logcounts(sce.test) <- log2(normcounts(sce.test) + 1)
    } else {
        logcounts(sce.test) <- normcounts(sce.test)
    }
    rowData(sce.test)$feature_symbol <- rownames(sce.test)
    sce.test <- sce.test[!duplicated(rownames(sce.test)), ] %>%
        selectFeatures(.) %>%
        indexCell(.)

    pre.res <- scmapCell(sce.test, list(metadata(sce)$scmap_cell_index))
    scmap.clusters <- scmapCell2Cluster(pre.res, list(as.character(colData(sce)[, 1])))
    obj.seu$scmapAnno <- scmap.clusters$combined_labs
    return(obj.seu)
}

#' chetahAnno

# ‘ Automated single-cell annotation using CHETAH (Kanter et al., 2019).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param ref.obj Reference Seurat object.
#' @param ref.ctype Please specify information about the cell type in the reference object, otherwise: defaults to NULL.
#' i.e. Idents are extracted set to the cell type.
#' @param slot.data Extract expression in Seurat slot. Default: "data".
#' @param ... More arguments can be assessed using the CHETAHclassifier function in the CHETAH package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export chetahAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' obj.seu <- chetahAnno(obj.seu, ref.obj)
chetahAnno <- function(obj.seu, ref.obj, ref.ctype = NULL, slot.data = "data", ...) {
    test.expr <- GetAssayData(obj.seu, slot = slot.data) %>% as.matrix()
    sce.test <- SingleCellExperiment(assays = list(counts = test.expr))
    if (slot.data == "count") sce.test <- scater::logNormCounts(sce.test)

    ref.expr <- GetAssayData(ref.obj, slot = slot.data) %>% as.matrix()
    if (!is.null(ref.ctype)) {
        label.ref <- ref.obj@meta.data[, ref.ctype]
    } else {
        label.ref <- Idents(ref.obj)
    }
    sce.ref <- SingleCellExperiment(
        assays = list(counts = as.matrix(ref.expr)),
        colData = DataFrame(celltypes = label.ref)
    )
    if (slot.data == "count") sce.ref <- scater::logNormCounts(sce.ref)
    pred.res <- CHETAH::CHETAHclassifier(input = sce.test, ref_cells = sce.ref, ...)
    obj.seu$CHETAH <- pred.res$celltype_CHETAH %>% data.frame()
    return(obj.seu)
}

#------------------------------------------------------------------------------------------
# Marker-based auto-annotation tools.

#' scsorterAnno

# ‘ Automated single-cell annotation using scsorter (Guo et al., 2021).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param marker.lst A list contained maker genes for each cell type.
#' @param ... More arguments can be assessed using the scSorter function in the scSorter package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export scsorterAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' marker.lst <- findMarkerToolsForSc(ref.obj, to.list = TRUE, top.k = 30)
#' obj.seu <- scsorterAnno(obj.seu, marker.lst)
scsorterAnno <- function(obj.seu, marker.lst, ...) {
    test.expr <- GetAssayData(obj.seu, slot = "data") %>% as.matrix()
    anno <- collapseLisToFrame(marker.lst)
    pred.res <- scSorter::scSorter(expr = test.expr, anno = anno, ...)
    obj.seu$scSorterAnno <- pred.res$Pred_Type
    return(obj.seu)
}


#' sctypeAnno

#' Automated single-cell annotation using sc.type (Ianevski et al., 2022).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param marker.lst A list contained maker genes for each cell type. Default: NULL.
#' @param tissue.type Please specify the tissue type when using the default marker database wrapped by the sc.type package. Default: Immune system.
#' @param ... More arguments can be assessed using the sctype_score function in the sc.type package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export sctypeAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' marker.lst <- findMarkerToolsForSc(ref.obj, to.list = TRUE, top.k = 30)
#' obj.seu <- sctypeAnno(obj.seu, marker.lst)
sctypeAnno <- function(obj.seu, marker.lst = NULL, tissue.type = "Immune system", ...) {
    system.file(package = "scAnnoX") %>%
        file.path(., "sc-type-master/R/auto_detect_tissue_type.R") %>%
        source(.)
    system.file(package = "scAnnoX") %>%
        file.path(., "sc-type-master/R/gene_sets_prepare.R") %>%
        source(.)
    system.file(package = "scAnnoX") %>%
        file.path(., "sc-type-master/R/sctype_score_.R") %>%
        source(.)

    marker.dbfile <- system.file("modules", "sc-type-master/ScTypeDB_full.xlsx", package = "scAnnoX")
    if (is.null(marker.lst)) {
        gs.list <- gene_sets_prepare(marker.dbfile, tissue.type)
    } else {
        gs.neg <- lapply(names(marker.lst), function(xx) character(0)) %>% `names<-`(names(marker.lst))
        gs.pos <- marker.lst
        gs.list <- list(gs_positive = gs.pos, gs_negative = gs.neg)
    }
    sc.expr <- GetAssayData(obj.seu) %>% as.matrix()
    es.max <- sctype_score(scRNAseqData = sc.expr, gs = gs.list$gs_positive, gs2 = gs.list$gs_negative, ...)
    pred.res <- apply(es.max, 2, which.max) %>% rownames(es.max)[.]
    obj.seu$scTypeAnno <- pred.res
    return(obj.seu)
}


#' cellIDAnno

#' Automated single-cell annotation using CellID (Ianevski et al., 2022).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param marker.lst A list contained maker genes for each cell type.
#' @param gset.len Cell types with a number of markers less than k need to be removed. Default: 5.
#' @param sig.cut Optimal cutoff value to determine significant prediction or not. Default: 2.
#' @param species Specify species of the query Seurat object data when using markes from PanglaoDB database. Default: Hs.
#' @param tissue.type Specify the tissue type when using markes from PanglaoDB database. Default: NULL.
#' @param ... More arguments can be assessed using the sctype_score function in the sc.type package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export cellIDAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' marker.lst <- findMarkerToolsForSc(ref.obj, to.list = TRUE, top.k = 30)
#' obj.seu <- cellIDAnno(obj.seu, marker.lst)
cellIDAnno <- function(obj.seu, marker.lst = NULL, gset.len = 5, sig.cut = 2, species = c("Hs", "Mm"), tissue.type = NULL, ...) {
    assay.names <- obj.seu@assays %>% names()
    if (!("data" %in% assay.names)) obj.seu@assays$data <- obj.seu@assays$RNA
    obj.seu <- CelliD::RunMCA(obj.seu)
    if (is.null(marker.lst)) {
        panglao <- system.file("data", "PanglaoDB_markers_27_Mar_2020.tsv.gz", package = "scAnnoX") %>% read_tsv()
        panglao.all <- panglao %>% filter(str_detect(species, match.arg(species)))
        if (!is.null(tissue.type)) panglao.all <- panglao.all %>% filter(organ == tissue.type)
        panglao.all <- panglao.all %>%
            group_by(`cell type`) %>%
            summarise(geneset = list(`official gene symbol`))
        all.gs <- setNames(panglao.all$geneset, panglao.all$`cell type`)
        all.gs <- all.gs[sapply(all.gs, length) >= gset.len]
    } else {
        all.gs <- marker.lst[sapply(marker.lst, length) >= gset.len]
    }

    HGT.all.gs <- CelliD::RunCellHGT(obj.seu, pathways = all.gs, dims = 1:50)
    all.gs.pred <- rownames(HGT.all.gs)[apply(HGT.all.gs, 2, which.max)]
    all.gs.pred.signif <- ifelse(apply(HGT.all.gs, 2, max) > 2, yes = all.gs.pred, "unassigned")
    obj.seu$CellIDAnno <- all.gs.pred
    obj.seu$CellIDAnno.adj <- all.gs.pred.signif
    return(obj.seu)
}


#' sccatchAnno

#' Automated single-cell annotation using scCATCH (Shao et al., 2020).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param qry.cluster Specify cluster information of query Seurat object in meta.data slot. Default: idents.
#' @param marker.lst A list contained maker genes for each cell type.
#' @param species Species of query Seurat object data which needs to be annotated. Default: Human.
#' @param ... More arguments can be assessed using the findcelltype function in the scCATCH package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export sccatchAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' marker.lst <- findMarkerToolsForSc(ref.obj, to.list = TRUE, top.k = 30)
#' obj.seu <- sccatchAnno(obj.seu, marker.lst)
#' obj.seu <- sccatchAnno(obj.seu, marker.lst = NULL)
sccatchAnno <- function(obj.seu, marker.lst = NULL, qry.cluster = "idents", species = c("Human", "Mouse"), ...) {
    system.file(package = "scAnnoX") %>%
        file.path(., "scCATCH-master") %>%
        load_all(.)

    test.expr <- GetAssayData(obj.seu, slot = "data") %>% as.matrix()
    if (qry.cluster == "idents") {
        clus.info <- Idents(obj.seu)
    } else {
        clus.info <- obj.seu@meta.data[, qry.cluster]
    }
    data <- scCATCH::createscCATCH(data = test.expr, cluster = clus.info %>% as.vector())

    if (is.null(marker.lst)) {
        cellmatch <- scCATCH::cellmatch
        marker.db <- cellmatch[cellmatch$species == match.arg(species), ]
    } else {
        marker.db <- collapseLisToFrame(marker.lst) %>%
            `colnames<-`(c("celltype", "gene")) %>%
            mutate(
                species = "Human",
                tissue = NA,
                cancer = "Normal",
                condition = "Normal",
                subtype1 = NA,
                subtype2 = NA,
                subtype3 = NA,
                resource = NA,
                pmid = NA
            )
    }
    obj <- scCATCH::findmarkergene(data, species = match.arg(species), if_use_custom_marker = TRUE, marker = marker.db)
    pred.res <- scCATCH::findcelltype(obj)
    obj.tmp <- RenameIdents(obj.seu, split(pred.res@celltype$cell_type, pred.res@celltype$cluster) %>% unlist(.))
    obj.seu$scCATCH <- Idents(obj.tmp) %>% as.vector()
    return(obj.seu)
}


#' scinaAnno

# ‘ Automated single-cell annotation using SCINA (Zhang et al., 2018).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param marker.lst A list contained maker genes for each cell type.
#' @param ... More arguments can be assessed using the scSorter function in the SCINA package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export scinaAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' marker.lst <- findMarkerToolsForSc(ref.obj, to.list = TRUE, top.k = 30)
#' obj.seu <- scinaAnno(obj.seu, marker.lst)
scinaAnno <- function(obj.seu, marker.lst, ...) {
    pred.res <- SCINA::SCINA(
        obj.seu %>% GetAssayData(., slot = "data") %>% as.matrix(),
        marker.lst,
        max_iter = 1000,
        convergence_n = 12,
        convergence_rate = 0.999,
        sensitivity_cutoff = 0.9,
        ...
    )
    obj.seu$SCINA <- pred.res$cell_labels
    return(obj.seu)
}
