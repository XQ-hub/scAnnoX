#' @title autoAnnoTools

#' @description Automated single-cell annotation using publicly available tools.
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param ref.obj Seurat object, only when used with refernce-based tools. Default: NULL.
#' @param marker.lst A list contained maker genes for each cell type for marker based tool. 
#' Notably, ref.obj and marker.lst are antagonistic, and must be provide one. Default: NULL.
#' @param method  A vector of automated annotation tools.
#' @param select.marker Specify method for inferring markers for each subset. Default: Seurat.
#' @param top.k Top k expressed genes of each subset remained. Default: NULL.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export autoAnnoTools
#'  

autoAnnoTools <- function(
	obj.seu, 
	ref.obj = NULL, 
	ref.ctype = NULL,
	marker.lst = NULL, 
	method = c('SingleR', 'Seurat', 'sciBet', 'scmap', 'CHETAH', 'scSorter', 'sc.type', 'cellID', 'scCATCH', 'SCINA'),
	select.marker = c('Seurat'),
	top.k = 30,
	strategy = c('refernce-based', 'marker-based'),
	...
) {
	multipleProcess(10)
	method.sc <- match.arg(method)
	if (match.arg(strategy) == 'marker-based' && is.null(marker.lst)) {
		if (!is.null(ref.obj) && !is.null(ref.ctype)) Idents(ref.obj) <- ref.obj@meta.data[, ref.ctype]
		if (!is.null(ref.obj)) marker.lst <- findMarkerToolsForSc(ref.obj, to.list = TRUE, top.k = top.k)
		else println('Marker list cannot be accessed', status = 'ERROR')	
	}
	if (match.arg(strategy) == 'refernce-based' && is.null(ref.obj)) {
		println('Reference Seurat object must be assigned', status = 'ERROR')
	}
	obj.seu <- switch(
		method.sc,				
		SingleR = singleRAnno(obj.seu, ref.obj, ref.ctype, ...),
		Seurat = seuratAnno(obj.seu, ref.obj, ref.ctype, ...),
		sciBet = sciBetAnno(obj.seu, ref.obj, ref.ctype, ...),
		scmap = scmapAnno(obj.seu, ref.obj, ref.ctype, ...), 
		CHETAH = chetahAnno(obj.seu, ref.obj, ref.ctype, ...),
		scSorter = scsorterAnno(obj.seu, marker.lst, ...),
		scCATCH =  sccatchAnno(obj.seu, marker.lst, qry.cluster = 'idents', ...),
		cellID = cellIDAnno(obj.seu, marker.lst, ...),
		sc.type = sctypeAnno(obj.seu, marker.lst, ...),
		SCINA = scinaAnno(obj.seu, marker.lst, ...)
	)
    return(obj.seu)
}


#' @title autoAnnoResult

#' @description Integrate the prediction results of the algorithm in the autoAnnoTools function.
#' @param anno The predicted results of the algorithm in the autoAnnoTools function.
#' @return Seurat object with consolidated results.
#' @export autoAnnoResult
#' 
 
autoAnnoResult <- function(anno) {
    if(! type(anno) == 'data.frame') anno <- as.data.frame(anno)
    anno$scAnnoX <- apply(anno.res, 1, function(row) names(which.max(table(row)/ncol(anno))))
    return(anno)
    }


#' @title listToolMethods

#' @description Select the direction function of the single cell annotation tool.
#' @return A vector functions of tools used for single-cell annotation.
#' @export listToolMethods
#' 

listToolMethods <- function(){
	return(
		SingleR = singleRAnno,
		Seurat = seuratAnno,
		sciBet = sciBetAnno,
		scmap = scmapAnno, 
		CHETAH = chetahAnno,
		scSorter = scsorterAnno,
		scCATCH =  sccatchAnno,
		cellID = cellIDAnno,
		sc.type = sctypeAnno,
		SCINA = scinaAnno
	)
}