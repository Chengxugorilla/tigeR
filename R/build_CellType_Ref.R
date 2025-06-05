#' @title Generate a custom reference matrix
#' @description Generate a custom reference matrix from single-cell RNA sequencing data for TME deconvolution.
#' @param Seurat_obj a Seurat object in which samples are clustered into different cell types.
#' @param logfc.threshold limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.1 Increasing logfc.threshold speeds up the function, but can miss weaker signals. If the slot parameter is "scale.data" no filtering is performed.
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.01.
#' @param only.pos only return positive markers (FALSE by default).
#' @export

build_CellType_Ref <- function(Seurat_obj,
                            logfc.threshold = 0.1,
                            min.pct=0.1,only.pos=TRUE){
  Markers <- Seurat::FindAllMarkers(object = Seurat_obj,
                                    assay = "RNA",
                                    logfc.threshold = 0.15,
                                    min.pct = min.pct,
                                    only.pos = only.pos)
  Markers <- Markers %>%
    dplyr::filter(p_val_adj < 0.01)
  Markers$d <- Markers$pct.1 - Markers$pct.2

  Markers <- Markers %>%
    dplyr::filter(d > 0.1)
  Markers <- Markers %>%
    dplyr::arrange(cluster,dplyr::desc(avg_log2FC))

  used_gene <- sort(unique(Markers$gene))
  tpm <- Seurat::GetAssayData(Seurat_obj, layer = "count") %>%
    as.matrix() %>%
    count2tpm()

  sig.list <-
    lapply(levels(Seurat_obj$celltype), function(x){
      apply(tpm[,Seurat_obj$celltype==x],1,mean)
    })
  names(sig.list) <- levels(Seurat_obj$celltype)
  sig_matrix <- do.call("cbind",sig.list)
  selected_gene <- intersect(rownames(sig_matrix),used_gene)
  return(sig_matrix[selected_gene,])
}
