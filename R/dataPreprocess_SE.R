#' @title Prepare expression matrix for down string analysis
#' @description dataPreprocess will remove missing genes. Then returns the sub-matrix of the genes whose SYMBOLs are in the signature.
#' @param SE d
#' @param Signature The aiming gene set(only Gene SYMBOL allowed).
#' @param turn2HL If TRUE, the expression value of a gene is divided to "HIGH" or "LOW" based on its median expression.
#' @export

dataPreprocess_SE <- function(SE, Signature = NULL, turn2HL = TRUE){
  if(is.null(Signature))
    Signature <- rownames(SE)

  genes <- S4Vectors::intersect(Signature, rownames(SE))

  absent_genes <- length(Signature) - length(genes)
  if(absent_genes)
    message(paste0(absent_genes," Signature genes are not found in expression matrix. The function can execute properly, but the performance of the model may be compromised."))

  exp_mtr <- assay(SE[genes,,drop=FALSE])
  idx_R <- which(SE$response_NR=="R")
  idx_N <- which(SE$response_NR=="N")

  exp_mtr <- t(apply(exp_mtr, 1, function(x, Batch) {
    for (i in unique(Batch)) {
      tmp <- x[which(Batch == i)]
      isNA <- is.na(tmp)
      if (all(isNA))
        next
      if (all(tmp[!isNA] == 0))
        tmp[!isNA] <- NA
      x[which(Batch == i)] <- tmp
    }
    x
  }, Batch = SE$dataset_id))

  if(is.vector(exp_mtr))
    exp_mtr <- matrix(exp_mtr, nrow = 1)

  if(turn2HL){
    threshold <-
      apply(exp_mtr,1,function(x){
        mean(mean(x[idx_R],na.rm=TRUE),mean(x[idx_N],na.rm=TRUE),na.rm=TRUE)
      })
    exp_mtr <- t(
      apply(exp_mtr,1,function(x){
        thres <- mean(mean(x[idx_R],na.rm=TRUE),mean(x[idx_N],na.rm=TRUE),na.rm=TRUE)
        x <- ifelse(x>=thres,'HIGH','LOW')

      }))
    if(is.vector(exp_mtr))
      exp_mtr <- matrix(exp_mtr, nrow = 1)
  }else{
    exp_mtr[is.na(exp_mtr)] <- 0
  }

  SE_final <-
    SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(exp_mtr),
                                               colData=S4Vectors::DataFrame(SummarizedExperiment::colData(SE)),
                                               checkDimnames=TRUE)
  if(turn2HL)
    return(list(SE_final, threshold))
  else
    return(exp_mtr)
}
