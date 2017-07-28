#' Run edgeR Exact Test
#'
#' Run edgeR Exact Test for differential expression analysis.
#'
#' @param counts A data.frame of raw sequence counts.
#' @param group A character vector of group labels.
#'
#' @export
quick.edgeR <- function(counts, group){

  if(!requireNamespace("edgeR", quietly = TRUE)){
    stop("Uh oh! This method depends on edgeR. ",
         "Try running: miSciTools::demand('edgeR')")
  }

  y <- edgeR::DGEList(counts = counts, group = group)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateCommonDisp(y)
  y <- edgeR::estimateTagwiseDisp(y)
  et <- edgeR::exactTest(y)
  tt <- as.data.frame(edgeR::topTags(et, n = nrow(et)))
  tt.fdr <- tt[tt$FDR < .05,]
  edgeR::plotSmear(et, de.tags = rownames(tt.fdr), cex = 0.5)
  return(tt.fdr)
}
