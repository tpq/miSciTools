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

#' Run DESeq2 LRT or Wald Test
#'
#' Run DESeq2 LRT or Wald Test for differential expression analysis.
#'
#' @param x A data.frame of raw sequence counts.
#' @param y A data.frame of annotations.
#' @param formula The formula used for differential expression analysis.
#' @param test The test to run. Choose between "LRT" and "Wald".
#'
#' @export
quick.DESeq2 <- function(x, y, formula, test = "Wald"){

  if(!requireNamespace("DESeq2", quietly = TRUE)){
    stop("Uh oh! This method depends on DESeq2. ",
         "Try running: miSciTools::demand('DESeq2')")
  }

  # Remove NAs
  f <- as.character(formula)
  f <- f[length(f)]
  terms <- unlist(strsplit(f, " \\+ "))
  print("Terms:")
  print(terms)
  termsNA <- apply(y[, terms, drop = FALSE], 2, is.na)
  keep <- rowSums(termsNA) == 0
  print("NAs removed:")
  print(sum(!keep))

  # Build DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(x[keep,] + 1),
                                        colData = y[keep,],
                                        design = formula)

  # Run DESeq2 analysis
  if(test == "LRT"){

    terms <- terms[length(terms)]
    reduced <- stats::update(formula, stats::as.formula(paste0("~.-", terms)))
    print("Reduced:")
    print(reduced)
    res <- DESeq2::DESeq(dds, test = "LRT", reduced = reduced)

  }else if(test == "Wald"){

    res <- DESeq2::DESeq(dds)

  }else{

    stop("Provided 'test' not supported.")
  }

  # Extract results
  res <- DESeq2::results(res)
  res[is.na(res$padj), "padj"] <- 1
  print(utils::head(res))
  print("Significant findings:")
  print(sum(res$padj < .05))
  as.data.frame(res)
}
