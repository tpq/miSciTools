#' Linc2GO GeneSetCollection
#'
#' These data expose the Linc2GO database as a \code{GeneSetCollection} object
#'  to use in conjunction with the \code{Category} and \code{GOstats}
#'  packages for gene set enrichment analysis. The gene ontology terms
#'  in this database stem from an application of the competing endogenous
#'  RNA (ceRNA) hypothesis which speculates that lncRNAs regulate expression,
#'  at least in part, by competing with mRNAs for microRNA binding sites.
#'  Since microRNAs bind mRNA transcriptS to silence their expression,
#'  a microRNA-lncRNA complex could free up that mRNA for expression.
#'  Therefore, the curators integrate microRNA-mRNA and microRNA-lncRNA
#'  binding patterns to predict potential lncRNA-mRNA collaborations.
#'  Based on the gene ontology of the mRNA transcript pairs, the curators can
#'  make an educated guess about the gene ontology of the lncRNA.
#'
#' This database uses the lincRNA ID archived by the Human lincRNA Catalog.
#'  Conveniently, these represent the same terms used by UCSC in the
#'  \code{TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts} \code{TxDb} object.
#'
#' @format A \code{GeneSetCollection} object from the \code{GSEABase} package.
#'
#' @source \url{http://www.bioinfo.tsinghua.edu.cn/~liuke/Linc2GO/}
#'
#' @example
#' \dontrun{
#' library(miSciTools)
#' data(linc2go)
#' show(linc2go)
#' universe <- unique(unlist(geneIds(linc2go))) # Pull out all features
#' genes <- geneIds(linc2go)[[2]] # Pull out all GO:0000003 terms
#' params <-
#'   Category::GSEAGOHyperGParams(name = "ExampleGSEA",
#'                                geneSetCollection = linc2go,
#'                                geneIds = genes,
#'                                universeGeneIds = universe,
#'                                ontology = "BP",
#'                                pvalueCutoff = 0.05,
#'                                conditional = FALSE,
#'                                testDirection = "over")
#' Over <- GOstats::hyperGTest(params)
#' head(summary(Over))
#' }
"linc2go"
