#' Linc2GO GeneSetCollection
#'
#' These data [placeholder].
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
NULL
