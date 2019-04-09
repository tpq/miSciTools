#' Comparative Toxicogenomics Database
#'
#' These data catalog chemical-gene interactions curated from the literature.
#'  Each interaction in this database contains useful information about the
#'  source of this knowledge and the type of interaction observed.
#'
#' @format A \code{data.frame} with >1,000,000 chemical-gene interactions (rows)
#'  and 11 annotations (columns). Retrieved on April 9, 2019.
#'
#' @source \url{http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz}
#'
#' @usage data(ctd)
#'
#' @examples
#' \dontrun{
#' library(miSciTools)
#' data(ctd)
#' head(ctd[, c("GeneID", "ChemicalName")])
#' genes <- ctd$GeneID[ctd$ChemicalName == "Platinum"]
#' pvals <-
#'   simpliGSEA(genes = genes, # Test "Platinum" interactions
#'              universe = unique(ctd$GeneID),
#'              annot.genes = ctd$GeneID,
#'              annot.terms = ctd$ChemicalName)
#' head(sort(p.adjust(pvals)))
#' }
"ctd"

#' Comparative Toxicogenomics Database (Wide Format)
#'
#' Unlike \code{data(ctd)}, this object explodes the "InteractionActions"
#'  column into multiple columns, each describing an association type.
#'
#' @format A \code{data.frame} with >1,000,000 chemical-gene interactions (rows)
#'  and 11 annotations (columns). Retrieved on April 9, 2019.
#'
#' @source \url{http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz}
#'
#' @usage data(ctd.wide)
"ctd.wide"

#' KEGGREST Database
#'
#' These data catalog KEGGREST-gene relationships from the KEGGREST database.
#'
#' @format A \code{data.frame} with 27,944 KEGGREST-gene relationships (rows).
#'  Retrieved on July 28, 2017.
#'
#' @source \url{http://rest.kegg.jp/link/pathway/hsa}
#'
#' @usage data(kegg)
"kegg"

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
#' @usage data(linc2go)
#'
#' @examples
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
