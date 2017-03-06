#' Perform Gene Set Enrichment Analysis
#'
#' This function simplifies gene set enrichment analysis for custom annotations.
#'  Do not use this function for annotations that have directed acyclic graph
#'  (DAG) relationships (e.g., GO). Gene set enrichment analysis performed
#'  using Fisher's exact test for count data. This function will automatically
#'  filter non-unique \code{annot.genes} and \code{annot.terms} pairs.
#'
#' @param genes A character vector. The gene IDs of interest.
#' @param universe A character vector. The complete set of gene IDs used to
#'  generate the gene IDs of interest.
#' @param annot.genes A character vector. The gene IDs to which \code{annot.terms}
#'  will correspond.
#' @param annot.terms A character vector. The annotations corresponding to
#'  \code{annot.genes}.
#' @param alternative A character string. Argument passed to \code{fisher.test}.
#' @return A vector of unadjusted p-values named as the tested annotation.
#'
#' @export
simpliGSEA <- function(genes, universe, annot.genes, annot.terms, alternative = "greater"){

  if(length(annot.genes) != length(annot.terms)){

    stop("The 'annot.genes' and 'annot'terms' vectors must have same length!")
  }

  if(class(genes) == "factor") genes <- as.character(genes)
  if(class(universe) == "factor") universe <- as.character(universe)
  if(class(annot.genes) == "factor") annot.genes <- as.character(annot.genes)
  if(class(annot.terms) == "factor") annot.terms <- as.character(annot.terms)

  if(!all(universe %in% annot.genes)){

    message("Some 'universe' not in 'annot.genes'. Pruning unused genes.")
    universe <- universe[universe %in% annot.genes]
  }

  if(!all(genes %in% universe)){

    message("Some 'genes' not found in 'universe'. Pruning unused genes.")
    genes <- genes[genes %in% universe]
  }

  if(length(unique(universe)) < length(universe)){

    message("Removing duplicate IDs in 'universe'.")
    universe <- unique(universe)
  }

  if(length(unique(genes)) < length(genes)){

    message("Removing duplicate IDs in 'genes'.")
    genes <- unique(genes)
  }

  # Subset database to contain only the universe
  annot.universe <- annot.genes %in% universe
  annot.genes <- annot.genes[annot.universe]
  annot.terms <- annot.terms[annot.universe]

  # Initialize result container
  annot.terms.unique <- unique(annot.terms)
  k <- length(annot.terms.unique)
  enrichment <- vector("numeric", k)
  names(enrichment) <- annot.terms.unique
  i <- 1

  for(term in annot.terms.unique){

    # Manage progress bar
    numTicks <- progress(i, k, numTicks)
    i <- i + 1

    # TRUE if has term
    A <- factor(universe %in% annot.genes[annot.terms == term], levels = c(TRUE, FALSE))

    # TRUE if in set
    G <- factor(universe %in% genes, levels = c(TRUE, FALSE))

    # See: https://www.biostars.org/p/111953/
    table <- table(G, A)
    enrichment[term] <- stats::fisher.test(as.matrix(table),
                                           alternative = alternative)$p.value
  }

  return(enrichment)
}
