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
#' @return A vector of unadjusted p-values named as the tested annotation.
#'
#' @export
simpliGSEA <- function(genes, universe, annot.genes, annot.terms){

  if(length(annot.genes) != length(annot.terms)){

    stop("The 'annot.genes' and 'annot'terms' vectors must have same length!")
  }

  if("factor" %in% c(class(genes), class(universe), class(annot.genes), class(annot.terms))){

    stop("This function will not accept any factors as input!")
  }

  if(!all(universe %in% annot.genes)){

    message("Some 'universe' not in 'annot.genes'. Pruning unused genes.")
    universe <- universe[universe %in% annot.genes]
  }

  if(!all(genes %in% universe)){

    message("Some 'genes' not found in 'universe'. Pruning unused genes.")
    genes <- genes[genes %in% universe]
  }

  if(length(unique(genes)) < length(genes)){

    cat("Removing duplicate IDs in 'genes'...")
    genes <- unique(genes)
  }

  if(length(unique(universe)) < length(universe)){

    cat("Removing duplicate IDs in 'universe'...")
    universe <- unique(universe)
  }

  # Subset database to contain only the universe
  annot.universe <- annot.genes %in% universe
  annot.genes <- annot.genes[annot.universe]
  annot.terms <- annot.terms[annot.universe]

  # Initialize result container
  enrichment <- vector("numeric", length(unique(annot.terms)))
  names(enrichment) <- unique(annot.terms)

  for(term in unique(annot.terms)){

    # TRUE if has term
    A <- factor(universe %in% annot.genes[annot.terms == term], levels = c(TRUE, FALSE))

    # TRUE if in set
    G <- factor(universe %in% genes, levels = c(TRUE, FALSE))

    # See: https://www.biostars.org/p/111953/
    table <- table(G, A)
    print(table)

    enrichment[term] <- stats::fisher.test(as.matrix(table))$p.value

    cat(term, "Fisher's exact test p-value:", enrichment[term], "\n\n")
  }

  return(enrichment)
}
