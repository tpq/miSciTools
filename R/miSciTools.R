###########################################################
### Functions to assist in genomic annotations

#' Lift Over Genome Build
#'
#' This function simplifies the process of "lifting over" coordinates from one
#'  genome build to another by automatically downloading the necessary
#'  \code{.over.chain.gz} file.
#'
#' This function will convert any non-UCSC coordinates to UCSC equivalents.
#'  Make sure the arguments \code{from} and \code{to} refer to valid UCSC
#'  genome builds. Feature names will get lost during 'lift over'. To retain,
#'  these labels, consider storing names as a separate column in the
#'  \code{GRanges} object.
#'
#' @param A \code{GRanges} object.
#' @param from A character string. The initial UCSC equivalent genome build
#'  (e.g., hg18).
#' @param to A character string. The target UCSC equivalent genome build
#'  (e.g., hg19).
#' @param flatGrl A logical scalar. Toggles whether to collapse the
#'  default \code{GRangesList} output to a \code{GRanges} output.
#' @return A \code{GRangesList} object or \code{GRanges} object.
#'
#' @export
lift <- function(object, from = "hg18", to = "hg19", flatGrl = TRUE){

  if(!all(genome(object) %in% from)){

    warning("Check 'genome(object)' to make sure you selected correct 'over.chain' directory.")
  }

  over.chain <- paste0(from, "To", capitalize(to), ".over.chain.gz")
  cat("Downloading", over.chain, "as temporary file...\n")
  url <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/", from, "/liftOver/", over.chain)
  file.gz <- tempfile(fileext = ".over.chain.gz")
  utils::download.file(url, file.gz)

  cat("Unzipping", over.chain, "temporary file...\n")
  file.oc <- tempfile(fileext = ".over.chain")
  R.utils::gunzip(filename = file.gz, destname = file.oc)

  cat("Importing unzipped file...\n")
  chain <- rtracklayer::import.chain(file.oc)

  cat("Convering coordinates to UCSC equivalents...\n")
  tempStyle <- GenomeInfoDb::seqlevelsStyle(object)
  GenomeInfoDb::seqlevelsStyle(object) <- "UCSC"

  cat("Using over.chain to 'lift over'...\n")
  lifted <- rtracklayer::liftOver(object, chain)
  GenomeInfoDb::genome(lifted) <- to

  cat("Returning coordinates to original style...\n")
  GenomeInfoDb::seqlevelsStyle(object) <- tempStyle

  if(flatGrl){

    cat("Transforming GRangesList to GRanges...\n")
    lifted <- biovizBase::flatGrl(lifted)
    return(lifted)

  }else{

    cat("Returning GRangesList...\n")
    return(lifted)
  }
}

#' Perform Gene Set Enrichment Analysis
#'
#' This function simplifies gene set enrichment analysis for custom annotations.
#'  Do not use this function for annotations that have directed acyclic graph
#'  (DAG) relationships (e.g., GO). Gene set enrichment analysis performed
#'  using Fisher's exact test for count data. This function will implicitly
#'  filter non-unique \code{annot.genes} and \code{annot.terms} pairs.
#'
#' @param genes A character vector. The gene IDs of interest.
#' @param universe A character vector. The complete set of gene IDs used to
#'  generate the gene IDs of interest.
#' @param annot.genes A character vector. The gene IDs to which \code{annot.terms}
#'  will correspond.
#' @param annot.terms A character vector. The annotations corresponding to
#'  \code{annot.genes}.
#' @return A vector of p-values named based on the tested annotation.
#'
#' @export
simpliGSEA <- function(genes, universe, annot.genes, annot.terms){

  if(length(annot.genes) != length(annot.terms)){

    stop("The 'annot.genes' and 'annot'terms' vectors must have same length!")
  }

  if("factor" %in% c(class(genes), class(universe), class(annot.genes), class(annot.terms))){

    stop("This function will not accept any factors as input!")
  }

  if(!all(genes %in% universe)){

    stop("Some provided 'genes' not found in 'universe'!")
  }

  if(!all(universe %in% annot.genes)){

    stop("Some 'universe' not in 'annot.genes'!")
  }

  if(length(unique(genes)) < length(genes)){

    cat("Removing duplicate IDs in 'genes'...")
    genes <- unique(genes)
  }

  if(length(unique(universe)) < length(universe)){

    cat("Removing duplicate IDs in 'universe'...")
    universe <- unique(universe)
  }

  # No need to test terms that do not show up in universe
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

    enrichment[term] <- fisher.test(as.matrix(table))$p.value

    cat(term, "Fisher's exact test p-value:", enrichment[term], "\n\n")
  }

  return(enrichment)
}

###########################################################
### Functions to assist in data visualization

#' Plot Multiple Graphs
#'
#' Easily plot multiple graphs within the same window. Code adapted from
#'  http://www.cookbook-r.com/.
#'
#' @param ... Multiple plots.
#' @param cols A numeric scalar. The number of plot columns.
#'
#' @export
multiplot <- function(..., cols = 1){

  library(grid)

  # Make a list of plots
  plots <- list(...)
  numPlots <- length(plots)

  # Layout the panel
  layout <- matrix(seq(from = 1,
                       to = cols * ceiling(numPlots/cols)),
                   ncol = cols,
                   nrow = ceiling(numPlots/cols)
  )

  if(numPlots == 1){

    print(plots[[1]])

  }else{

    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow(layout),
                                   ncol(layout))
      )
    )

    # Place each plot, in the correct location
    for(i in 1:numPlots){

      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(
        plots[[i]],
        vp = grid::viewport(layout.pos.row = matchidx$row,
                            layout.pos.col = matchidx$col)
      )
    }
  }

  return(TRUE)
}

###########################################################
### Functions to assist in cluster computing

#' Load or Install Package
#'
#' This function finds and installs a package, if not already installed.
#'
#' @param package A character string. The package to load or install.
#'
#' @export
demand <- function(package){

  if(!is.element(package, installed.packages()[, 1])){

    try({

      cat("Looking for package at CRAN...\n")
      install.packages(package)
    })
  }

  if(!is.element(package, installed.packages()[, 1])){

    cat("Looking for package at Bioconductor...\n")
    source("https://bioconductor.org/biocLite.R")
    biocLite(package, suppressUpdates = TRUE)
  }

  if(!is.element(package, installed.packages()[, 1])){

    try({

      cat("Looking for package at R-Forge...\n")
      install.packages(package, repos = "http://R-Forge.R-project.org")
    })
  }

  if(is.element(package, installed.packages()[,1])){

    library(package, character.only = TRUE)

  }else{

    stop("Could not find requested package!")
  }
}

#' Calculate Peak RAM Used
#'
#' This function calculates the peak amount of RAM used during a function call.
#'
#' @param fx A function to call (e.g., \code{peakRAM(function() phit(counts))}).
#'
#' @export
peakRAM <- function(fx){

  start <- gc(reset = TRUE)
  start <- start["Vcells", 6]
  run <- fx
  run()
  end <- gc(TRUE)
  end <- end["Vcells", 6]
  final <- end - start
  names(final) <- "max.Mb"
  gc(reset = TRUE)
  return(final)
}

#' Write R Script
#'
#' This function saves a character string as an R script. By default, this
#'  script gets stored in a temporary directory as a temporarily file.
#'  This function returns the file path for the saved R script for []
#'  use (e.g., via \code{\link{qsub}}).
#'
#' @param R A character string. The R script to save.
#' @param folder A character string. The folder in which to save the script.
#'  Defaults to a temporary directory name.
#' @param file A character string. The file name to assign to the script.
#'  Defaults to a temporary file name.
#' @return The file path for the saved R script.
#'
#' @export
writeR <- function(R, folder = tempdir(), file = paste0(basename(folder), ".R")){

  script <- paste0(folder, "/", file)
  file.create(script)
  fileConn <- file(script)
  writeLines(R, fileConn)
  close(fileConn)

  return(script)
}

#' Qsub Linux Command
#'
#' This function sends a Linux command to the PBS queue via \code{qsub}.
#'  Specifically,
#'
#' @param cmd A character string. The Linux command to \code{qsub}.
#' @param ... Any additional PBS argument(s). Each argument should
#'  get named according to the bash character argument (e.g., set
#'  \code{-M thom@tpq.me} with \code{qsub(command, M = "thom@tpq.me")}).
#'
#' @export
qsub <- function(cmd, ...){

  # Prepare #PBS args
  args <- as.list(substitute(list(...)))[-1]
  args.bash <- paste0("#PBS -", names(args), " ", unlist(args), "\n",
                      collapse = "")

  # Write script to execute cmd
  bash <- paste0("#! /bin/bash", "\n",
                 args.bash, "\n", cmd)

  # Send script to queue
  run <- paste0("echo \"", bash, "\" | qsub")
  system(run)

  return(TRUE)
}
