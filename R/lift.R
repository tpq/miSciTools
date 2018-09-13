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
#' @param object A \code{GRanges} object.
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

  if(!requireNamespace("R.utils", quietly = TRUE)){
    stop("Uh oh! This method depends on R.utils. ",
         "Try running: miSciTools::demand('R.utils')")
  }

  if(!requireNamespace("rtracklayer", quietly = TRUE)){
    stop("Uh oh! This method depends on rtracklayer. ",
         "Try running: miSciTools::demand('rtracklayer')")
  }

  if(!requireNamespace("GenomeInfoDb", quietly = TRUE)){
    stop("Uh oh! This method depends on GenomeInfoDb. ",
         "Try running: miSciTools::demand('GenomeInfoDb')")
  }

  if(!all(GenomeInfoDb::genome(object) %in% from)){

    warning("Check 'genome(object)' to make sure you selected correct 'over.chain' directory.")
  }

  over.chain <- paste0(from, "To", R.utils::capitalize(to), ".over.chain.gz")
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

    if(!requireNamespace("biovizBase", quietly = TRUE)){
      stop("Uh oh! This method depends on biovizBase. ",
           "Try running: miSciTools::demand('biovizBase')")
    }

    cat("Transforming GRangesList to GRanges...\n")
    lifted <- biovizBase::flatGrl(lifted)
    return(lifted)

  }else{

    cat("Returning GRangesList...\n")
    return(lifted)
  }
}
