#' Randomly Assign Class Labels
#'
#' This function randomly assigns class labels.
#'
#' @param n.samples The number of samples to label.
#' @param n.groups The number unique classes.
#' @param prob Argument passed to \code{sample}.
#' @export
randSamps <- function(n.samples, n.groups = 2, prob = NULL){

  sample(letters[1:n.groups], replace = TRUE, size = n.samples, prob = prob)
}

#' Randomly Generate Bootstraps
#'
#' This function samples each sample index (from 1 to \code{n.samples})
#'  into a group (from 1 to \code{b}).
#'
#' @param n.samples The number of samples to group.
#' @param b The number of folds to use.
#' @param percent.include The number of samples to include.
#' @export
randBoots <- function(n.samples, b = 10, percent.include = 67){

  len <- 1:n.samples
  lapply(1:b, function(x) sample(len, size = percent.include/100 * n.samples))
}

#' Randomly Generate Folds
#'
#' This function assigns each sample index (from 1 to \code{n.samples})
#'  into a fold (from 1 to \code{f}).
#'
#' @param n.samples The number of samples to group.
#' @param f The number of folds to use.
#' @export
randFolds <- function(n.samples, f = 10){

  randomize <- sample(1:n.samples)
  splits <- split(randomize, 1:f)
  names(splits) <- NULL
  splits
}
