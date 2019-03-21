#' Reconstruct Data from RDA Object
#'
#' This function reconstructs the original data from the
#'  \code{$CA} slot of a \code{vegan::rda} object. This should
#'  work for partial RDA also, returning the data "X" that
#'  excludes the constrained part "Z".
#'
#' @param rda An object from \code{vegan::rda}.
#' @return Returns a \code{data.frame} object.
#' @export
rda.reconstruct <- function(rda){

  if(!"rda" %in% class(rda)){
    stop("This function expects an rda object from vegan::rda().")
  }

  zz <- rda$CA
  CONSTANT <- sqrt(nrow(zz$u) - 1) * sqrt(zz$eig)
  U <- t(apply(zz$u, 1, function(x) x*CONSTANT)) # rda scales zz$u
  xhat <- U %*% t(zz$v) # zz$v is same as rotation
  mu <- attr(zz$Xbar, "scaled:center") # original mean center
  final <- scale(xhat, center = -mu, scale = FALSE)
  as.data.frame(final)
}
