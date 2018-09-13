#' Make Long Data from Wide Data
#'
#' @param wide A data set in wide format.
#'
#' @return A data set in long format.
#'
#' @export
wide2long <- function(wide){

  # Force column names
  if(is.null(colnames(wide))){
    colnames(wide) <- paste0("Col", 1:ncol(wide))
  }

  # Force row names
  if(is.null(rownames(wide))){
    rownames(wide) <- as.character(rownames(wide))
  }

  df <- data.frame("value" = as.vector(as.matrix(wide)))
  df$variable <- unlist(lapply(colnames(wide), function(x) rep(x, nrow(wide))))
  df$id <- rownames(wide)
  return(df)
}

#' Make Wide Data from Long Data
#'
#' @param formula A formula to guide expansion. The y will be the cell.
#'  The first x will be the row. The last x will be the column.
#' @param long A data set in long format.
#'
#' @return A data set in wide format.
#'
#' @export
long2wide <- function(formula, long){

  fchar <- as.character(formula)
  if(length(fchar) != 3) stop("Please provide a y and at least two x to the formula.")
  as.val <- fchar[2] # leftside as value
  rightside <- unlist(strsplit(fchar[3], split = " \\+ "))
  N <- length(rightside)
  if(N < 2) stop("Please provide a y and at least two x to the formula.")
  as.row <- rightside[1:(N-1)]
  as.col <- rightside[N]

  stats::reshape(long, v.names = as.val, idvar = as.row,
                 timevar = as.col, direction = "wide")
}
