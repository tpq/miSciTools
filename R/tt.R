#' Make a Table of t-test Results
#'
#' @param df A tidy long-form \code{data.frame}.
#' @param y A string. The column name of the outcome variable.
#' @param x A string. The column name of the contrast variable.
#' @param subset A boolean of length \code{nrow(df)}. How to
#'  subset the \code{data.frame} prior to \code{t.test}.
#' @param digits An integer. The number of significant digits.
#' @param ... Arguments passed to \code{t.test}.
#'
#' @return A table of t-test results.
#'
#' @export
tttest <- function(df, y, x, subset, digits = 2, ...){

  df.sub <- df[subset,]
  uniq <- unique(df[,x])
  formula1 <- stats::as.formula(paste0(y, "~", x))

  li <- lapply(1:length(uniq), function(i){
    lj <- lapply(1:length(uniq), function(j){
      if(i != j){
        x.i <- uniq[i]
        x.j <- uniq[j]
        a <- df.sub[df.sub[,x] == x.i, y]
        b <- df.sub[df.sub[,x] == x.j, y]
        tt <- stats::t.test(a, b, ...)
        paste(format(tt$conf.int, scientific = FALSE, digits = digits), collapse = " to ")
      }else{
        "---"
      }
    })
    as.matrix(lj)
  })

  final <- do.call("cbind", li)
  colnames(final) <- paste(uniq, "vs.")
  rownames(final) <- uniq
  final
}

#' Make a Table of Wilcoxon Rank Sum test Results
#'
#' @param df A tidy long-form \code{data.frame}.
#' @param y A string. The column name of the outcome variable.
#' @param x A string. The column name of the contrast variable.
#' @param subset A boolean of length \code{nrow(df)}. How to
#'  subset the \code{data.frame} prior to \code{t.test}.
#' @param digits An integer. The number of significant digits.
#' @param ... Arguments passed to \code{t.test}.
#'
#' @return A table of Wilcoxon Rank Sum test results.
#'
#' @export
wwtest <- function(df, y, x, subset, digits = 2, ...){

  df.sub <- df[subset,]
  uniq <- unique(df[,x])
  formula1 <- stats::as.formula(paste0(y, "~", x))

  li <- lapply(1:length(uniq), function(i){
    lj <- lapply(1:length(uniq), function(j){
      if(i != j){
        x.i <- uniq[i]
        x.j <- uniq[j]
        a <- df.sub[df.sub[,x] == x.i, y]
        b <- df.sub[df.sub[,x] == x.j, y]
        tt <- stats::wilcox.test(a, b, conf.int = TRUE, ...)
        paste(format(tt$conf.int, scientific = FALSE, digits = digits), collapse = " to ")
      }else{
        "---"
      }
    })
    as.matrix(lj)
  })

  final <- do.call("cbind", li)
  colnames(final) <- paste(uniq, "vs.")
  rownames(final) <- uniq
  final
}
