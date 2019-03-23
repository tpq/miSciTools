#' Class Check
#'
#' Checks whether an object belongs to a specified class.
#'  For back-end use only.
#'
#' @param x An object.
#' @param what A character vector. The classes any of which \code{x} should have.
#' @param msg A string. An error message if \code{x} is not \code{what}.
classCheck <- function(x, what, msg){

  if(class(x) %in% what){ return(TRUE)
  }else if(any(sapply(what, function(cl) inherits(x, cl)))){ return(TRUE)
  }else{ stop(msg) }
}

#' Select Features by Differential Proportionality Analysis
#'
#' \code{fsPropd} selects features using the \code{propd} function
#'  from the \code{propr} package.
#'
#' @param object An \code{ExprsArray} object to undergo feature selection.
#' @param top A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should undergo feature selection. A character vector
#'  indicates specifically which features by name should undergo feature selection.
#'  Set \code{top = 0} to include all features. A numeric vector can also be used
#'  to indicate specific features by location, similar to a character vector.
#' @param keep A numeric scalar. Specifies the number of top features that should get
#'  returned by the feature selection method. Use of \code{keep} is generally not
#'  recommended, but can speed up analyses of large data.
#' @param modRatios A logical scalar. Toggles whether to compute theta from
#'  the feature ratios as provided. Set \code{modRatios = TRUE} if data were
#'  recasted by a prior \code{modRatios} call. If \code{TRUE}, the \code{alpha}
#'  and \code{weighted} arguments will not work.
#' @param ... Arguments passed to the detailed function.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPropd <- function(object, top = 0, keep = 0, modRatios = FALSE, ...){ # args to propd

  packageCheck("propr")
  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  if(modRatios){

    exprso::fs.(object, top,
                uniqueFx = function(data, outcome, top, ...){

                  # Set up groups and sizes
                  grp1 <- as.character(outcome) == "Control"
                  p1 <- sum(grp1) - 1
                  grp2 <- !grp1
                  p2 <- sum(grp2) - 1
                  p <- p1 + p2 + 1

                  # Calculate theta
                  thetas <- apply(data, 2, function(x){
                    (p1 * stats::var(x[grp1]) + p2 * stats::var(x[grp2])) / (p * stats::var(x))
                  })

                  # Sort results
                  top <- names(thetas[order(thetas)])
                  top
                }, keep, ...)

  }else{

    exprso::fs.(object, top,
                uniqueFx = function(data, outcome, top, ...){

                  # Order pairs by theta
                  pd <- suppressMessages(propr::propd(data, outcome))
                  pd@results <- pd@results[order(pd@results$theta),]

                  # Index features by when they first appear
                  nrows <- nrow(pd@results)
                  index <- floor(seq(1, nrows+.5, .5))
                  odds <- as.logical(1:(nrows*2) %% 2)
                  index[odds] <- index[odds] + nrows
                  join <- c(pd@results$Partner, pd@results$Pair)
                  join <- join[index]

                  # Rank features by first appearance
                  rankedfeats <- unique(join)
                  top[rankedfeats]
                }, keep, ...)
  }
}

#' Select Features by Principal Ratio Analysis
#'
#' \code{fsPRA} selects features using the \code{pra} function
#'  from the \code{propr} package.
#'
#' @inheritParams fsPropd
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPRA <- function(object, top = 0, keep = 0, ...){ # args to selbal

  packageCheck("balance")
  packageCheck("propr")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  exprso::fs.(object, top,
              uniqueFx = function(data, outcome, top, ...){

                args <- as.list(substitute(list(...)))[-1]
                args <- append(list("counts" = data), args)

                # Run pra on data
                res <- do.call(propr::pra, args)

                # Describe log-ratio results as an SBP object
                best <- res$best
                sbp <- matrix(0, ncol(data), nrow(best))
                rownames(sbp) <- colnames(data)
                for(i in 1:nrow(best)){
                  sbp[best$Partner[i],i] <- 1
                  sbp[best$Pair[i],i] <- -1
                }
                colnames(sbp) <- paste0("z", 1:ncol(sbp))

                # Re-compute log-ratios as balances
                sbp <- balance::sbp.subset(sbp, ratios = TRUE)
                balances <- t(balance::balance.fromSBP(data, sbp))
                colnames(balances) <- rownames(data)
                class(sbp) <- "SBP"

                list(balances,
                     sbp)

              }, keep, ...)
}

#' Select a Discriminative Balance
#'
#' \code{fsSelbal} selects features using the \code{selbal} function
#'  from the \code{selbal} package.
#'
#' @inheritParams fsPropd
#' @return Returns an \code{ExprsArray} object.
#' @export
fsSelbal <- function(object, top = 0, keep = 0, ...){ # args to selbal

  packageCheck("balance")
  packageCheck("selbal")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(class(object) == "RegrsArray"){

    exprso::fs.(object, top,
                uniqueFx = function(data, outcome, top, ...){

                  args <- as.list(substitute(list(...)))[-1]
                  args <- append(list("x" = data,
                                      "y" = outcome),
                                 args)

                  # Run selbal on data
                  res <- do.call(selbal::selbal, args)

                  # Build an SBP from selbal object
                  n <- as.character(res[[6]]$NUMERATOR)
                  d <- as.character(res[[6]]$DENOMINATOR)
                  sbp <- matrix(0, ncol(data), 1)
                  rownames(sbp) <- colnames(data)
                  sbp[n[n!="-"],] <- 1
                  sbp[d[d!="-"],] <- -1
                  colnames(sbp) <- "z1"
                  class(sbp) <- "SBP"

                  # Compute balances
                  balances <- t(balance::balance.fromSBP(data, sbp))
                  colnames(balances) <- rownames(data)

                  list(balances,
                       sbp)
                }, keep, ...)

  }else if(class(object) %in% c("ExprsBinary", "ExprsMulti")){

    exprso::fs.(object, top,
                uniqueFx = function(data, outcome, top, ...){

                  args <- as.list(substitute(list(...)))[-1]
                  args <- append(list("x" = data,
                                      "y" = factor(outcome, levels = c("Control", "Case"))),
                                 args)

                  # Run selbal on data
                  res <- do.call(selbal::selbal, args)

                  # Build an SBP from selbal object
                  n <- as.character(res[[6]]$NUMERATOR)
                  d <- as.character(res[[6]]$DENOMINATOR)
                  sbp <- matrix(0, ncol(data), 1)
                  rownames(sbp) <- colnames(data)
                  sbp[n[n!="-"],] <- 1
                  sbp[d[d!="-"],] <- -1
                  colnames(sbp) <- "z1"
                  class(sbp) <- "SBP"

                  # Compute balances
                  balances <- t(balance::balance.fromSBP(data, sbp))
                  colnames(balances) <- rownames(data)

                  list(balances,
                       sbp)
                }, keep, ...)
  }
}

#' Select Features by Recursive Feature Elimination
#'
#' \code{fsPathClassRFE} selects features using the \code{fit.rfe} function
#'  from the \code{pathClass} package.
#'
#' @inheritParams fsPropd
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPathClassRFE <- function(object, top = 0, keep = 0, ...){ # args to fit.rfe

  packageCheck("pathClass")
  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  exprso::fs.(object, top,
              uniqueFx = function(data, outcome, top, ...){

                # Set up "make.names" key for improper @exprs row.names
                labels <- factor(outcome, levels = c("Control", "Case"))
                key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

                # NOTE: RFE as assembled by pathClass is via a linear kernel only
                # NOTE: By default, fit.rfe iterates through C = 10^c(-3:3)
                # Run fit.rfe()
                rfe <- pathClass::fit.rfe(x = data, y = labels, ...)

                # Use "make.names" key to return to original row.names
                final <- merge(data.frame("new" = rfe$features), key, sort = FALSE)$old
                if(length(final) < 2) stop("Uh oh! fsPathClassRFE did not find enough features!")
                as.character(final)
              }, keep, ...)
}
