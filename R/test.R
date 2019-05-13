#' Test All Equal Within List
#'
#' This function tests whether all elements in a list are identical.
#'  Works like an iterative \code{all.equal}.
#'
#' @param list A list.
#' @export
test.lequal <- function(list){

  for(i in 1:length(list)){
    if(!identical(list[1], list[i])){
      stop("Element", i, "does not match the first.")
    }}
}
