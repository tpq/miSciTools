#' Replace NA Values
#'
#' This function replaces all NAs in a vector or \code{data.frame}.
#'  If the vector (or column) is numeric, NAs get replaced by the median.
#'  If the vector (or column) is a character or factor, NAs get
#'  replaced by the mode.
#'
#' @param x A vector or \code{data.frame}.
#' @export
replaceNA <- function(x){

  if(class(x) == "data.frame"){

    # Recurse through each column in list...
    return(as.data.frame(lapply(x, replaceNA)))

  }else{

    if(any(is.na(x))){

      # Replace with MEDIAN if numeric
      # Replace with MODE if factor
      # Else ERROR
      if(class(x) == "numeric"){

        message("Numeric: Using median replacement.")
        replace <- stats::median(x[!is.na(x)])

      }else if(class(x) == "factor" | class(x) == "character"){

        message("Factor: Using mode replacement.")
        replace <- names(sort(table(x[!is.na(x)]), decreasing = TRUE)[1])

      }else{

        stop("This function expects a vector with NA values.")
      }

      x[is.na(x)] <- replace

    }else{

      message("No NAs to replace!")
    }

    return(x)
  }
}
