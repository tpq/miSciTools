#' Remove GEOM From \code{ggplot} Object
#'
#' This function removes a GEOM from a \code{ggplot} object,
#'  something that is useful in conjunction with automated
#'  \code{gg} tools like \code{ggtern}, \code{ggbiplot},
#'  or \code{ggvegan}.
#'
#' Use like \code{remove_geom(g1, "GeomText")}.
#'
#' @author Kamil Slowikowski from Stack Overflow.
#'
#' @param ggplot2_object A \code{ggplot} object.
#' @param geom_type A string. The GEOM.
#' @export
remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}
