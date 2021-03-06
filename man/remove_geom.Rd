% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geom_remove.R
\name{remove_geom}
\alias{remove_geom}
\title{Remove GEOM From \code{ggplot} Object}
\usage{
remove_geom(ggplot2_object, geom_type)
}
\arguments{
\item{ggplot2_object}{A \code{ggplot} object.}

\item{geom_type}{A string. The GEOM.}
}
\description{
This function removes a GEOM from a \code{ggplot} object,
 something that is useful in conjunction with automated
 \code{gg} tools like \code{ggtern}, \code{ggbiplot},
 or \code{ggvegan}.
}
\details{
Use like \code{remove_geom(g1, "GeomText")}.
}
\author{
Kamil Slowikowski from Stack Overflow.
}
