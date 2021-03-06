% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fig.R
\name{fig.fromAdj}
\alias{fig.fromAdj}
\title{Plot Adjacency Matrix}
\usage{
fig.fromAdj(A, n.group = 1)
}
\arguments{
\item{A}{An adjacency matrix.}

\item{n.group}{A character vector or factor. Groupings for the nodes.
Used to color the nodes in viridis hues.}
}
\value{
A \code{ggplot} object.
}
\description{
This function uses \code{ggraph} to make a \code{ggplot} object.
}
\examples{
A <- matrix(sample(0:1, 25, replace = TRUE), 5, 5)
fig.fromAdj(A, n.group = c(1, 2, 1, 2, 1))
}
