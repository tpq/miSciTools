#' @importFrom graphics plot
NULL

#' Paste with Sum
#'
#' Use %+% to paste0 two strings.
#'
#' @param a,b A string.
#' @return A string.
#' @examples
#' "Now you can paste() like " %+% " in Python."
#' @export
`%+%` <- function(a, b) paste0(a, b)

#' Save Figure
#'
#' This function saves a figure.
#'
#' @param plot The figure.
#' @param file A string. The file to save.
#' @param dir A string. The directory.
#' @param width,height The figure dimensions.
#' @param res The figure resolution.
#' @examples
#' \dontrun{
#' fig.save(plot(1:5, 1:5), file = "Figure1")
#' }
#' @export
fig.save <- function(plot, file = "fig", dir = getwd(), width = 5, height = 5, res = 600){

  oldwd <- getwd()
  setwd(dir)
  grDevices::png(file %+% ".png", width = width, height = height,
                 res = res, units = "in")
  plot
  grDevices::dev.off()
  setwd(oldwd)
}

#' Vertical Join Figures
#'
#' @param plot1,plot2 The figures.
#' @inheritParams fig.save
#' @examples
#' \dontrun{
#' fig.cbind(plot(1:5, 1:5), plot(1:5, 5:1), file = "Figure2")
#' }
#' @export
fig.cbind <- function(plot1, plot2, file = "fig", dir = getwd(), width = 5, height = 5, res = 600){

  oldwd <- getwd()
  setwd(dir)
  fig.save(plot1, file = "out9000-temp1" %+% file, width = width, height = height, res = res)
  fig.save(plot2, file = "out9000-temp2" %+% file, width = width, height = height, res = res)
  system("convert +append out9000-temp1* out9000-temp2* " %+% file %+% ".png")
  system("rm out9000-temp*")
  setwd(oldwd)
}

#' Horizontal Join Figures
#'
#' @param plot1,plot2 The figures.
#' @inheritParams fig.save
#' @examples
#' \dontrun{
#' fig.rbind(plot(1:5, 1:5), plot(1:5, 5:1), file = "Figure3")
#' }
#' @export
fig.rbind <- function(plot1, plot2, file = "fig", dir = getwd(), width = 5, height = 5, res = 600){

  oldwd <- getwd()
  setwd(dir)
  fig.save(plot1, file = "out9000-temp1" %+% file, width = width, height = height, res = res)
  fig.save(plot2, file = "out9000-temp2" %+% file, width = width, height = height, res = res)
  system("convert -append out9000-temp1* out9000-temp2* " %+% file %+% ".png")
  system("rm out9000-temp*")
  setwd(oldwd)
}

#' Tile Join Figures
#'
#' @param ... The figures.
#' @param tile A string. Describes how to tile the figures (e.g., "1x3").
#' @inheritParams fig.save
#' @examples
#' \dontrun{
#' fig.tile(plot(1:5), plot(6:10), plot(11:15), tile = "1x3", file = "Figure 4")
#' }
#' @export
fig.tile <- function(..., tile = "2x2", file = "fig", dir = getwd(), width = 5, height = 5, res = 600){

  oldwd <- getwd()
  setwd(dir)
  args <- as.list(substitute(list(...)))[-1]
  for(i in 1:length(args)){
    fig.save(eval(args[[i]]), file = "out9000-temp" %+% i %+% file,
             width = width, height = height, res = res)
  }
  system("montage out9000-temp* -tile " %+% tile %+% " -mode concatenate " %+% file %+% ".png")
  system("rm out9000-temp*")
  setwd(oldwd)
}

#' Make Figure from Table
#'
#' @param mat The matrix.
#' @param colmat A matrix of colors. Specifies the colors
#'  for each cell in the matrix.
#' @return A grob.
#' @examples
#' colmat <- matrix("red", 5, 5)
#' colmat[,1] <- "blue"
#' fig.fromTable(iris[1:5, 1:5], colmat)
#' @export
fig.fromTable <- function(mat, colmat = matrix("black", nrow(mat), ncol(mat))){

  if(!identical(dim(mat), dim(colmat))){
    stop("'mat' and 'colmat' must have same dimensions.")
  }

  table_theme <- gridExtra::ttheme_minimal(
    core = list(fg_params = list(col = colmat)))
  grid::grid.newpage()
  fig <- gridExtra::tableGrob(mat, theme = table_theme)
  plot(fig)
  fig
}

#' Label Figure Panel
#'
#' This function adds a label to the top-left of the figure.
#'
#' @param plot The figure.
#' @param label The label.
#' @return A grob.
#' @examples
#' grob <- fig.fromTable(iris[1:5,1:5])
#' fig.label(grob, "1A) An example")
#' @export
fig.label <- function(plot, label = "1) Place Label Here"){

  tryCatch(
    gridExtra::grid.arrange(
      gridExtra::arrangeGrob(
        plot, top = grid::textGrob(
          label,
          x = grid::unit(0.05, "npc"),
          y = grid::unit(0, "npc"),
          just = c("left","top"),
          gp = grid::gpar(col = "black",
                          fontsize = 18))
      )
    ),
    error = function(e){
      stop("This method does not work with base R graphics.")
    }
  )
}

#' Plot Adjacency Matrix
#'
#' This function uses \code{ggraph} to make a \code{ggplot} object.
#'
#' @param A An adjacency matrix.
#' @param n.group A character vector or factor. Groupings for the nodes.
#'  Used to color the nodes in viridis hues.
#' @return A \code{ggplot} object.
#' @examples
#' A <- matrix(sample(0:1, 25, replace = TRUE), 5, 5)
#' fig.fromAdj(A, n.group = c(1, 2, 1, 2, 1))
#' @export
fig.fromAdj <- function(A, n.group = 1){

  packageCheck("igraph")
  packageCheck("ggraph")

  # Turn adjacency matrix into igraph object
  graph <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  igraph::V(graph)$class <- as.character(igraph::V(graph))
  igraph::V(graph)$grp <- rep_len(n.group, nrow(A))

  # Get x-lim and y-lim (used later to make space for panel legend)
  set.seed(1)
  lay <- ggraph::create_layout(graph, "igraph", algorithm = "nicely")
  Ymax <- max(lay$y)
  Ymin <- min(lay$y)

  # Plot ggraph object from igraph object
  set.seed(1)
  ggplot <-
    ggraph::ggraph(graph) +
    ggraph::geom_edge_link(colour = viridis::viridis(2)[1], alpha = .2) +
    ggraph::geom_node_point(ggplot2::aes_string(colour = "grp"), size = 10, alpha = .1) +
    ggraph::geom_node_text(ggplot2::aes_string(label = "class")) +
    ggplot2::ylim(Ymin, Ymax + abs(Ymax * .05)) + # increases a negative or positive Ymax
    ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      legend.position = "none",
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )

  ggplot <- tryCatch(
    plot(ggplot + viridis::scale_color_viridis(discrete = FALSE)),
    error = function(e) plot(ggplot + viridis::scale_color_viridis(discrete = TRUE))
  )
}

#' Make \code{ggplot} Pretty
#'
#' @param ggplot A \code{ggplot} object.
#' @param col A boolean. Toggles whether to color the figure with viridis.
#' @param rotate.x A boolean. Toggles whether to rotate the x-axis.
#' @examples
#' \dontrun{
#' p <- ggplot(mtcars, aes(wt, mpg, col = cyl))
#' p <- p + geom_point()
#' fig.pretty(p, col = TRUE)
#' }
#' @export
fig.pretty <- function(ggplot, col = FALSE, rotate.x = FALSE){

  ggplot <- ggplot +
    ggplot2::theme(text = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(size = 24)) +
    ggplot2::theme_bw()

  if(col) ggplot <- tryCatch(
    plot(ggplot + viridis::scale_color_viridis(discrete = FALSE)),
    error = function(e) plot(ggplot + viridis::scale_color_viridis(discrete = TRUE))
  )

  if(rotate.x) ggplot <- ggplot +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  return(ggplot)
}

#' Make 3D Plot
#'
#' @param data A \code{data.frame}.
#' @param x,y,z A numeric. The x, y, and z dimensions to plot.
#' @param n.group A character vector or factor. Groupings for the points.
#'  Used to color the figure.
#' @param main A string. The plot title.
#' @examples
#' \dontrun{
#' fig.quick3D(mtcars, x = 1, y = 2, z = 3, n.group = mtcars$cyl)
#' }
#' @export
fig.3D <- function(data, x = 1, y = 2, z = 3, n.group = 1, main = ""){

  packageCheck("lattice")

  if(is.null(colnames(data))){
    stop("Provided data must have column names.")
  }

  data <- as.data.frame(data)
  xlab <- colnames(data)[x]
  ylab <- colnames(data)[y]
  zlab <- colnames(data)[z]
  colnames(data)[c(x, y, z)] <- c("x", "y", "z")
  grp <- rep_len(n.group, nrow(data))

  lattice::cloud(z ~ x * y, data = data,
                 group = grp, main = main, pch = 3, alpha = .8,
                 xlab = list(xlab, rot = 29),
                 ylab = list(ylab, rot = 321),
                 zlab = list(zlab, rot = 94+180))
}
