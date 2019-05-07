#' Check for Packages
#'
#' Triggers an error if a package is not installed.
#'
#' @param packages A character vector of packages.
#'
#' @export
packageCheck <- function(packages){

  for(package in packages){

    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Uh oh! This method depends on ", package,
           ". ", "Please install.")
    }
  }

  return(TRUE)
}

#' Plot Multiple Graphs
#'
#' Easily plot multiple graphs within the same window. Code adapted from
#'  http://www.cookbook-r.com/.
#'
#' @param ... Multiple plots.
#' @param cols A numeric scalar. The number of plot columns.
#'
#' @export
multiplot <- function(..., cols = 1){

  if(!requireNamespace("grid", quietly = TRUE)){
    stop("Uh oh! This method depends on grid. ",
         "Try running: miSciTools::demand('grid')")
  }

  # Make a list of plots
  plots <- list(...)
  numPlots <- length(plots)

  # Layout the panel
  layout <- matrix(seq(from = 1,
                       to = cols * ceiling(numPlots/cols)),
                   ncol = cols,
                   nrow = ceiling(numPlots/cols)
  )

  if(numPlots == 1){

    print(plots[[1]])

  }else{

    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow(layout),
                                   ncol(layout))
      )
    )

    # Place each plot, in the correct location
    for(i in 1:numPlots){

      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(
        plots[[i]],
        vp = grid::viewport(layout.pos.row = matchidx$row,
                            layout.pos.col = matchidx$col)
      )
    }
  }

  return(TRUE)
}

#' Load or Install Package
#'
#' This function installs a package if not already installed, then loads it.
#'
#' @param packages A character vector. The package(s) to load or install.
#'
#' @export
demand <- function(packages){

  for(package in packages){

    if(!is.element(package, utils::installed.packages()[, 1])){

      try({

        cat("Looking for", package, "at CRAN...\n")
        suppressWarnings(
          utils::install.packages(package, repos = "http://cran.us.r-project.org")
        )
      })
    }

    if(!is.element(package, utils::installed.packages()[, 1])){

      cat("Looking for", package, "at Bioconductor...\n")
      source("https://bioconductor.org/biocLite.R")
      suppressWarnings(
        BiocManager::install(package, update = FALSE)
      )
    }

    if(!is.element(package, utils::installed.packages()[, 1])){

      try({

        cat("Looking for", package, "at R-Forge...\n")
        suppressWarnings(
          utils::install.packages(package, repos = "http://R-Forge.R-project.org")
        )
      })
    }

    if(is.element(package, utils::installed.packages()[,1])){

      library(package, character.only = TRUE)

    }else{

      warning("Could not find ", package, "!")
    }
  }
}

#' Set Temporary Object
#'
#' Set a given value to a temporary object in a specified environment.
#'
#' @param value An object to assign as a temporary object.
#' @param pos Where to do the assignment (relative to the \code{asTempObj}
#'  environment). Defaults to \code{pos = 1}.
#' @param envir The environment to use.
#' @return The name of the temporary object as assigned in the provided
#'  environment.
#'
#' @export
asTempObj <- function(value, pos = 1, envir = as.environment(pos)){

  alreadyInUse <- TRUE
  while(alreadyInUse){

    id <- paste(sample(c(letters, LETTERS), 8), collapse = "")
    alreadyInUse <- id %in% ls(pos = pos, envir = envir)
  }

  assign(id, value, pos = pos, envir = envir)
  return(id)
}

#' Write R Script
#'
#' This function saves a character string as an R script. By default, this
#'  script gets stored in a temporary directory as a temporarily file.
#'  This function returns the file path for the saved R script for later
#'  use (e.g., via \code{\link{qsub}}).
#'
#' @param ... Any number of character strings or R expressions to join
#'  together and save as an R script.
#' @param file A character string. The file path, including directory, where
#'  to save the new script. Defaults to a temporary file name.
#' @param preview A logical scalar. Toggles whether to preview the script
#'  in the console before saving it. Defaults to \code{FALSE}.
#' @return The file path for the saved R script.
#'
#' @export
writeR <- function(..., file = paste0(tempfile(), ".R"), preview = FALSE){

  # Save parent environment to temporary directory
  file.wd <- paste0(file, "Data")
  save.image(file = file.wd)

  # Combine strings and expressions into an R script
  R <- paste0(..., collapse = "")

  # Load parent environment from within R script
  R <- paste0("load(\"", file.wd, "\")\n", R)

  # Save R script in a temporary directory
  if(preview) cat(R)
  file.create(file)
  fileConn <- file(file)
  writeLines(R, fileConn)
  close(fileConn)
  return(file)
}

#' Qsub Linux Command
#'
#' This function sends a Linux command to the PBS queue via \code{qsub}.
#'
#' @param cmd A character string. A Linux command to \code{qsub} or the
#'  location of an R script to \code{qsub}.
#' @param ... Any additional PBS argument(s). Each argument should
#'  get named according to the bash character argument (e.g., set
#'  \code{-M thom@tpq.me} with \code{qsub(command, M = "thom@tpq.me")}).
#'
#' @export
qsub <- function(cmd, ...){

  # Make sure 'cmd' is character input
  if(class(cmd) != "character"){

    stop("Please supply the 'cmd' argument as a character string.")
  }

  # Check if 'cmd' is an R file
  if(file.exists(cmd) & substr(cmd, nchar(cmd)-1, nchar(cmd)) == ".R"){

    cmd <- paste0("R CMD BATCH ", cmd)
  }

  # Prepare #PBS args
  args <- as.list(substitute(list(...)))[-1]
  if(length(args) > 0){

    args.bash <- paste0("#PBS -", names(args), " ", unlist(args), "\n",
                        collapse = "")

  }else{

    args.bash <- "#PBS\n"
  }

  # Write script to execute cmd
  bash <- paste0("#! /bin/bash", "\n",
                 args.bash, "\n", cmd)

  # Send script to queue
  run <- paste0("echo \"", bash, "\" | qsub")
  system(run)

  return(TRUE)
}

#' Insert New Line
#'
#' Inserts a new line of text into an existing file at a specified breakpoint.
#'
#' @param file A character string. The text file to modify.
#' @param after A character string. A fragment of existing text that marks
#'  the line after which to insert the new text.
#' @param what A character string. The new text to insert.
#'
#' @export
insert <- function(file, after, what){

  # Read lines from connection
  con <- file(file)
  lines <- readLines(con)

  # Find where to insert the new line
  i <- which(grepl(after, lines))[1]

  # Insert the new line
  if(i == length(lines)){

    new <- c(lines[1:i], as.character(what))

  }else{
    new <- c(lines[1:i], as.character(what), lines[(i+1):length(lines)])
  }

  # Write lines to connection
  writeLines(new, con = con)
  close(con)
}

#' Make Progress Bar
#'
#' @param i The current iteration.
#' @param k Total iterations.
#' @param numTicks The result of \code{progress}.
#' @return The next \code{numTicks} argument.
#'
#' @export
progress <- function(i, k, numTicks){

  if(i == 1) numTicks <- 0

  if(numTicks == 0) cat("|-")

  while(i > numTicks*(k/40)){

    cat("-")
    if(numTicks == 10) cat("(25%)")
    if(numTicks == 20) cat("(50%)")
    if(numTicks == 30) cat("(75%)")
    numTicks <- numTicks + 1
  }

  if(i == k) cat("-|\n")

  return(numTicks)
}
