#' Get recount2 Data
#'
#' This function retrieves data from the recount2 server.
#'
#' @param ID A character string. The recount2 ID.
#' @param characteristic An integer. The characteristic
#'  used to split the data.
#'
#' @return Saves data to working directory.
#'
#' @export
getRecount <- function(ID, characteristic){

  catchwd <- getwd()
  packageCheck("recount")
  packageCheck("SummarizedExperiment")
  url <- recount::download_study(ID)
  `%+%` <- function(x, y) paste0(x, y)

  # Make folder structure
  system("mkdir recount2-" %+% ID)
  load(getwd() %+% "/" %+% ID %+% "/rse_gene.Rdata")
  setwd(getwd() %+% "/recount2-" %+% ID)
  system("mkdir 0-make")
  system("mkdir annot")
  system("mkdir data")
  gc()

  # Retrieve data from RSE object
  colData <- SummarizedExperiment::colData
  assays <- SummarizedExperiment::assays
  colnames(colData(rse_gene))
  annot <- colData(rse_gene)

  # Get V1-VN characteristics
  projSummary <- function(data){
    annots <- lapply(data, t)
    annots <- lapply(annots, as.data.frame)
    do.call(plyr::rbind.fill, annots)
  }

  # Ask user which characteristics to use
  projects <- projSummary(annot$characteristics)
  if(missing(characteristic)){
    # Ask for user input
    lapply(projects, function(x){
      print(table(x))
    })
    characteristic <- readline(paste(
      "Which characteristic do you want to use to split files?",
      "[Input Numeric]: "))
    characteristic <- as.numeric(characteristic)
  }

  projects <- projects[, characteristic]
  project <- table(projects)
  if (!identical(sum(project), nrow(colData(rse_gene))))
    stop()
  prj <- names(project)
  for (i in 1:length(prj)) {

    # Subset the i-th project, scale counts
    rse.i <- rse_gene[, projects == prj[i]]
    rse.i <- recount::scale_counts(rse.i)
    scale_counts <- assays(rse.i)[[1]]
    annot <- colData(rse.i)

    # Summarize projects for annot.i
    annot.split <- projSummary(annot$characteristics)
    annot <- cbind(annot, annot.split)
    annot$characteristics <- "See V1-VN"

    # Remove all NA annot
    allisNA <- apply(annot, 2, function(x) all(is.na(x)))
    annot <- annot[, !allisNA]
    for (j in 1:ncol(annot)) {
      if (any(grepl(",", annot[, j]))) {
        annot[, j] <- gsub(",", "-", annot[, j])
      }
    }

    # Check counts.i and annot.i are same
    if (!identical(rownames(annot), colnames(scale_counts)))
      stop("")

    # Save counts.i and annot.i
    id <- gsub(" ", "_", prj[i])
    utils::write.csv(scale_counts, file = paste0("data/",
                                                 id, "-scale_counts.csv"))
    utils::write.csv(annot, file = paste0("annot/", id, "-annot.csv"))
  }

  # Leave miSciTools note
  utils::write.table(paste("see miSciTools::getRecount v.",
                           utils::packageVersion("miSciTools")),
                     file = "0-make/0-make.R")
  setwd(catchwd)
}
