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
getRecount <- function(ID, characteristic = 1){

  catchwd <- getwd()

  miSciTools::packageCheck("recount")
  miSciTools::packageCheck("SummarizedExperiment")

  ## Download data summarized at gene-level
  url <- recount::download_study(ID)

  ## Build directory structure
  `%+%` <- function(x, y) paste0(x, y)
  system("mkdir recount2-" %+% ID)
  load(getwd() %+% "/" %+% ID %+% "/rse_gene.Rdata")
  setwd(getwd() %+% "/recount2-" %+% ID)
  system("mkdir 0-make")
  system("mkdir annot")
  system("mkdir data")

  ## Find per-disease label
  gc()
  colData <- SummarizedExperiment::colData
  assays <- SummarizedExperiment::assays
  colnames(colData(rse_gene))
  annot <- colData(rse_gene)
  projects <- do.call("rbind", lapply(annot$characteristics, t))
  projects <- projects[, characteristic]
  project <- table(projects)
  if(!identical(sum(project), nrow(colData(rse_gene)))) stop()
  prj <- names(project)

  ## Split data by per-disease label
  for(i in 1:length(prj)){

    # Get counts and annot
    rse.i <- rse_gene[, projects == prj[i]]
    rse.i <- recount::scale_counts(rse.i)
    scale_counts <- assays(rse.i)[[1]]
    annot <- colData(rse.i)
    annot.split <- do.call("rbind", lapply(annot$characteristics, t))
    annot <- cbind(annot, annot.split)
    annot$characteristics <- "See V1-VN"

    # Remove NA columns and commas
    allisNA <- apply(annot, 2, function(x) all(is.na(x)))
    annot <- annot[, !allisNA]
    for(j in 1:ncol(annot)){
      if(any(grepl(",", annot[,j]))){
        annot[,j] <- gsub(",", "-", annot[,j])
      }
    }

    if(!identical(rownames(annot),
                  colnames(scale_counts))) stop("")

    # Save counts and annot
    id <- gsub(" ", "_", prj[i])
    utils::write.csv(scale_counts, file = paste0("data/", id, "-scale_counts.csv"))
    utils::write.csv(annot, file = paste0("annot/", id, "-annot.csv"))
  }

  # Add miSciTools note
  utils::write.table("see miSciTools::getRecount", "0-make/0-make.R")

  setwd(catchwd)
}
