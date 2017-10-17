#' Sequence Dissimilarity by Secondary Structure
#'
#' This function calculates the dissimilarity between sequences based on secondary structure.
#'  If \code{estimateConvergence = TRUE}, this function also calculates the dissimilarity
#'  based on the primary sequence and saves a convergence measure as output. If
#'  \code{groupBy > 0}, this function also clusters sequences by secondary structure and
#'  saves additional files as output (including secondary structure figures).
#'
#' This function requires a multi-sequence FASTA file as input. It also depends on the UNIX
#'  programs RNAfold and imagemagick. Install RNAfold through the ViennaRNA package.
#'
#' @param fasta A multi-sequence FASTA file.
#' @param rmTails Toggles whether to remove trailing "..." from start and end
#'  of the secondary structure bracket notation.
#' @param estimateConvergence Toggles whether to estimate secondary convergence
#'  (log of the secondary alignment score divided by primary alignment score).
#'  Produces intermediate files.
#' @param gapOpening,gapExtension Penalty for alignment gaps.
#'  Passed to \code{Biostrings::pairwiseAlignment}.
#' @param groupBy The height of the tree at which to cut when defining clusters.
#'  If \code{groupBy >= 1}, defines the total number of clusters to extract.
#'  Disable by setting \code{groupBy = 0}.
#'  Produces intermediate files.
#' @param cores The number of cores to use.
#'
#' @return A dissimilarity measure from [0, 1].
#'
#' @importFrom foreach %:% %dopar%
#' @export
fastafolder <- function(fasta, rmTails = FALSE, estimateConvergence = FALSE,
                        gapOpening = 10, gapExtension = 4,
                        groupBy = .5, cores = 1){

  packageCheck(
    c("magick", "seqinr", "LncFinder", "Biostrings",
      "parallel", "doParallel", "foreach")
  )

  if(!groupBy == 0) packageCheck(c("GeneR", "GeneRfold"))

  seqplot <- function(seq, file){
    temp <- tempfile()
    GeneRfold::rnaPlot(seq, file = temp)
    ps <- magick::image_read(temp)
    if(!missing(file)){
      magick::image_write(ps, file)
    }else{
      graphics::plot(ps)
    }
  }

  message("* Reading FASTA...")
  Seqs <- seqinr::read.fasta(file = fasta)

  message("* Folding RNA...")
  SS.seq_2 <- LncFinder::run_RNAfold(Seqs, RNAfold.path = "RNAfold", parallel.cores = cores)
  system("rm rna.ps")

  if(rmTails){
    message("* Removing unbound tails...")
    for(i in 1:length(SS.seq_2)){
      SS.seq_2[[i]][2] <- gsub("^\\.*", "", SS.seq_2[[i]][2], perl = TRUE)
      SS.seq_2[[i]][2] <- gsub("\\.*$", "", SS.seq_2[[i]][2], perl = TRUE)
    }
  }

  message("* Aligning secondary structures...")
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  out2 <- foreach::foreach(i = 1:length(SS.seq_2), .combine = cbind) %:%
    foreach::foreach(j = 1:length(SS.seq_2), .combine = c) %dopar% {
      if(i > j){ ij <- 0
      }else{
        ij <- Biostrings::pairwiseAlignment(
          pattern = SS.seq_2[[i]][2], subject = SS.seq_2[[j]][2],
          gapOpening = gapOpening, gapExtension = gapExtension,
          type = "local", scoreOnly = TRUE)
      }}
  colnames(out2) <- names(SS.seq_2)
  rownames(out2) <- names(SS.seq_2)
  parallel::stopCluster(cl)

  message("* Calculating distance...")
  a <- do.call("rbind", rep(list(diag(out2)), nrow(out2)))
  b <- do.call("cbind", rep(list(diag(out2)), nrow(out2)))
  sim2 <- out2 / pmax(a, b)
  sim2[upper.tri(sim2)] <- sim2[lower.tri(sim2)]
  dis <- 1 - abs(sim2)

  system("mkdir fastafolder")
  utils::write.csv(dis, "fastafolder/dis-2ary.csv")

  if(estimateConvergence){

    message("* Aligning primary sequences...")
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    out1 <- foreach::foreach(i = 1:length(SS.seq_2), .combine = cbind) %:%
      foreach::foreach(j = 1:length(SS.seq_2), .combine = c) %dopar% {
        if(i > j){ ij <- 0
        }else{
          ij <- Biostrings::pairwiseAlignment(
            pattern = SS.seq_2[[i]][1], subject = SS.seq_2[[j]][1],
            gapOpening = gapOpening, gapExtension = gapExtension,
            type = "local", scoreOnly = TRUE)
        }}
    colnames(out1) <- names(SS.seq_2)
    rownames(out1) <- names(SS.seq_2)
    parallel::stopCluster(cl)

    message("* Estimating convergence...")
    a <- do.call("rbind", rep(list(diag(out1)), nrow(out1)))
    b <- do.call("cbind", rep(list(diag(out1)), nrow(out1)))
    sim1 <- out1 / pmax(a, b)
    sim1[upper.tri(sim1)] <- sim1[lower.tri(sim1)]
    utils::write.csv(1 - abs(sim1), "fastafolder/dis-1ary.csv")

    converge <- log(sim2/sim1)
    utils::write.csv(converge, "fastafolder/convergence.csv")
  }

  if(!groupBy == 0){

    message("* Clustering secondary structures...")
    hc <- stats::hclust(stats::as.dist(dis))
    if(groupBy >= 1){ c <- stats::cutree(hc, k = groupBy)
    }else{ c <- stats::cutree(hc, h = groupBy) }
    utils::write.csv(table(c), "fastafolder/cutree.csv")

    message("* Making figures...")
    for(i in 1:length(c)){
      if(!dir.exists(paste0("fastafolder/", c[i]))) system(paste0("mkdir fastafolder/", c[i]))
      seqplot(SS.seq_2[[i]][1], file = paste0("fastafolder/", c[i], "/", names(c)[i], ".png"))
    }
  }

  return(dis)
}

#' Aptly Characterize Aptamers
#'
#' This function procedurally characterizes the secondary structure of RNA sequences
#'  using regular expressions. The characterizations provided may help identify
#'  potentially functional aptamers.
#'
#' This function requires a multi-sequence FASTA file as input. It also depends on the UNIX
#'  programs RNAfold and imagemagick. Install RNAfold through the ViennaRNA package.
#'
#' @inheritParams fastafolder
#' @return A data.frame of characteristics.
#'
#' @export
aptly <- function(fasta, cores = 1){

  packageCheck(c("seqinr", "LncFinder", "Biostrings", "ggplot2", "reshape2"))

  message("* Reading FASTA...")
  Seqs <- seqinr::read.fasta(file = fasta)

  message("* Folding RNA...")
  SS.seq_2 <- LncFinder::run_RNAfold(Seqs, RNAfold.path = "RNAfold", parallel.cores = cores)
  system("rm rna.ps")

  # Define regex terms
  lopen <- "\\("
  lonly <- "\\(+"
  lcont <- "(\\(|\\.)+"
  pin <- "\\.+"
  ropen <- "\\)"
  ronly <- "\\)+"
  rcont <- "(\\)|\\.)+"

  message("* Characterizing structures...")
  out <- lapply(SS.seq_2, function(SEQobj){

    note <- SEQobj[2]

    find <- gregexpr("^\\.*", note, perl = TRUE)[[1]]
    LeadLength <- attr(find, "match.length")[1]

    find <- gregexpr("\\.*$", note, perl = TRUE)[[1]]
    LagLength <- attr(find, "match.length")[1]

    find <- gregexpr(paste0(lopen, pin, ropen), note)[[1]]
    HairpinTipsTotal <- length(find) # number hairpin tips
    HairpinTipsLength <- mean(attr(find, "match.length") - 2) # mean hairpin tip length
    HairpinTipsPercent <- sum(attr(find, "match.length") - 2) / nchar(note) # % hairpin tip

    find <- gregexpr(paste0(lopen, lcont, pin, rcont, ropen), note)[[1]]
    HairpinsTotal <- length(find) # number hairpins
    HairpinsLength <- mean(attr(find, "match.length") - 2) # mean hairpin length
    HairpinsPercent <- sum(attr(find, "match.length") - 2) / nchar(note) # % hairpin

    find <- gregexpr(paste0(lonly, pin, lonly, pin, ronly, pin, ronly), note)[[1]]
    BulgedPinsTotal <- length(find) # number pre-hairpin bulges
    BulgedPinsLength <- mean(attr(find, "match.length") - 2) # mean pre-hairpin bulge length
    BulgedPinsPercent <- sum(attr(find, "match.length") - 2) / nchar(note) # % pre-hairpin bulges

    data.frame(
      LeadLength, LagLength,
      HairpinTipsTotal, HairpinTipsLength, HairpinTipsPercent,
      HairpinsTotal, HairpinsLength, HairpinsPercent,
      BulgedPinsTotal, BulgedPinsLength, BulgedPinsPercent
    )
  })

  table <- do.call("rbind", out)
  if(!is.null(names(SS.seq_2))) rownames(table) <- names(SS.seq_2)

  message("* Creating hybrid pseudo-sequences...")
  out <- lapply(SS.seq_2, function(SEQobj){

    s1 <- toupper(SEQobj[1])
    s1 <- gsub("U", "T", s1)
    s1 <- strsplit(s1, "")[[1]]

    s2 <- gsub("\\(", ".", SEQobj[2])
    s2 <- gsub("\\)", "+", s2)
    s2 <- strsplit(s2, "")[[1]]

    dots <- gregexpr("\\.", SEQobj[2])[[1]]

    s2[dots] <- s1[dots]
    return(paste0(s2, collapse = ""))
  })

  message("* Aligning pseudo-sequences to reference...")
  a <- Biostrings::DNAStringSet(unlist(out))
  a <- a[order(Biostrings::width(a), decreasing = TRUE)]
  b <- Biostrings::pairwiseAlignment(a[2:length(a)], a[1])
  c <- Biostrings::BStringSet(b)

  message("* Finding consensuses...")
  pwm <- Biostrings::consensusMatrix(c)
  pwm <- pwm[c("A", "T", "G", "C", ".", "+"), ]
  rownames(pwm) <- c("A", "T", "G", "C", "(", ")")

  string <- Biostrings::consensusString(c)
  string <- gsub("\\.", "(", string)
  string <- gsub("\\+", ")", string)
  string <- gsub("\\?", "N", string)
  title <- paste0("Visualization of Consensus Sequence: [", string, "]")

  ggplot2::ggplot(reshape2::melt(pwm), ggplot2::aes(x = Var2, y = value, fill = Var1)) +
    ggplot2::geom_bar(stat = "identity") + ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::xlab("Distance from Reference Origin") + ggplot2::ylab("Frequency of Base") +
    ggplot2::labs("fill" = "Base") + ggplot2::theme_bw() + ggplot2::ggtitle(title)

  return(list(table, pwm, string))
}
