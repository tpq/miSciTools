#' Aptly Characterize Aptamers
#'
#' This function procedurally characterizes the secondary structure of RNA sequences
#'  using regular expressions. The characterizations provided may help identify
#'  potentially functional aptamers.
#'
#' This function requires a multi-sequence FASTA file as input. It also depends on the UNIX
#'  programs RNAfold and imagemagick. Install RNAfold through the ViennaRNA package.
#'
#' @param fasta A multi-sequence FASTA file.
#' @param RNAfold The commandline program to call. Append RNAfold arguments here.
#' @param cores The number of cores to use.
#' @param select A numeric or character vector. Selects which sequences from
#'  the multi-sequence FASTA file to include in the analysis.
#' @return A list of characteristics.
#'
#' @export
aptly <- function(fasta, RNAfold = "RNAfold", cores = 1, select){

  packageCheck(c("seqinr", "LncFinder"))

  message("* Reading FASTA...")
  Seqs <- seqinr::read.fasta(file = fasta)
  if(!missing(select)) Seqs <- Seqs[select]

  message("* Folding RNA...")
  SS.seq_2 <- LncFinder::run_RNAfold(Seqs, RNAfold.path = RNAfold, parallel.cores = cores)
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
    if(sum(find != -1) > 0){
      HairpinTipsTotal <- length(find) # number hairpin tips
      HairpinTipsLengthMean <- mean(attr(find, "match.length") - 2) # mean hairpin tip length
      HairpinTipsLengthTotal <- sum(attr(find, "match.length") - 2) # total hairpin tip
    }else{
      HairpinTipsTotal <- 0
      HairpinTipsLengthMean <- 0
      HairpinTipsLengthTotal <- 0
    }

    find <- gregexpr(paste0(lopen, lcont, pin, rcont, ropen), note)[[1]]
    if(sum(find != -1) > 0){
      HairpinsTotal <- length(find) # number hairpins
      HairpinsLengthMean <- mean(attr(find, "match.length") - 2) # mean hairpin length
      HairpinsLengthTotal <- sum(attr(find, "match.length") - 2) # total hairpin
    }else{
      HairpinsTotal <- 0
      HairpinsLengthMean <- 0
      HairpinsLengthTotal <- 0
    }

    find <- gregexpr(paste0(lopen, pin, lopen), note)[[1]]
    if(sum(find != -1) > 0){
      LeftBulgesTotal <- length(find) # number hairpins
      LeftBulgesLengthMean <- mean(attr(find, "match.length") - 2) # mean bulge length
      LeftBulgesLengthTotal <- sum(attr(find, "match.length") - 2) # total bulge
    }else{
      LeftBulgesTotal <- 0
      LeftBulgesLengthMean <- 0
      LeftBulgesLengthTotal <- 0
    }

    find <- gregexpr(paste0(ropen, pin, ropen), note)[[1]]
    if(sum(find != -1) > 0){
      RightBulgesTotal <- length(find) # number hairpins
      RightBulgesLengthMean <- mean(attr(find, "match.length") - 2) # mean bulge length
      RightBulgesLengthTotal <- sum(attr(find, "match.length") - 2) # total bulge
    }else{
      RightBulgesTotal <- 0
      RightBulgesLengthMean <- 0
      RightBulgesLengthTotal <- 0
    }

    find <- gregexpr(paste0(lonly, pin, lonly, pin, ronly, pin, ronly), note)[[1]]
    if(sum(find != -1) > 0){
      BulgedPinsTotal <- length(find) # number pre-hairpin bulges
      BulgedPinsLengthMean <- mean(attr(find, "match.length") - 2) # mean pre-hairpin bulge length
      BulgedPinsLengthTotal <- sum(attr(find, "match.length") - 2) # pre-hairpin bulges
    }else{
      BulgedPinsTotal <- 0
      BulgedPinsLengthMean <- 0
      BulgedPinsLengthTotal <- 0
    }

    data.frame(
      LeadLength, LagLength,
      HairpinTipsTotal, HairpinTipsLengthMean, HairpinTipsLengthTotal,
      HairpinsTotal, HairpinsLengthMean, HairpinsLengthTotal,
      LeftBulgesTotal, LeftBulgesLengthMean, LeftBulgesLengthTotal,
      RightBulgesTotal, RightBulgesLengthMean, RightBulgesLengthTotal,
      BulgedPinsTotal, BulgedPinsLengthMean, BulgedPinsLengthTotal
    )
  })

  table <- do.call("rbind", out)
  if(!is.null(names(SS.seq_2))) rownames(table) <- names(SS.seq_2)
  utils::write.csv(table, file = paste0(RNAfold, "-summary.csv"))
  return(table)

  # message("* Creating hybrid pseudo-sequences...")
  # out <- lapply(SS.seq_2, function(SEQobj){
  #
  #   s1 <- toupper(SEQobj[1])
  #   s1 <- gsub("U", "T", s1)
  #   s1 <- strsplit(s1, "")[[1]]
  #
  #   s2 <- gsub("\\(", ".", SEQobj[2])
  #   s2 <- gsub("\\)", "+", s2)
  #   s2 <- strsplit(s2, "")[[1]]
  #
  #   dots <- gregexpr("\\.", SEQobj[2])[[1]]
  #
  #   s2[dots] <- s1[dots]
  #   return(paste0(s2, collapse = ""))
  # })
  #
  # message("* Aligning pseudo-sequences to reference...")
  # a <- Biostrings::DNAStringSet(unlist(out))
  # a <- a[order(Biostrings::width(a), decreasing = TRUE)]
  # b <- Biostrings::pairwiseAlignment(a, a[1])
  # c <- Biostrings::BStringSet(b)
  # names(c) <- names(a)
  #
  # message("* Finding consensuses...")
  # pwm <- Biostrings::consensusMatrix(c)
  # pwm <- pwm[c("A", "T", "G", "C", ".", "+"), ]
  # rownames(pwm) <- c("A", "T", "G", "C", "(", ")")
  #
  # string <- Biostrings::consensusString(c)
  # string <- gsub("\\.", "(", string)
  # string <- gsub("\\+", ")", string)
  # string <- gsub("\\?", "N", string)
  # title <- paste0("Visualization of Consensus Sequence: [", string, "]")
  #
  # g <- ggplot2::ggplot(reshape2::melt(pwm), ggplot2::aes(x = Var2, y = value, fill = Var1)) +
  #   ggplot2::geom_bar(stat = "identity") + ggplot2::scale_fill_brewer(palette = "Set2") +
  #   ggplot2::xlab("Distance from Reference Origin") + ggplot2::ylab("Frequency of Base") +
  #   ggplot2::labs("fill" = "Base") + ggplot2::theme_bw() + ggplot2::ggtitle(title)
  # plot(g)
  #
  # return(
  #   list(
  #     "table" = table, "hybrid" = a, "aligned" = c,
  #     "g" = g, "pwm" = pwm, "string" = string
  #   )
  # )
}
