###########################################################
### Functions to assist in proportionality analyses

# Calculate phi from a data.frame of feature counts
phit <- function(counts, symmetrize = TRUE){
  
  require(compositions)
  
  # Centered log-ratio transform count data matrix
  counts.clr <- as.data.frame(clr(t(counts)))
  
  # Replace zeroes with next smallest number
  counts[counts == 0] <- unique(sort(as.matrix(counts)))[2]
  
  # Calculate the variance of the log-ratio
  counts.vlr <- variation(acomp(t(counts)))
  colnames(counts.vlr) <- rownames(counts)
  rownames(counts.vlr) <- rownames(counts)
  
  # Calculate var of clr transformed treatments
  counts.clr.var <- apply(counts.clr, 2, var)
  
  # Sweep out clr transformed variance from vlr
  counts.phi <- sweep(counts.vlr, 2, counts.clr.var, FUN = "/")
  
  # Symmetrize matrix if symmetrize = TRUE
  if(symmetrize) counts.phi <- phitSym(counts.phi)
  
  return(counts.phi)
}

# Calculate p(phi) based on a NULL distribution
phitDistr <- function(counts, iter = 10, iterSize = nrow(counts), returnPval = TRUE){
  
  distr <- vector("numeric", iter * iterSize * (iterSize - 1) / 2)
  
  for(i in 1:iter){
    
    cat(paste0("Calculating all phi for iter ", i, "...\n"))
    null.i <- apply(counts, 2, sample, iterSize)
    phi.i <- phit(null.i)
    begin <- (i - 1) * iterSize * (iterSize - 1) / 2 + 1
    end <- (i - 1) * iterSize * (iterSize - 1) / 2 + iterSize * (iterSize - 1) / 2
    distr[begin:end] <- phitTri(phi.i)
    rm(phi.i)
  }
  
  cat("Fitting phi to distribution...\n")
  fit <- ecdf(distr)
  rm(distr)
  
  if(returnPval){
    
    cat("Calculating all phi for actual counts...\n")
    phi <- phit(counts)
    raw <- phitRaw(phi)
    
    cat("Using 'fit' to convert phi into pval...\n")
    pval <- fit(raw$phi)
    cat("Correcting for multiple testing...\n")
    fdr <- p.adjust(pval, method = "BH")
    
    cat("Building results...\n")
    result <- data.frame(raw, "pval" = pval, "fdr.BH" = fdr, stringsAsFactors = FALSE)
    final <- result[order(result$pval),]
    rownames(final) <- 1:nrow(final)
    
    return(final)
    
  }else{
    
    return(fit)
  }
}

# Retrieve phi for each feature pair
phitRaw <- function(phi){
  
  index.i <- vector("numeric", length = (nrow(phi) - 1)*nrow(phi)/2)
  index.j <- vector("numeric", length = (nrow(phi) - 1)*nrow(phi)/2)
  index.phi <- vector("numeric", length = (nrow(phi) - 1)*nrow(phi)/2)
  counter <- 1
  
  for(j in 2:nrow(phi)){
    
    for(i in 1:(j-1)){
      
      index.i[counter] <- i
      index.j[counter] <- j
      index.phi[counter] <- phi[j, i]
      counter <- counter + 1
    }
  }
  
  result <- data.frame("Feature.1" = rownames(phi)[index.i],
                       "Feature.1.index" = index.i,
                       "Feature.2" = rownames(phi)[index.j],
                       "Feature.2.index" = index.j,
                       "phi" = index.phi,
                       stringsAsFactors = FALSE)
  
  final <- result[order(result$phi),]
  rownames(final) <- 1:nrow(final)
  
  return(final)
}

# Retrieve the lower triangle of a phi matrix
phitTri <- function(phi){
  
  result <- vector("numeric", length = (nrow(phi) - 1)*nrow(phi)/2)
  counter <- 1
  
  for(j in 2:nrow(phi)){
    
    for(i in 1:(j-1)){
      
      result[counter] <- phi[j, i]
      counter <- counter + 1
    }
  }
  
  return(result)
}

# Symmetrize a phi matrix
phitSym <- function(phi){
  
  for(j in 2:nrow(phi)){
    
    for(i in 1:(j-1)){
      
      phi[i, j] <- phi[j, i]
    }
  }
  
  return(phi)
}

###########################################################
### Functions to assist in genomic annotations

# 'Lifts over' coordinates from one genome build to another
# NOTE: Arguments 'from' and 'to' must refer to valid UCSC builds
# NOTE: Non-UCSC coordinates converted to UCSC equivalents
lift <- function(object, from = "hg18", to = "hg19", flatGrl = TRUE){
  
  demand("rtracklayer")
  demand("biovizBase")
  demand("R.utils")
  
  warning("Object names will get lost during 'lift over'. Consider storing names as column.\n")
  
  if(!all(genome(object) %in% from)){
    
    warning("Check 'genome(object)' to make sure you selected correct 'over.chain' directory.")
  }
  
  over.chain <- paste0(from, "To", capitalize(to), ".over.chain.gz")
  cat("Downloading", over.chain, "as temporary file...\n")
  url <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/", from, "/liftOver/", over.chain)
  file.gz <- tempfile(fileext = ".over.chain.gz")
  download.file(url, file.gz)
  
  cat("Unzipping", over.chain, "temporary file...\n")
  file.oc <- tempfile(fileext = ".over.chain")
  gunzip(filename = file.gz, destname = file.oc)
  
  cat("Importing unzipped file...\n")
  chain <- import.chain(file.oc)
  
  cat("Convering coordinates to UCSC equivalents...\n")
  tempStyle <- seqlevelsStyle(object)
  seqlevelsStyle(object) <- "UCSC"
  
  cat("Using over.chain to 'lift over'...\n")
  lifted <- liftOver(object, chain)
  genome(lifted) <- to
  
  cat("Returning coordinates to original style...\n")
  seqlevelsStyle(object) <- tempStyle
  
  if(flatGrl){
    
    cat("Transforming GRangesList to GRanges...\n")
    lifted <- flatGrl(lifted)
    return(lifted)
    
  }else{
    
    cat("Returning GRangesList...\n")
    return(lifted)
  }
}

# Performs GSEA by Fischer's exact test based on vector inputs
simpliGSEA <- function(genes, universe, annot.genes, annot.terms){
  
  if(length(annot.genes) != length(annot.terms)){
    
    stop("The 'annot.genes' and 'annot'terms' vectors must have same length!")
  }
  
  if("factor" %in% c(class(genes), class(universe), class(annot.genes), class(annot.terms))){
    
    stop("This function will not accept any factors as input!")
  }
  
  if(!all(genes %in% universe)){
    
    stop("Some provided 'genes' not found in 'universe'!")
  }
  
  if(!all(universe %in% annot.genes)){
    
    stop("Some 'universe' not in 'annot.genes'!")
  }
  
  if(length(unique(genes)) < length(genes)){
    
    cat("Removing duplicate IDs in 'genes'...")
    genes <- unique(genes)
  }
  
  if(length(unique(universe)) < length(universe)){
    
    cat("Removing duplicate IDs in 'universe'...")
    universe <- unique(universe)
  }
  
  # No need to test terms that do not show up in universe
  annot.universe <- annot.genes %in% universe
  annot.genes <- annot.genes[annot.universe]
  annot.terms <- annot.terms[annot.universe]
  
  # Initialize result container
  enrichment <- vector("numeric", length(unique(annot.terms)))
  names(enrichment) <- unique(annot.terms)
  
  for(term in unique(annot.terms)){
    
    # TRUE if has term
    A <- factor(universe %in% annot.genes[annot.terms == term], levels = c(TRUE, FALSE))
    
    # TRUE if in set
    G <- factor(universe %in% genes, levels = c(TRUE, FALSE))
    
    # See: https://www.biostars.org/p/111953/
    table <- table(G, A)
    print(table)
    
    enrichment[term] <- fisher.test(as.matrix(table))$p.value
    
    cat(term, "Fisher's exact test p-value:", enrichment[term], "\n\n")
  }
  
  return(enrichment)
}

###########################################################
### Functions to assist in data visualization

# Easily plot multiple graphs within the same window
# NOTE: Adapted from: http://www.cookbook-r.com/
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL){
  
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots <- length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if(is.null(layout)){
    
    # Make the panel
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if(numPlots == 1){
    
    print(plots[[1]])
    
  }else{
    
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for(i in 1:numPlots){
      
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}

###########################################################
### Functions to assist in cluster computing

# Finds package if not already installed
demand <- function(package){
  
  if(!is.element(package, installed.packages()[,1])){
    
    try({
      
      cat("Looking for package at CRAN...\n")
      install.packages(package)
    })
  }
  
  if(!is.element(package, installed.packages()[,1])){
    
    cat("Looking for package at Bioconductor...\n")
    source("https://bioconductor.org/biocLite.R")
    biocLite(package, suppressUpdates = TRUE)
  }
  
  if(!is.element(package, installed.packages()[,1])){
    
    try({
      
      cat("Looking for package at R-Forge...\n")
      install.packages(package, repos = "http://R-Forge.R-project.org")
    })
  }
  
  if(is.element(package, installed.packages()[,1])){
    
    require(package, character.only = TRUE)
    
  }else{
    
    stop("Could not find requested package!")
  }
}

# Saves a character string as R script
# NOTE: Returns script location
writeR <- function(R, folder = tempdir(), file = paste0(basename(folder), ".R")){
  
  script <- paste0(folder, "/", file)
  file.create(script)
  fileConn <- file(script)
  writeLines(R, fileConn)
  close(fileConn)
  
  return(script)
}

# Sends a Linux command to PBS queue via 'qsub'
# NOTE: More qsub arguments to come
qsub <- function(command, N = "name", M = "thom@tpq.me", m = "abe"){
  
  # Write script to execute command
  bash <- paste0("#! /bin/bash", "\n",
                 "#PBS -N ", N, "\n",
                 "#PBS -M ", M, "\n",
                 "#PBS -m ", m, "\n",
                 "\n", command)
  
  # Send script to queue
  run <- paste0("echo \"", bash, "\" | qsub")
  system(run)
}
