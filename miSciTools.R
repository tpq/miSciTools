###########################################################
### Functions to assist in proportionality analyses

# Calculate phi from a data.frame of feature counts
phit <- function(counts, symmetrize = TRUE){
  
  require(compositions)
  
  # Replace zeroes with next smallest number
  counts[counts == 0] <- unique(sort(as.matrix(counts)))[2]
  
  # Calculate the variance of the log-ratio ("variation array")
  counts.vlr <- variation(acomp(t(counts)))
  colnames(counts.vlr) <- rownames(counts)
  rownames(counts.vlr) <- rownames(counts)
  
  # Calculate feature variance across clr transformed treatments
  counts.clr <- as.data.frame(clr(t(counts)))
  counts.clr.var <- apply(counts.clr, 2, var)
  
  # Sweep out feature clr variance from the variation array
  counts.phi <- sweep(counts.vlr, 2, counts.clr.var, FUN = "/")
  
  # Symmetrize matrix if symmetrize = TRUE
  if(symmetrize) counts.phi <- propSym(counts.phi)
  
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
    distr[begin:end] <- propTri(phi.i)
    rm(phi.i)
  }
  
  cat("Fitting phi to distribution...\n")
  fit <- ecdf(distr)
  rm(distr)
  
  if(returnPval){
    
    cat("Calculating all phi for actual counts...\n")
    phi <- phit(counts)
    raw <- propRaw(phi)
    
    cat("Using 'fit' to convert phi into pval...\n")
    pval <- fit(raw$prop)
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

# Calculates proportionality using p coefficient from Erb 2016
# NOTE: IF 'ivar' = NULL, divide variation array by clr transformed variance
# NOTE: ELSE, divide variation array by alr transformed variance
perb <- function(counts, ivar = NULL){
  
  require(compositions)
  
  # Replace zeroes with next smallest number
  counts[counts == 0] <- unique(sort(as.matrix(counts)))[2]
  
  # Calculate the variance of the log-ratio ("variation array")
  counts.vlr <- variation(acomp(t(counts)))
  colnames(counts.vlr) <- rownames(counts)
  rownames(counts.vlr) <- rownames(counts)
  
  if(!is.null(ivar)){
    
    # Calculate feature variance across alr transformed treatments
    counts.vlr <- counts.vlr[-ivar, -ivar] # returns one less dimension
    counts.alr <- as.data.frame(alr(t(counts), ivar = ivar)) # returns one less dimension
    counts.var <- apply(counts.alr, 2, var)
    
  }else{
    
    # Calculate feature variance across clr transformed treatments
    counts.clr <- as.data.frame(clr(t(counts)))
    counts.var <- apply(counts.clr, 2, var)
  }
  
  # Divide variation array by sum of feature variances
  for(i in 1:ncol(counts.vlr)){
    for(j in 1:nrow(counts.vlr)){
      counts.vlr[i, j] <- counts.vlr[i, j] / (counts.var[i] + counts.var[j])
    }
  }
  
  # Calculate: p = 1 - (var(x - y))/(var(x) + var(y))
  counts.prop <- 1 - counts.vlr
  
  return(counts.prop)
}

# Calculate p(x) based on a NULL distribution of p from Erb 2016
perbDistr <- function(counts, ivar = NULL, iter = 10, iterSize = nrow(counts) - !is.null(ivar), returnPval = TRUE){
  
  distr <- vector("numeric", iter * iterSize * (iterSize - 1) / 2)
  
  for(i in 1:iter){
    
    # Sample counts matrix
    cat(paste0("Calculating all perb for iter ", i, "...\n"))
    
    # Handle alr properly
    if(!is.null(ivar)){
      
      null.i <- apply(counts[-ivar, ], 2, sample, iterSize)
      fixed <- counts[ivar, ]
      null.i <- rbind(null.i, fixed)
      ivar.i <- nrow(null.i)
      
    }else{
      
      null.i <- apply(counts, 2, sample, iterSize)
      ivar.i <- NULL
    }
    
    # Calculate perb matrix
    perb.i <- perb(null.i, ivar.i)
    begin <- (i - 1) * iterSize * (iterSize - 1) / 2 + 1
    end <- (i - 1) * iterSize * (iterSize - 1) / 2 + iterSize * (iterSize - 1) / 2
    
    # Save lower triangle
    distr[begin:end] <- propTri(perb.i)
    rm(perb.i)
  }
  
  cat("Fitting phi to distribution...\n")
  fit <- ecdf(distr)
  rm(distr)
  
  if(returnPval){
    
    cat("Calculating all phi for actual counts...\n")
    prop <- perb(counts, ivar = ivar)
    raw <- propRaw(prop)
    
    cat("Using 'fit' to convert phi into pval...\n")
    pval <- fit(raw$prop)
    pval[pval >= .5] <- 1 - pval[pval >= .5] # make 2-tails
    pval <- pval * 2 # scale 0 to 1
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

# Retrieve phi (or p) for each feature pair as data.frame
propRaw <- function(prop){
  
  index.i <- vector("numeric", length = (nrow(prop) - 1)*nrow(prop)/2)
  index.j <- vector("numeric", length = (nrow(prop) - 1)*nrow(prop)/2)
  index.prop <- vector("numeric", length = (nrow(prop) - 1)*nrow(prop)/2)
  counter <- 1
  
  for(j in 2:nrow(prop)){
    
    for(i in 1:(j-1)){
      
      index.i[counter] <- i
      index.j[counter] <- j
      index.prop[counter] <- prop[j, i]
      counter <- counter + 1
    }
  }
  
  result <- data.frame("feature1" = rownames(prop)[index.i],
                       "feature2" = rownames(prop)[index.j],
                       "prop" = index.prop,
                       stringsAsFactors = FALSE)
  
  final <- result[order(result$prop),]
  rownames(final) <- 1:nrow(final)
  
  return(final)
}

# Retrieve the lower triangle of a phi (or p) matrix
propTri <- function(prop){
  
  result <- vector("numeric", length = (nrow(prop) - 1)*nrow(prop)/2)
  counter <- 1
  
  for(j in 2:nrow(prop)){
    
    for(i in 1:(j-1)){
      
      result[counter] <- prop[j, i]
      counter <- counter + 1
    }
  }
  
  return(result)
}

# Symmetrize the asymmetric phi matrix
propSym <- function(prop){
  
  for(j in 2:nrow(prop)){
    
    for(i in 1:(j-1)){
      
      prop[i, j] <- prop[j, i]
    }
  }
  
  return(prop)
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

# NOTE: Supply fx as e.g. maxRAM(function() vlr(t(counts)))
maxRAM <- function(fx){
  
  start <- gc(reset = TRUE)
  start <- start["Vcells", 6]
  run <- fx
  run()
  end <- gc(TRUE)
  end <- end["Vcells", 6]
  final <- end - start
  names(final) <- "max.Mb"
  gc(reset = TRUE)
  return(final)
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
