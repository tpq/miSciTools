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
  
  cat("Downloading", over.chain, "as temporary file...\n")
  url <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/", from, "/liftOver/",
                from, "To", capitalize(to), ".over.chain.gz")
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
