library(miSciTools)
data(ctd)

# For each chemical-gene link
res <- lapply(1:nrow(ctd), function(i){

  numTicks <<- miSciTools::progress(i = i, k = nrow(ctd), numTicks = numTicks)

  s.i <- strsplit(ctd$InteractionActions[i], "\\|")
  s.i <- lapply(s.i, strsplit, "\\^")[[1]]

  # For each InteractionAction (of form x^y)
  l.i <- lapply(s.i, function(x){

    df.i <- data.frame(x[1]) # name entry by x^
    colnames(df.i) <- x[2] # name column by ^y
    df.i
  })

  do.call("cbind", l.i)
})

# Save
ctd.actions <- do.call(plyr::rbind.fill, res)
ctd.wide <- cbind(ctd, ctd.actions)
devtools::use_data(ctd.wide)
