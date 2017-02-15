## miSciTools 1.0.2
---------------------
* Have `migraph` functions coerce arguments as character input
* Fix errors in documentation and package imports

## miSciTools 1.0.1
---------------------
* Modified `simpliGSEA` function
    * Now prunes genes missing in universe or database
* Implemented `migraph` functions
    * New `migraph.add` adds vertices or edges to graph
    * New `migraph.color` colors vertices or edges
    * New `migraph.clean` cleans graph to plot

## miSciTools 1.0.0
---------------------
* Implemented `lift` function
    * "Lifts over"" coordinates from one genome build to another
* Implemented `simpliGSEA` function
    * Performs gene set enrichment analysis for non-DAG
* Implemented `multiplot` function
    * Quickly plots multiple figures in the same figure
    * Source: http://www.cookbook-r.com/
* Implemented `demand` function
    * Looks for and installs package from multiple sources
* Implemented `asTempObj` function
    * Sets value as a random but unique variable name
* Implemented `writeR` function
    * Writes an R script as a temporary file
* Implemented `qsub` function
    * Sends a system command through qsub
* Implemented `insert` function
    * Inserts a new line of text into an existing file
* Implemented `cache.key` function
    * Generates a unique hash key from a list of objects
* Implemented `cache.save` function
    * Saves a value as raw to an SQLite database
* Implemented `cache.load` function
    * Loads a value as raw from an SQLite database
* Added `ctd` data
    * Details chemical-gene interactions
* Added `linc2go` data
    * A `GeneSetCollection` for lincRNA
