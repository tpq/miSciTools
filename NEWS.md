## miSciTools 1.1.3
---------------------
* Add `fs` methods
    * Migrate `fsPathClassRFE` and `fsSelbal` from `exprso`

## miSciTools 1.1.2
---------------------
* Add `fs` methods
    * New `fsPropd` module selects features using `propr::propd`
    * New `fsPRA` module selects features using `propr::pra`
* Add methods
    * New `rda.reconstruct` reconstructs data from RDA object
    * New `wwtest` resembles `tttest` but non-parametric

## miSciTools 1.1.1
---------------------
* Revise `getRecount` function
    * Make `getRecount` more reliable

## miSciTools 1.1.0
---------------------
* Add methods
    * New `getRecount` retrieves data from recount2 server

## miSciTools 1.0.9
---------------------
* Add methods
    * New `quick.DESeq2` function performs formula-based DE analysis
    * New `wide2long` and `long2wide` functions convert between data styles
    * New `tttest` function to make t-test tables from long data

## miSciTools 1.0.8
---------------------
* Revise `aptly` function
    * Delineate mean and total lengths for found structures
    * Allow user to append additional RNAfold arguments
    * Characterize left and right bulges

## miSciTools 1.0.7
---------------------
* Implemented `aptly` function
    * Procedurally characterizes RNA sequences based on secondary structure
    * Fixed bug caused by absence of any secondary structure

## miSciTools 1.0.6
---------------------
* Implemented `fastafolder` function
    * Calculates the dissimilarity between sequences based on secondary structure
* Added `packageCheck` function
    * Triggers an error if a package is not installed

## miSciTools 1.0.5
---------------------
* Update `ctd` data and remove `NA` GeneIDs

## miSciTools 1.0.4
---------------------
* Implemented `quick.edgeR` function
    * Run `edgeR` Exact Test for differential expression analysis
* Added `kegg` data
    * Details KEGGREST-gene relationships

## miSciTools 1.0.3
---------------------
* Have `simpliGSEA` pass `alternative` argument to `fisher.test`

## miSciTools 1.0.2
---------------------
* Have `migraph` functions coerce arguments as character input
* Fix errors in documentation and package imports
* Fix `simpliGSEA` pruning so it actually works
* Add code for `progress` bar function

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
