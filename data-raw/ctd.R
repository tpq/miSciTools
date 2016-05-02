###########################################################
### Build Data for Comparative Toxicogenomics Database

# Link to many toxigenomic databases:
# http://alttox.org/resource-center/databases/

# Link to CTD:
# nar.oxfordjournals.org/content/41/D1/D1104.full
# http://ctdbase.org/

library(utils)
library(R.utils)
library(devtools)

# Download
url <- "http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz"
file.gz <- tempfile(fileext = ".gz")
utils::download.file(url, file.gz)

# Unzip
file.csv <- tempfile(fileext = ".csv")
R.utils::gunzip(filename = file.gz, destname = file.csv)

# Read
ctd <- read.csv(file.csv, skip = 28, header = FALSE, stringsAsFactors = FALSE)
colnames(ctd) <-
  c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID", "GeneForms",
    "Organism", "OrganismID", "Interaction", "InteractionActions", "PubMedIDs")

# Check
if(!all(ctd$GeneID[1:8] == c(4149, 4149, 4609, 4609, 4609, 4609, 4609, 3784))){

  stop("Import check failed. A 'CTD' update may have caused this?")
}

# Save
devtools::use_data(ctd)
