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

# Remove NAs
ctd <- ctd[!is.na(ctd$GeneID), ]
rownames(ctd) <- 1:nrow(ctd)

# Save
devtools::use_data(ctd)
