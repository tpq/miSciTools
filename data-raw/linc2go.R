###########################################################
### Build Linc2GO GeneSetCollection Object for GSEA

# Link to Linc2GO
# http://www.bioinfo.tsinghua.edu.cn/~liuke/Linc2GO/

# Build GeneSetCollection object, use by:
# (1) Category::GSEAGOHyperGParams()
# (2) GOstats::hyperGTest()

library(utils)
library(AnnotationDbi)
library(GSEABase)
library(devtools)

# Download BP
url <- "http://www.bioinfo.tsinghua.edu.cn/~liuke/Linc2GO/LincNetwork_BP_Enrichment.csv"
file.bp <- tempfile(fileext = ".csv")
utils::download.file(url, file.bp)

# Read BP and build GOFrame
linc2go.bp <- read.csv(file.bp, header = FALSE)
linc2go.bp.frame <- linc2go.bp[, c(1, 2, 3)]
colnames(linc2go.bp.frame) <- c("gene_id", "go_id", "Evidence")
linc2go.bp.frame <- cbind(linc2go.bp.frame, data.frame("Ontology" = "BP"))

# Download MF
url <- "http://www.bioinfo.tsinghua.edu.cn/~liuke/Linc2GO/LincNetwork_MF_Enrichment.csv"
file.mf <- tempfile(fileext = ".csv")
utils::download.file(url, file.mf)

# Read MF and build GOFrame
linc2go.mf <- read.csv(file.mf, header = FALSE)
linc2go.mf.frame <- linc2go.mf[, c(1, 2, 3)]
colnames(linc2go.mf.frame) <- c("gene_id", "go_id", "Evidence")
linc2go.mf.frame <- cbind(linc2go.mf.frame, data.frame("Ontology" = "MF"))

# Combine BP, MF GOFrames
linc2go.frame <- rbind(linc2go.bp.frame, linc2go.mf.frame)
linc2go.frame <- linc2go.frame[, c(2, 3, 1)]
linc2go.frame$Evidence <- "RCA"
linc2go.goframe <- AnnotationDbi::GOFrame(linc2go.frame, organism = "Homo sapiens")

# Turn GOFrame into GOAllFrame
linc2go.goall <- AnnotationDbi::GOAllFrame(linc2go.goframe)

# Turn GOAllFrame into GeneSetCollection
# See also: ?GSEABase::KEGGCollection
linc2go <- GSEABase::GeneSetCollection(linc2go.goall, setType = GOCollection())

devtools::use_data(linc2go)
