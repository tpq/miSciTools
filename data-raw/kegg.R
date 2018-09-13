###########################################################
### Build Data for KEGGREST

list <- read.delim("http://rest.kegg.jp/list/pathway/hsa", header = FALSE)
head(list)
colnames(list) <- c("path", "kegg")

link <- read.delim("http://rest.kegg.jp/link/pathway/hsa", header = FALSE)
head(link)
colnames(link) <- c("entrez", "path")

merge <- merge(list, link)
merge$path <- unlist(lapply(strsplit(as.character(merge$path), split = ":"), function(x) x[2]))
merge$entrez <- unlist(lapply(strsplit(as.character(merge$entrez), split = ":"), function(x) x[2]))
merge$kegg <- unlist(lapply(strsplit(as.character(merge$kegg), split = " - "), function(x) x[1]))
kegg <- merge[, c("entrez", "path", "kegg")]

devtools::use_data(kegg)
