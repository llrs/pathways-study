library("GSEABase")

# Load the files
gmt <- list.files(path = "data/", full.names = TRUE, pattern = "entrez.gmt")
db <- sapply(gmt, function(x){geneIds(getGmt(x))})
names(db) <- gsub("data//(.*)\\.v6\\.1\\.entrez\\.gmt",
                         replacement =  "\\1", names(db))
databases <- sapply(gmt, function(x){
  paths2Genes <- geneIds(getGmt(x))

  genes <- unlist(paths2Genes, use.names = FALSE)
  pathways <- rep(names(paths2Genes), lengths(paths2Genes))
  split(pathways, genes)
})
names(databases) <- gsub("data//(.*)\\.v6\\.1\\.entrez\\.gmt",
                         replacement =  "\\1", names(databases))

genesPerPathway <- sapply(db, lengths)

all <- sapply(names(databases), function(y){
  sapply(databases[[y]], function(x){
    genesPerPathway[[y]][x]
  })
})

pathwaysPerGene <- sapply(databases, lengths)

x <- "c2.cp.reactome"
plot(table(pathwaysPerGene[[x]]))
plot(table(genesPerPathway[[x]]))



# plot(density(genesPerPathway[["c2.cp"]]), col = "red")
lines(density(genesPerPathway[["c2.cp.reactome"]]), col = "green", xlab = "Pathways",
      main = "Comparison of the size of genes per Pathways for multiple databases")
# lines(density(genesPerPathway[["c2.all"]]))
lines(density(genesPerPathway[["c2.cp.kegg"]]), col = "blue")
lines(density(genesPerPathway[["c2.cp.biocarta"]]), col = "pink") # Out of charts
# lines(density(genesPerPathway[["c2.cgp"]]), col = "brown")
maxs <- sapply(genesPerPathway[c("c2.cp.biocarta", "c2.cp.reactome", "c2.cp.biocarta")],
       function(x){max(table(x))})
abline(v = maxs, col = c("pink", "green", "blue")[order(maxs)])
range(maxs) # The most expected site of the pathway regardless the database
plot(function(x) dt(x, df = 25, ncp = 10), -3, 800,
     main = "Non-central t - Density", yaxs = "i")


lines(density(pathwaysPerGene[["c2.cp.reactome"]]), col = "green",  xlab = "Genes",
     main = "Comparison of the number of Pathways per gene for multiple databases")
lines(density(pathwaysPerGene[["c2.cp.kegg"]]), col = "blue") # Out of charts
# lines(density(pathwaysPerGene[["c2.cp"]]), col = "red")
# lines(density(pathwaysPerGene[["c2.all"]]))
lines(density(pathwaysPerGene[["c2.cp.biocarta"]]), col = "pink") # Out of charts
# lines(density(pathwaysPerGene[["c2.cgp"]]), col = "brown")
maxs <- sapply(pathwaysPerGene[c("c2.cp.biocarta", "c2.cp.reactome", "c2.cp.biocarta")],
               function(x){a <- table(x)
                 as.numeric(names(a)[which.max(a)])})
abline(v = maxs, col = c("pink", "green", "blue")[order(maxs)])
range(maxs) # If a gene has a pathway the most common are between 1 and 3
# It looks like a t student distribution
plot(function(x) dt(x, df = 3, ncp = 4), -3, 100,
     main = "Non-central t - Density", yaxs = "i")

# Total ~pathways
(pathways <- length(unique(unlist(sapply(databases, function(x){
  names(x)
})))))

# Total genes
(genes <- length(unique(unlist(sapply(databases, function(x){
  unlist(x)
})))))

# Analyse the pathways
library("BioCor")

for (f in seq_along(databases)){
  o <- mpathSim(unique(unlist(databases[[f]])),
                info = databases[[f]], method = NULL)
  write.csv(o, file =  paste0(names(databases[f]), ".csv"))
}

paths <- list.files(pattern = "*.csv")
pathwaySim <- lapply(paths, read.csv)

gse <- read.csv("c7.all.csv", row.names = 1)
gse <- as.matrix(gse)
library("Hmisc")
Ecdf(gse[upper.tri(gse)])

combined <- combineSources(databases)

o <- mpathSim(unique(unlist(combined)), info = combined, method = NULL)
