library("GSEABase")

gmt <- list.files(path = "data/", full.names = TRUE, pattern = ".gmt")
sapply(gmt, function(x){
  paths2Genes <- geneIds(getGmt(x,
                                geneIdType=EntrezIdentifier()))

  genes <- unlist(paths2Genes, use.names = FALSE)
  pathways <- rep(names(paths2Genes), lengths(paths2Genes))
  split(pathways, genes)
})
