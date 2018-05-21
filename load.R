library("GSEAdv")
library("tidyverse")

theme_set(theme_bw())
# Load the files
gmt <- list.files(path = "data/", full.names = TRUE, pattern = "entrez.gmt")
db <- sapply(gmt, function(x){getGmt(x)})

names(db) <- gsub("data//(.*)\\.v6\\.1\\.entrez\\.gmt",
                         replacement =  "\\1", names(db))


do_all <- function(x, data) {
  dir.create(x, showWarnings = FALSE)
  # cd <- setwd(x)
  pdf(file.path(x, "figures.pdf"))
  gsc <- data[[x]]
  gpp <- genesPerPathway(gsc)
  ppg <- pathwaysPerGene(gsc)

  plot(table(gpp), main = x)
  plot(table(ppg), main = x)

  sPG <- sizesPerGene(gsc)
  sPP <- sizesPerPathway(gsc)

  plot(table(sPG), main = x)
  plot(table(sPP), main = x)

  pg <- data.frame("ppg"= ppg, "sPG" = sPG)
  p <- ggplot(pg) +
    geom_count(aes(ppg, sPG)) +
    geom_abline(slope = 1, intercept = 0)
  print(p)
  hist(sPG/ppg, main = "Unique sizes per")

  gp <- data.frame("ppg"= gpp, "sPG" = sPP)
  p <- ggplot(gp) +
    geom_count(aes(gpp, sPP)) +
    geom_abline(slope = 1, intercept = 0)
  print(p)
  hist(sPP/gpp, main = "Unique sizes per")

  sG <- sizeGenes(gsc)
  sP <- sizePathways(gsc)

  sG_d <- sweep(sG, 2, ppg, "/")
  uscPG <- unlist(sG_d)
  uscPG <- uscPG[uscPG != 0]
  hist(uscPG)

  sP_d <- sweep(sP, 2, gpp, "/")
  uscPP <- unlist(sP_d)
  uscPP <- uscPP[uscPP != 0]
  hist(uscPP, main = "Distribution of probabilities by pathway")
  diff0 <- function(x) {sum(x != 0, na.rm = TRUE)}

  cPG <- condPerGenes(gsc)
  cPG0 <- apply(cPG, 2, diff0)
  sG_d <- sweep(cPG, 2, cPG0, "*")
  uscPG <- unlist(sG_d)
  uscPG <- uscPG[uscPG != 0]
  hist(uscPG)

  cPG <- condPerPathways(gsc)
  cPG0 <- apply(cPG, 2, diff0)
  sG_d <- sweep(cPG, 2, cPG0, "*")
  uscPG <- unlist(sG_d)
  uscPG <- uscPG[uscPG != 0]
  hist(uscPG, main = "Distribution of probabilities (of sizes) by pathway")
  dev.off()
  # on.exit(setwd(cd))
}

sapply(names(db), do_all, data = db)
