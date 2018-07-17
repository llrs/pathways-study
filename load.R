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
  cd <- setwd(x)
  # pdf(file.path(x, "figures.pdf"))
  gsc <- data[[x]]
  gpp <- genesPerPathway(gsc)
  ppg <- pathwaysPerGene(gsc)
  png("gpp.png")
  plot(table(gpp), main = "Genes per pathways")
  dev.off()
  png("ppg.png")
  plot(table(ppg), main = "Pathways per genes")
  dev.off()

  sPG <- sizesPerGene(gsc)
  sPP <- sizesPerPathway(gsc)

  png("sPG.png")
  plot(table(sPG), main = "Different size of pathways per gene")
  dev.off()
  png("sPP.png")
  plot(table(sPP), main = "Different size of genes per pathway")
  dev.off()
  gpg <- genesPerGene(gsc)
  pg <- data.frame("ppg"= ppg, "sPG" = sPG, "gpg" = gpg)
  p <- ggplot(pg) +
    geom_count(aes(ppg, sPG)) +
    xlab("Pathways per Gene") +
    geom_abline(slope = 1, intercept = 0) +
    ylab("Pathways of different size")
  ggsave(p, filename = "pg1.png")
  p2 <- ggplot(pg) +
    geom_count(aes(ppg, gpg)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Pathways per Gene") +
    ylab("Genes related to gene")
  ggsave(p2, filename = "pg2.png")
  p3 <- ggplot(pg) +
    geom_count(aes(sPG, gpg)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Pathways of different size")
    ylab("Genes related to gene")
  ggsave(p3, filename = "pg1.png")

  png("ratio_ppg.png")
  hist(sPG/ppg, main = "Unique sizes per gene")
  dev.off()

  gp <- data.frame("ppg"= gpp, "sPG" = sPP)
  p <- ggplot(gp) +
    geom_count(aes(gpp, sPP)) +
    geom_abline(slope = 1, intercept = 0) +
    ylab("size pathways") +
    xlab("Genes per pathway")
  ggsave(p, filename = "gp.png")

  png("ratio_gpp.png")
  hist(sPP/gpp, main = "Unique sizes per")
  dev.off()

  sG <- sizeGenes(gsc)
  sP <- sizePathways(gsc)

  sG_d <- sweep(sG, 2, ppg, "/")
  uscPG <- unlist(sG_d)
  uscPG <- uscPG[uscPG != 0]

  png("uscPG.png")
  hist(uscPG, main = "Distribution of probabilities of number of genes in pathways")
  dev.off()

  sP_d <- sweep(sP, 2, gpp, "/")
  uscPP <- unlist(sP_d)
  uscPP <- uscPP[uscPP != 0]

  png("uscPP.png")
  hist(uscPP, main = "Distribution of probabilities of number of pathways per gene")
  dev.off()

  diff0 <- function(x) {sum(x != 0, na.rm = TRUE)}

  cPG <- condPerGenes(gsc)
  png("cPG.png")
  hist(cPG, main = "condPerGenes")
  dev.off()

  cPG0 <- apply(cPG, 2, diff0)
  sG_d <- sweep(cPG, 2, cPG0, "*")
  uscPG <- unlist(sG_d)
  uscPG <- uscPG[uscPG != 0]
  png("uscPG.png")
  hist(uscPG, main = "Distribution of probabilities (of sizes) by gene")
  dev.off()

  cPP <- condPerPathways(gsc)
  png("cPP.png")
  hist(cPP, main = "condPerPathways")
  dev.off()
  cPP0 <- apply(cPP, 2, diff0)
  sG_d <- sweep(cPP, 2, cPP0, "*")
  uscPP <- unlist(sG_d)
  uscPP <- uscPP[uscPP != 0]
  png("uscPP.png")
  hist(uscPP, main = "Distribution of probabilities (of sizes) by pathway")
  dev.off()
  setwd(cd)
}

sapply(names(db), do_all, data = db)
