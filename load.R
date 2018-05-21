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
gpp <- lapply(db, genesPerPathway)
ppg <- lapply(db, pathwaysPerGene)

l2df <- function(x){
  l <- unlist(x, use.names = FALSE)
  data.frame(x = l,
                   db = rep(names(x), lengths(x)),
                   pathway = unlist(sapply(x, names), use.names = FALSE))
}
df_gpp <- l2df(gpp)
df_ppg <- l2df(pp)
colnames(df_ppg) <- c("ppg", "db", "gene")

ggplot(df[df$db == "c1.all", ]) +
  geom_density(aes(seq(1, max(df[df$db == "c1.all", "gpp"])), gpp, color = db)) +
  theme_bw()

gpp_t <- lapply(gpp, function(x){
  table(x)
})


t_num <- as.numeric(unique(unlist(sapply(gpp_t, names), use.names = FALSE)))
t_num <- t_num[order(t_num)]
m <- matrix(NA, nrow = max(t_num), ncol = length(gpp),
            dimnames = list(seq_len(max(t_num)), names(gpp)))

for (db in names(gpp_t)){
  db_m <- as.matrix(gpp_t[[db]])
  m[rownames(db_m), db] <- db_m
}
m <- as.data.frame(m)
m$size <- as.numeric(rownames(m))
ggplot(m) +
  geom_point(aes(size, c1.all), col = "blue") +
  geom_point(aes(size, c2.all), col = "red") +
  geom_point(aes(size, c3.all), col = "pink") +
  geom_point(aes(size, c4.all), col = "green") +
  geom_point(aes(size, c6.all), col = "grey") +
  # geom_point(aes(size, c7.all), col = "lightblue") +
  geom_point(aes(size, h.all), col = "orange") +
  xlab("Size of pathways") +
  ylab("Number of pathways") +
  theme_bw()

ggplot(m) +
  # geom_point(aes(size, c2.all), col = "blue") +
  geom_point(aes(size, c2.cp.biocarta), col = "red") +
  geom_point(aes(size, c2.cp.kegg), col = "pink") +
  geom_point(aes(size, c2.cp.reactome), col = "green") +
  geom_point(aes(size, c2.cp), col = "grey") +
  scale_x_continuous(limits = c(0, 1150)) +
  xlab("Size of pathways") +
  ylab("Number of pathways") +
  theme_bw()

paths <- as.data.frame(lengths(gpp))
colnames(paths) <- "Pathways"
genes <- as.data.frame(lengths(ppg))
colnames(genes) <- "Genes"
db2 <- cbind(paths, genes)
db2$DB <- rownames(db2)

ggplot(db2) +
  geom_bar(aes(DB, Pathways), stat = "identity") +
  theme_classic()

gpp_m <- sapply(gpp, '[', seq(max(lengths(gpp))))
gpp_m
