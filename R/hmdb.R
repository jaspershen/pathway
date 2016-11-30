##抓取HMDB通路信息
library(XML)
library(RCurl)

hmdb.main <- "http://www.hmdb.ca/pathways?page="
hmdb.main <- paste(hmdb.main, c(2:46), sep = "")
hmdb.main <- c("http://www.hmdb.ca/pathways", hmdb.main)

##从HMDB主页上抓取代谢通路的url
path.name <- list()
metabolite.id <- list()
spec <- list()
path.class <- list()
for (i in 40:length(hmdb.main)) {
  cat(paste("page",i))
  cat(":")
  URL = getURL(hmdb.main[i])
  doc<-htmlParse(URL,encoding="utf-8")
  xpath1 <- "//div[@class='panel-heading']"
  node1 <- getNodeSet(doc, xpath1)
  pathway.name <- sapply(node1, xmlValue)

  cat(paste(length(pathway.name), "pathways"))
  cat("\n")

  path.name[[i]] <- pathway.name

  xpath2 <- "//div[@class='panel-body']"
  node2 <- getNodeSet(doc, xpath2)

  metabolite <- sapply(node2, xmlValue)
  metabolite <- unname(sapply(metabolite, function(x) {gsub("Show", " ", x)}))

  idx <- sapply(metabolite, function(x) {gregexpr("HMDB[0-9]{5}", x)})

  met.id <- list()
  for (j in 1:length(idx)) {
    id <- NULL
    for (k in 1:length(idx[[j]])) {
      id[k] <- substr(metabolite[j], idx[[j]][k], idx[[j]][k]+8)
    }
    met.id[[j]] <- id
  }

  metabolite.id[[i]] <- met.id

  xpath.a <- "//a[@class='link-out']/@href"
  node<-getNodeSet(doc, xpath.a)

  url1 <- sapply(node, as.character)
  url1 <- substr(url1, start = 1, stop = 29)
  url1 <- url1[!duplicated(url1)]


  ###获取通路的人种和类别
  species <- NULL
  metabolic <- NULL
  for (t in 1:length(url1)) {
    cat(paste("t:",t));cat(" ")
    URL = getURL(url1[t])
    doc <- htmlParse(URL,encoding="utf-8")
    xpath <- "//div[@class='species']/text()"
    node <- getNodeSet(doc, xpath)
    species[t] <- xmlValue(node[[1]])

    xpath <- "//div[@id='des_subject']/text()"
    node <- getNodeSet(doc, xpath)
    metabolic[t] <- xmlValue(node[[1]])

  }

  spec[[i]] <- species
  path.class[[i]] <- metabolic

}

getwd()

save(path.name, file = "path.name")
save(metabolite.id, file = "metabolite.id")
save(spec, file = "spec")
save(path.class, file = "path.class")


pathway.name <- NULL
metabolite.ID <- list()
species <- NULL
pathway.class <- NULL
for (i in 1:length(path.name)) {
  pathway.name <- c(pathway.name, path.name[[i]])
  metabolite.ID <- c(metabolite.ID, metabolite.id[[i]])
  species <- c(species, spec[[i]])
  pathway.class <- c(pathway.class, path.class[[i]])
}


pathway.class <- substr(x = pathway.class, 1, regexpr("\\\n", pathway.class)-1)



setwd("F:/pathway")

metabolite.name <- list()
for (i in 1:length(metabolite.ID)) {
  id <- metabolite.ID[[i]]
  idx <- match(id, hmdbdatabase[,1])
  name <- hmdbdatabase[idx,2]
  metabolite.name[[i]] <- name
}

metabolite.name
metabolite.ID
pathway.name
pathway.class
species

a <- unlist(lapply(metabolite.name, function(x) {paste(x, collapse = ";")}))
b <- unlist(lapply(metabolite.ID, function(x) {paste(x, collapse = ";")}))

idx <- grep("Metabolic", pathway.class)


metabolite.name <- metabolite.name[idx]
metabolite.ID <- metabolite.ID[idx]
pathway.name <- pathway.name[idx]
pathway.class <- pathway.class[idx]
species <- species[idx]



hmdb.pathway <- data.frame(pathway.name, pathway.class,a, b)[idx,]
write.csv(hmdb.pathway, "hmdb.pathway.csv")

a <- list()
for (i in 1:length(pathway.name)) {
  a[[i]] <- pathway.name[i]
}

pathway.name <- a

hmdb.met <- list(gs = metabolite.name, pathwaynames = pathway.name, id = metabolite.ID)
save(hmdb.met, file = "hmdb.met")




###抓取HMDB代谢物信息
library(XML)
library(RCurl)

hmdb.main <- "http://www.hmdb.ca/metabolites?c=hmdb_id&d=up&page="
hmdb.main <- paste(hmdb.main, c(2:1681), sep = "")
hmdb.main <- c("http://www.hmdb.ca/metabolites", hmdb.main)

##从HMDB主页上抓取代谢物的url
url <- NULL
for (i in 1:length(hmdb.main)) {
  cat(i)
  cat(" ")
  URL = getURL(hmdb.main[i])
  doc<-htmlParse(URL,encoding="utf-8")
  xpath <- "//a[@href]/@href"
  node<-getNodeSet(doc, xpath)
  url1 <- sapply(node, as.character)
  url1 <- url1[grep("metabolites/HMDB", url1)]
  url1 <- unique(url1)
  url <- c(url, url1)
}

url1 <- paste("http://www.hmdb.ca/",url, sep = "")
save(url1, file = "url1")

HMDB_ID <- NULL
Common_Name <- NULL
Synonyms <- NULL
Chemical_Formula <- NULL
Monoisotopic_Molecular_Weight <- NULL
IUPAC_Name <- NULL
Traditional_Name <- NULL
CAS_Registry_Number <- NULL
Origin <- NULL
Pathways <- NULL

for (i in 1:14000) {
  cat(paste(i,"/", length(url1)))
  cat("\n")
  URL = getURL(url1[i])
  doc <- htmlParse(URL,encoding="utf-8")
  xpath1 <- "//tr"
  node1 <- getNodeSet(doc, xpath1)
  node1 <- sapply(node1, xmlValue)

  HMDB_ID[i] <- gsub(pattern = "HMDB ID", replacement = "",node1[grep("HMDB ID", node1)])

  Common_Name[i] <- gsub("Common Name", "",node1[grep("Common Name", node1)])

  Synonyms <- gsub("SynonymsValueSource", "",node1[grep("Synonyms", node1)])
  Synonyms <- gsub("Generator", ";",Synonyms)
  Synonyms[i] <- gsub("HMDB", ";",Synonyms)

  Chemical_Formula[i] <- gsub("Chemical Formula", "",node1[grep("Chemical Formula", node1)])

  Monoisotopic_Molecular_Weight[i] <- gsub("Monoisotopic Molecular Weight", "",node1[grep("Monoisotopic Molecular Weight", node1)])

  IUPAC_Name[i] <- gsub("IUPAC Name", "",node1[grep("IUPAC Name", node1)])

  Traditional_Name[i] <- gsub("Traditional Name", "",node1[grep("Traditional Name", node1)])

  CAS_Registry_Number[i] <- gsub("CAS Registry Number", "",node1[grep("CAS Registry Number", node1)])

  Origin[i] <- gsub("Origin", "",node1[grep("Origin", node1)])

  path <- gsub("PathwaysNameSMPDB LinkKEGG Link", "",node1[grep("Pathways", node1)])
  Pathways[i] <-substr(path, 1, stop = regexpr("SMP", path)-1)

  ##每1000次保存一次
  if (i %in% seq(1000, 60000, by = 1000)) {
    cat("save data...")
    save(HMDB_ID,
         Common_Name,
         Synonyms,
         Chemical_Formula,
         Monoisotopic_Molecular_Weight,
         IUPAC_Name,
         Traditional_Name,
         CAS_Registry_Number,
         Origin,
         Pathways,
         file = "hmdb.data")
  }

}
