library(XML)
library(RCurl)
##从kegg主页上抓取代谢通路的url
URL = getURL("http://www.genome.jp/kegg/pathway.html#global")
doc<-htmlParse(URL,encoding="utf-8")
xpath.a <- "//a/@href"
node<-getNodeSet(doc, xpath.a)
url1 <- sapply(node, as.character)

xpath.b <- "//a[@href]"
name <- getNodeSet(doc, xpath.b)
name <- sapply(name, xmlValue)

name2 <- name[59:247]
url2 <- url1[59:247]

url3 <- url2[grep("show", url2)]

pathwat.name <- NULL
metabolite.id <- list()
metabolite.name <- list()
for (i in 1:length(url3)) {
  cat(paste(i,"/",length(url3)))
  cat("\n")
  URL <- paste("http://www.genome.jp", url3[i], sep = "")
  URL = getURL(URL)
  doc<-htmlParse(URL,encoding="utf-8")
  xpath <- "//option[@value='hsa']"
  node<-getNodeSet(doc, xpath)
  if (length(node) ==0 ) {
    cat("No human pathwat.")
    next()
  }else{
    URL <- paste("http://www.genome.jp", url3[i], sep = "")
    URL <- gsub(pattern = "map=map", replacement = "map=hsa", x = URL)
    doc<-htmlParse(URL,encoding="utf-8")
    xpath1 <- "//title"
    node<-getNodeSet(doc, xpath1)
    pathway.name[i] <- xmlValue(node[[1]])
    pathway.name[i] <- substr(pathway.name[i], start = 2, stop = nchar(pathway.name[i])-1)

    xpath2 <- "//area[@shape='circle']/@title"
    node<-getNodeSet(doc, xpath2)
    metabolite <- lapply(node, function(x) as.character(x))
    metabolite.name[[i]] <- substr(metabolite, start = 9, nchar(metabolite)-1)
    metabolite.id[[i]] <- substr(metabolite, start = 1, stop = 6)
  }
}



##
idx <- which(!is.na(pathway.name))
pathway.name1 <- pathway.name[idx]
metabolite.id1 <- metabolite.id[idx]
metabolite.name1 <- metabolite.name[idx]

pathway.name2 <- pathway.name1[-c(83,84)]
metabolite.id2 <- metabolite.id1[-c(83,84)]
metabolite.name2 <- metabolite.name1[-c(83,84)]


met.name <- NULL
met.id <- NULL
path.name <- NULL
for(i in 1:length(pathway.name2)) {
  met.name[i] <- paste(metabolite.name2[[i]], collapse = ";")
  met.id[i] <- paste(metabolite.id2[[i]], collapse = ";")
  path.name[i] <- gsub(pattern = "KEGG PATHWAY: ", "", pathway.name2[i])
  path.name[i] <- substr(path.name[i], start = 1, stop = nchar(path.name[i])-23)
}


kegg <- data.frame(path.name, met.name, met.id)
write.csv(kegg, "kegg.csv", row.names = F)

save(path.name, file = "path.name")
save(met.name, file = "met.name")
save(met.id, file = "met.id")

kegg.met <- list()
kegg.met[[2]] <- sapply(path.name, list)
kegg.met[[1]] <- metabolite.name2
kegg.met[[3]] <- metabolite.id2

names(kegg.met) <- c("gs", "pathwaynames", "metid")

save(kegg.met, file = "kegg.met")


