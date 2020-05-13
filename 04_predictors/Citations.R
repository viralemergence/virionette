#install.packages("easyPubMed")
library("easyPubMed")
library(tidyverse)

tree <- readRDS("C:/Users/cjcar/Documents/GitHub/virionette/04_predictors/Full Supertree.rds")
tree$tip.label <- gsub("_"," ",tree$tip.label)

counter <- function(name) {
  as.numeric(as.character(get_pubmed_ids(gsub(' ','-',name))$Count))
}

citations <- c()

for(i in 1:length(tree$tip.label)) {
  citations[i] <- counter(tree$tip.label[i])
  print(i)
}

cites <- data.frame(name = tree$tip.label, cites = citations)

####

read_csv('~/GitHub/virionette/04_predictors/Han-BatTraits.csv') -> traits
traits$Pan <- gsub("_"," ",traits$Pan)

traitnames <- traits$Pan[!traits$Pan %in% cites$name]

citations2 <- c()

for(i in 1:length(traitnames)) {
  citations2[i] <- counter(traitnames[i])
  print(i)
}

cites2 <- data.frame(name = traitnames, cites = citations2)
cites <- rbind(cites, cites2)

write.csv(cites, 'C:/Users/cjcar/Documents/GitHub/virionette/04_predictors/Citations.csv')
