install.packages("easyPubMed")
library("easyPubMed")

tree <- readRDS("C:/Users/cjcar/Documents/GitHub/virionette/04_predictors/Full Supertree.rds")

counter <- function(name) {as.numeric(as.character(get_pubmed_ids(name)$Count))}

citations <- c()

for(i in 1:length(tree$tip.label)) {
  citations[i] <- counter(tree$tip.label[i])
  print(i)
}

cites <- data.frame(name = tree$tip.label, cites = citations)
write.csv(cites, 'C:/Users/cjcar/Documents/GitHub/virionette/04_predictors/Citations.csv')
