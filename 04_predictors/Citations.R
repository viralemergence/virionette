install.packages("easyPubMed")
library("easyPubMed")

Han <- read.csv('~/Github/cleanbats_betacov/clean data/Han-BatTraits_compatible.csv')

counter <- function(name) {as.numeric(as.character(get_pubmed_ids(name)$Count))}

Han$citations <- 0

for(i in 1:nrow(Han)) {
  Han$citations[i] <- counter(Han$Pan[i])
  print(i)
}

Han <- Han[,c('Pan','citations')]
write.csv(Han, 'Citations.csv')
