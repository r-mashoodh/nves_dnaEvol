## Testing GO terms for divergent genes 

library(topGO)

## read in diverged gene list
xl <- readxl::read_xlsx("fet_genes.xlsx")

## get genes of interest in a list
fet_genes <- unique(xl$LOCid)

## Go term file
## new GO annotation, supplemented by orthofinder
geneID2GOdm <- readMappings("fly_orthoF_merged.go.tab") # mf_crowdGO.tab or cc_crowdGO.tab; remake files above

geneUniversedm <- names(geneID2GOdm)
genesOfInterest <- as.character(fet_genes)

## identify where your genes of interest are within the gene Universe
geneList <- factor(as.integer(geneUniversedm %in% genesOfInterest))
names(geneList) <- geneUniversedm

myGOdataBP<- new("topGOdata", 
                 description="Evol_Me", 
                 ontology="BP", 
                 allGenes=geneList,
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GOdm)

resultFisherBP <- runTest(myGOdataBP, algorithm="weight01", statistic="fisher")


allResBP <- GenTable(myGOdataBP, 
                     classicFisher = resultFisherBP,
                     ranksOf = "classicFisher", 
                     topNodes = 200)

allResBP<- allResBP %>% 
  filter(classicFisher < 0.05)

write.csv(allResBP, "GO_terms.csv")

## go terms of interest
row_numbers <- c(1, 12, 28, 31, 36, 47, 51, 84, 91, 131)

library(tidyverse)

allResBP %>% 
  #rownames_to_column(var="rowname") %>% 
  slice(row_numbers) %>% 
  #mutate(Term = paste(Term, rowname, sep="_")) %>% 
  mutate(p = as.numeric(classicFisher)) %>%
  mutate(name = fct_reorder(Term, desc(p))) %>%
  ggplot(., aes(x=name, y=-log(p))) +
  geom_bar(stat="identity", width=0.7, fill="#166F64", alpha=0.7) +
  theme_classic() +
  coord_flip()
