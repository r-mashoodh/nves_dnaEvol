## Looking for significant FET windows 

library(tidyverse)
library(patchwork)

fet <- data.table::fread("all_pops.idf.500bp.tef.fet")
head(fet)

# -log(product of p-value of all SNPs in window) for comparision of population 1 and 2

fet <- fet %>% 
  filter(V4 >= 0.8) %>% 
  mutate(V6 = as.numeric(substr(V6, 5,14)),
         V7 = as.numeric(substr(V7, 5,14)),
         V8 = as.numeric(substr(V8, 5,14)),
         V9 = as.numeric(substr(V9, 5,14)),
         V10 = as.numeric(substr(V10, 5,14)),
         V11 = as.numeric(substr(V11, 5,14))) %>% 
  rename(b_F1N1 = V6, 
         w_F1F2 = V7,
         b_F1N2 = V8,
         b_F2N1 = V9,
         w_N1N2 = V10, 
         b_F2N2 = V11) %>% 
  mutate(pos = paste(V1, V2, sep="_"),
         prod = b_F1N1*b_F2N2)


prd.fet <- quantile(fet$prod, probs = c(0.99)) #0.049

fet.rep1 <- quantile(fet$b_F1N1, probs = c(0.99)) #143

fet.rep2 <- quantile(fet$b_F2N2, probs = c(0.99)) #173

## look at differences between replicates
p1 <- ggplot(fet, aes(x=prod)) +
  geom_histogram() 
p2 <- ggplot(fet, aes(x=b_F1N1)) +
  geom_histogram() 
p3 <- ggplot(fet, aes(x=b_F2N2)) +
  geom_histogram()

(p1 / p2 / p3)

### get genes

## select Fst ratios that are at the top 0.5%

fet %>% 
  filter(prod > quantile(fet$prod, probs = c(0.995))) %>% 
  select(V1,V2,prod) %>% 
  mutate(start = V2 - 250,
         stop = V2 + 250) %>%
  select(V1,start,stop,prod) %>% 
  write.table(., file="prod_fet_of_interest.txt", quote=F, row.names=F, col.names=F, sep="\t")

system("bedtools intersect -a prod_fet_of_interest.txt -b Nves.features.sorted.bed -wao > prod_fet_intersect.bed")

#system("bedtools merge -i prod_fet_of_interest.txt > prod_fet_merged.bed")

wins <- data.table::fread("prod_fet_intersect.bed")

## here you have a bed file with all windows and what they intersect with in the annotation     
wins <- data.table::fread("prod_fet_intersect.bed")

### lets just look at lists of genes

## rename columns
col_names <- c("contig",
               "win_start",
               "win_stop", 
               "Fst_prod",
               "contig2",
               "feat_start",
               "feat_stop",
               "feature",
               "V9",
               "strand",
               "gene_name",
               "win_overlap")

colnames(wins) <- col_names

## lets look genes
fet_genes <- wins %>% 
  filter(feature!="RepeatModeler") %>% 
  filter(feature!=".") %>% 
  mutate(LOCid = substr(gene_name, 1, 12)) %>% 
  group_by(contig, gene_name, LOCid) %>% 
  summarise(avg_Fst_prod = mean(Fst_prod),
            num_windows = n(),
            overlap = sum(win_overlap))

length(unique(fet_genes$LOCid)) #648

fet_genes %>% 
  mutate(gene_name2 = str_replace(gene_name, ".*?;", "")) %>% #remove stuff before semicolon
  mutate(gene_name2 = str_replace(gene_name2, "isoform.*", "")) %>% # Remove 'isoform' and text after
  mutate(gene_name = ifelse(str_detect(LOCid, "Trna"), str_extract(LOCid, ".*?(?=;)"), gene_name2)) %>% #edit tRNA entries
  mutate(LOCid = gsub(";", "", LOCid)) %>% 
  mutate(gene_name = ifelse(gene_name2 == "gene_biotype=lncRNA", str_extract(gene_name2, "(?<=\\=).*"), gene_name)) %>% 
  select(-gene_name2) %>% 
  rename(geneID = "LOCid") %>% 
  select(contig, geneID, gene_name, Avg_score, num_windows, Category) %>% 
  write_tsv(., "TableS3.tsv")


### make a file for just 5'UTR genes
wins %>% 
  filter(feature == "5_prime_UTR") %>% 
  mutate(LOCid = substr(gene_name, 1, 12)) %>% 
  group_by(contig, gene_name, LOCid) %>% 
  summarise(avg_Fet_prod = mean(Fst_prod),
            num_windows = n()) %>% 
  mutate(gene_name2 = str_replace(gene_name, ".*?;", "")) %>% #remove stuff before semicolon
  mutate(gene_name2 = str_replace(gene_name2, "isoform.*", "")) %>% # Remove 'isoform' and text after
  mutate(gene_name = ifelse(str_detect(LOCid, "Trna"), str_extract(LOCid, ".*?(?=;)"), gene_name2)) %>% #edit tRNA entries
  mutate(LOCid = gsub(";", "", LOCid)) %>% 
  mutate(gene_name = ifelse(gene_name2 == "gene_biotype=lncRNA", str_extract(gene_name2, "(?<=\\=).*"), gene_name)) %>% 
  select(-gene_name2) %>% 
  write.csv(., "5_prime_UTR_genes.csv", quote=F, row.names = F)

#### Some plots ####

## -log(p) plot for Fig 2

ggplot(fet, aes(x=pos, y=prod)) +
  geom_point(color="black", alpha=0.4) +
  geom_hline(yintercept = 24500,
             color="red", linetype = "dashed") +
  #geom_hline(yintercept = quantile(dat400$ratio, probs = c(0.995)),
  #           color="blue", linetype = "dashed") +
  theme_classic() +
  #ggtitle("500kb windows") +
  xlab("Position") +
  ylab("-log10(p)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(
  "manhattan_ggsave.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 7,
  height = 2.25,
  units = "in",
  dpi = 300
)

### Fst plots - replicates separated 

maxdp400.tmp <- read.table("windows_500bp.maxDP600.idf.tef.fst")

# col1: reference contig (chromosome)
# col2: mean position of the sliding window
# col3: number of SNPs found in the window (not considering sites with a deletion) 
# col4: fraction of the window which has a sufficient coverage (min. coverage <= cov <= max. coverage) in every population;
# col5: average minimum coverage in all populations

# Cleanup data and compute ratio
dat400 <- maxdp400.tmp %>% 
  filter(V4 >= 0.80) %>% 
  mutate(V6 = as.numeric(substr(V6, 5,14)),
         V7 = as.numeric(substr(V7, 5,14)),
         V8 = as.numeric(substr(V8, 5,14)),
         V9 = as.numeric(substr(V9, 5,14)),
         V10 = as.numeric(substr(V10, 5,14)),
         V11 = as.numeric(substr(V11, 5,14))) %>% 
  rename(b_F1N1 = V6, 
         w_F1F2 = V7,
         b_F1N2 = V8,
         b_F2N1 = V9,
         w_N1N2 = V10, 
         b_F2N2 = V11) %>% 
  mutate(Fst_prod = b_F1N1*b_F2N2,
         pos = paste(V1, V2, sep="_")) %>% 
  rename(contig = "V1", 
         site = "V2")

rep1 <- ggplot(dat400, aes(x=pos, y=b_F1N1)) +
  geom_point(color="black", alpha=0.4) +
  #geom_hline(yintercept = quantile(dat400$ratio, probs = c(0.99)),
  #           color="red", linetype = "dashed") +
  #geom_hline(yintercept = quantile(dat400$ratio, probs = c(0.995)),
  #           color="blue", linetype = "dashed") +
  theme_classic() +
  ggtitle("FST FC1;NC1") +
  xlab("Position") +
  ylab("Fst prod") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(0,1)

rep2 <- ggplot(dat400, aes(x=pos, y=b_F2N2)) +
  geom_point(color="black", alpha=0.4) +
  #geom_hline(yintercept = quantile(dat400$ratio, probs = c(0.99)),
  #           color="red", linetype = "dashed") +
  #geom_hline(yintercept = quantile(dat400$ratio, probs = c(0.995)),
  #           color="blue", linetype = "dashed") +
  theme_classic() +
  ggtitle("FST FC2;NC2") +
  xlab("Position") +
  ylab("Fst prod") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(0,1)


(rep1 / rep2)

ggsave(
  "fst_reps_ggsave.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)
