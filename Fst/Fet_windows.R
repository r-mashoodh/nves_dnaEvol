## FET
setwd("/mnt/beegfs/home5/zoology/rm786/beetle_PoolSeq/map")

plot(dat400$prod, fet$prod)

nrow(dat400)
nrow(fet)

setwd("/mnt/beegfs/home5/zoology/rm786/beetle_PoolSeq/map")
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

#fet$padj1 <- p.adjust(exp(-fet$b_F1N1), method = "bonferroni", n=length(fet$b_F1N1))
#fet$padj2 <- p.adjust(exp(-fet$b_F2N2), method = "bonferroni", n=length(fet$b_F2N2))

hist(fet$prod)
prd.fet <- quantile(fet$prod, probs = c(0.99)) #0.049

fet.rep1 <- quantile(fet$b_F1N1, probs = c(0.99)) #143

fet.rep2 <- quantile(fet$b_F2N2, probs = c(0.99)) #173

## more extremes in block 2, but this gets better as you relax coverage
sum((fet$b_F1N1 > .3))
sum((fet$b_F2N2 > .3))

p1 <- ggplot(fet, aes(x=prod)) +
  geom_histogram() 
p2 <- ggplot(fet, aes(x=b_F1N1)) +
  geom_histogram() 
p3 <- ggplot(fet, aes(x=b_F2N2)) +
  geom_histogram()

library(patchwork)

(p1 / p2 / p3)

### get genes

## select Fst ratios that are at the top 1%

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

## here you have a bed file with all windows that are at the top 1&       
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

## lets look at windows

tmp <- wins %>% 
  filter(feature!="RepeatModeler") %>% 
  mutate(LOCid = substr(gene_name, 1, 12)) %>% 
  group_by(feature, gene_name, LOCid) %>% 
  summarise(avg_Fst_prod = mean(Fst_prod),
            num_windows = n())

xl <- readxl::read_xlsx("fet_genes_regulatory.xlsx")
xl <- xl %>% 
  mutate(Category = case_when(
    Category %in% c("gene expression", "chromatin") ~ "gene expression, chromatin",
    Category %in% c("txf", "transcription factor") ~ "transcription factor",
    is.na(Category) ~ "other protein coding",  # Replace NAs with "other protein coding"
    TRUE ~ Category  # Keep other categories unchanged
  )) %>% 
  select(LOCid, Category)

tmp %>% 
  group_by(feature) %>% 
  summarise(n=(sum(num_windows)/sum(tmp$num_windows)*100))

tx_win <- tmp %>% 
  left_join(xl) %>% 
  mutate(new = ifelse(feature == "gene", Category, feature)) %>% 
  group_by(new) %>% 
  summarise(n=(sum(num_windows)/sum(tmp$num_windows)*100),
            n_genes = n()) %>% 
  mutate(new = ifelse(new == ".", "no overlap", new))

ggplot(tx_win, aes(x = "", y = n, fill = new)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  #geom_text(aes(y = round(lab.ypos,2), label = n), color = "white")+
  scale_fill_brewer(palette="Set1") +
  theme_void()

## what is only hitting 5' UTR
## only hitting UTRs
x <- tmp %>% 
  left_join(xl) %>% 
  mutate(new = ifelse(feature == "gene", Category, feature)) %>% filter(new == "5_prime_UTR")
# only hitting genes
y <- tmp %>% 
  left_join(xl) %>% 
  mutate(new = ifelse(feature == "gene", Category, feature)) %>% filter(new != "5_prime_UTR")
# unique
z <- x %>% filter(!LOCid %in% y$LOCid)

fet_genes <- wins %>% 
  filter(feature!="RepeatModeler") %>% 
  filter(feature!=".") %>% 
  mutate(LOCid = substr(gene_name, 1, 12)) %>% 
  group_by(contig, gene_name, LOCid) %>% 
  summarise(avg_Fst_prod = mean(Fst_prod),
            num_windows = n(),
            overlap = sum(win_overlap))

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

####

## What if we increased to 1% how many more genes

fet %>% 
  filter(prod > quantile(fet$prod, probs = c(0.99))) %>% 
  select(V1,V2,prod) %>% 
  mutate(start = V2 - 250,
         stop = V2 + 250) %>%
  select(V1,start,stop,prod) %>% 
  write.table(., file="1perc_prod_fet_of_interest.txt", quote=F, row.names=F, col.names=F, sep="\t")
system("bedtools intersect -a 1perc_prod_fet_of_interest.txt -b Nves.features.sorted.bed -wao > 1perc_prod_fet_intersect.bed")

wins2 <- data.table::fread("1perc_prod_fet_intersect.bed") %>% 
  filter(V8!="RepeatModeler")

table(wins2$V8)

length(unique(wins2$V11)) ##1137

## plot each rep

F1N1.m <- fet %>% 
  filter(V4 >= 0.75) %>% 
  ggplot(., aes(x=pos, y=b_F1N1)) +
  geom_point(color="black", alpha=0.4) +
  geom_hline(yintercept = quantile(fet$b_F1N1, probs = c(0.995)),
             color="red", linetype = "dashed") +
  # geom_hline(yintercept = quantile(dat400$ratio, probs = c(0.995)),
  #           color="blue", linetype = "dashed") +
  theme_classic() +
  ggtitle("F1;N1") +
  xlab("Position") +
  ylab("-log10(p)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

F2N2.m <- fet %>% 
  filter(V4 >= 0.75) %>% 
  ggplot(., aes(x=pos, y=b_F2N2)) +
  geom_point(color="black", alpha=0.4) +
  geom_hline(yintercept = quantile(fet$b_F2N2, probs = c(0.995)),
             color="red", linetype = "dashed") +
  # geom_hline(yintercept = quantile(dat400$ratio, probs = c(0.995)),
  #           color="blue", linetype = "dashed") +
  theme_classic() +
  ggtitle("F2;N2") +
  xlab("Position") +
  ylab("-log10(p)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

(F1N1.m / F2N2.m)






########################################################################
fet %>% 
  mutate(prod = b_F1N1*b_F2N2) %>% 
  arrange(-prod) %>% 
  head()

dat400 %>% 
  filter(V4 >= 0.8) %>% 
  arrange(-prod) %>% 
  head()

tmp<-dat400 %>% 
  filter(prod > prd)
fet %>% 
  mutate(prod = b_F1N1*b_F2N2) %>% 
  filter(prod > quantile(prod, probs=c(0.99))) %>% 
  filter(pos %in% tmp$pos)

hist(dat400$prod)
prd <- quantile(dat400$prod, probs = c(0.99)) #0.049

p1 <- ggplot(fet, aes(x=b_F1N1*b_F2N2)) +
  geom_histogram() 
p2 <- ggplot(fet, aes(x=b_F1N1)) +
  geom_histogram()
p3 <- ggplot(fet, aes(x=b_F2N2)) +
  geom_histogram()

(p1 / p2 / p3)

# yintercept = quantile(dat400$ratio, probs = c(0.995))
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