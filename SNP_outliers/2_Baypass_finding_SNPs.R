#### Script to process baypass output

library(tidyverse)
library(poolfstat)

## Import sync file and extract data using poolfstat
psizes <- c(41,52,52,59)
psizes.haploid <- psizes*2
pnames <- c("F1","N1","F2","N2")

summary(gene.fst$V5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 20.00   49.00   58.00   61.94   66.00  389.00

evol_pops <- popsync2pooldata(sync.file = "all_pops.idf.tef.sync",
                              poolsizes = psizes.haploid,
                              poolnames = pnames,
                              min.rc = 4,
                              min.cov.per.pool = 30,
                              max.cov.per.pool = 700,
                              min.maf = 0.05,
                              noindel = TRUE,
                              nlines.per.readblock = 2e+06,
                              nthreads=6)

#write_rds(evol_pops, "evol_pops700.rds", compress="gz")

# write SNP data into baypass format
pooldata2genobaypass(evol_pops, writing.dir = getwd(), subsamplesize = -1)

# rename files
system("mv genobaypass Evol700.genobaypass")
system("mv snpdet Evol700.snpdet")
system("mv poolsize Evol700.snpdet")

# run this one first to get covariance  matrix etc
# baypass_2.3/sources/i_baypass -npop 4 -gfile Evol700.genobaypass -poolsizefile Evol700.poolsize -d0yij 8 -outprefix Evol700 -nthreads 64

# then run this to look at group differences with pops as covariates (0 vs 1)
#baypass_2.3/sources/i_baypass  -nthreads 64 -npop 4 -gfile Evol700.genobaypass -efile groups -poolsizefile Evol700.poolsize -d0yij 8 -auxmodel -omegafile Evol700_mat_omega.out -outprefix Evol700.grps.aux

# http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.3.pdf

## extract relationships between populations from covariance matrix
source("baypass_utils.R")
require(corrplot) ; require(ape)

omega=as.matrix(read.table("baypass_out/Evol700_mat_omega.out"))

pop.names <- c("FC1","NC1","FC2","NC2")

colnames(omega) <- pop.names
rownames(omega) <- pop.names

#Omega output as heatmap
cor.mat=cov2cor(omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))

#Omega output as tree
tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

#### identify outliers associated with environmental variables.
### (i.e. SNPs with dB > 20, from the dB column of the *_summary_betai.out output file).

# Gautier M. Genome-wide scan for adaptive divergence and association with population-specific covariates.
# Genetics. 2015;4:1555–79

# Jefferys H. Theory of probability (3rd edition). New York: Oxford university press; 1961

# the estimated Bayes Factor (column BF(dB)) in dB units (i.e., 10×log10(BF)) comparing the models
# with (βi ̸= 0) and without (βi = 0) association of the SNP with the given covariable.

## Get SNP locations
SNPs <- data.table::fread("baypass_out/Evol700.snpdet")
SNPs$MRK <- c(1:nrow(SNPs))

head(SNPs)


## Get BFs and filter
BayesFactors <- data.table::fread("baypass_out/Evol700.grps.aux_summary_betai.out.gz")

## Get SNP locations
SNPs <- data.table::fread("baypass_out/Evol700.snpdet")
SNPs$MRK <- c(1:nrow(SNPs))

sigSNPs <- SNPs %>%
  left_join(BayesFactors %>%
              filter(`BF(dB)` > 21)) %>%
  drop_na() %>%
  mutate(stop = V2+1) %>%
  rename(seqnames = V1,
         start = V2,
         ref = V3,
         alt = V4) %>%
  select(seqnames, start, stop, MRK, ref, alt, `BF(dB)`)

write.table(sigSNPs, "Evol.SigSNPs_Baypass.txt", quote=F, row.names = F, col.names = F, sep="\t")

### Get SNP frequency differences
## combine SNP ids and freq data
raw <- cbind(SNPs, data.table::fread("baypass_out/Evol700.genobaypass"))
colnames(raw) <- c("contig", "pos", "ref", "alt", "MRK", "F1.ref", "F1.alt", "N1.ref", "N1.alt",
                    "F2.ref", "F2.alt", "N2.ref", "N2.alt")
freq <- raw %>%
  mutate(F1.tot = F1.ref + F1.alt,
         N1.tot = N1.ref + N1.alt,
         F2.tot = F2.ref + F2.alt,
         N2.tot = N2.ref + N2.alt) %>%
  mutate(F1 = round(F1.ref/F1.tot, 2),
         N1 = round(N1.ref/N1.tot, 2),
         F2 = round(F2.ref/F2.tot, 2),
         N2 = round(N2.ref/N2.tot, 2)) %>%
  select(contig, pos, MRK, F1, N1, F2, N2) %>%
  filter(MRK %in% sigSNPs$MRK)

head(freq)

### Where are these SNPs?

#system("bedtools intersect -a Evol.SigSNPs_Baypass.txt -b ../../Gen30_PoolSeq/Popoolation2_data/Nves.features.sorted.bed -wao > Evol.SigSNPs_Baypass_Intersect.bed")

system("bedtools window -a Evol.SigSNPs_Baypass.txt -b genes.bed -w 500 > Evol.SigSNPs_Baypass_Intersect.bed")

nrow(sigSNPs) #3086 SNPs

snp.intersect <- data.table::fread("Evol.SigSNPs_Baypass_Intersect.bed")

colnames(snp.intersect)[1] <- c("contig")
colnames(snp.intersect)[2] <- c("pos")
colnames(snp.intersect)[4] <- c("MRK")

length(unique(snp.intersect$MRK)) #1644 of these were within 500bp of gene coordinates

tableS7 <- sigSNPs %>%
  select(seqnames,start,MRK) %>%
  rename(contig = seqnames,
         pos = start) %>%
  left_join(freq) %>%
  left_join(snp.intersect %>%
              select(contig, pos, MRK, V14)) %>%
  rename(gene_name = V14) %>%
    mutate(gene_name = str_replace(gene_name, "isoform.*", "")) %>% # Remove 'isoform' and text after
  group_by(contig, pos, MRK, F1, N1, F2, N2) %>%
  summarize(gene_name = paste(gene_name, collapse = ","))

write_tsv(tableS7, "TableS7_SNP_freq.tsv")



## What is within 500bp of SNPs?
bp_hits <- data.table::fread("Evol.SigSNPs_Baypass_Intersect.bed") %>%
  filter(V11=="gene") %>%
  select(V1, V2) %>%
  distinct()

nrow(bp_hits) #only 1644 SNPs are within 500bp of genes
nrow(sigSNPs) - nrow(bp_hits) #1442 SNPs near nothing

### What genes are near these SNPs?
bp_genes <- data.table::fread("Evol.SigSNPs_Baypass_Intersect.bed") %>%
  filter(V11=="gene") %>%
  #filter(V10!=".") %>%
  mutate(LOCid = substr(V14, 1, 12)) %>%
  group_by(V14,LOCid) %>%
  summarise(n=n(),
            BFmean=mean(V7))

length(unique(bp_genes$LOCid)) #1150

## How much overlap with other methods?
## top 0.5% of FET windows
xl <- readxl::read_xlsx("fet_genes.xlsx")

# get overlap
overlap.05<- bp_genes %>%
  filter(LOCid %in% xl$LOCid)

# count overlaps
length(unique(overlap.05$LOCid)) # 220 if you use the top 0.05% of FET windows

overlap.05 %>%
  rename(gene_name = V14,
         SNPs = n) %>%
  select(LOCid, gene_name, SNPs, BFmean) %>%
    mutate(gene_name2 = str_replace(gene_name, ".*?;", "")) %>% #remove stuff before semicolon
    mutate(gene_name2 = str_replace(gene_name2, "isoform.*", "")) %>% # Remove 'isoform' and text after
    mutate(gene_name = ifelse(str_detect(LOCid, "Trna"), str_extract(LOCid, ".*?(?=;)"), gene_name2)) %>% #edit tRNA entries
    mutate(LOCid = gsub(";", "", LOCid)) %>%
    mutate(gene_name = ifelse(gene_name2 == "gene_biotype=lncRNA", str_extract(gene_name2, "(?<=\\=).*"), gene_name)) %>%
    select(-gene_name2) %>%
  write.csv(., "TableS8_overlap.csv", quote=F, row.names = F)

### Lets look at how many of these SNPs occur in windows that were significant (FET windows)
## prod_fet_merged.bed <- merged intervals

system("bedtools intersect -a prod_fet_merged.bed -b Evol.SigSNPs_Baypass.txt -wao > SNPs_in_fet_windows.bed")

fet_window_SNPs <- read.table("SNPs_in_fet_windows.bed")

head(fet_window_SNPs)

## 138 out of 1576 were directly in the windows
table(fet_window_SNPs$V11)
# 0    1
# 1576  138
