## piN/piS

setwd("/mnt/beegfs/home5/zoology/rm786/beetle_PoolSeq/intrapop/tef_out")

# Within species you can compute the Pn/Ps (relative abundance of non-synonymous and 
#                                           synonymous polymorphisms) that measure the 
# direct effect of natural selection removing slighlty deleterious non-synonymous variants.
# 
# http://www.plospathogens.org/article/info%3Adoi%2F10.1371%2Fjournal.ppat.1002893

# col1: gene name
# col2: non-synonymous length
# col3: synonymous length
# col4: non-synonymous SNPs
# col5: synonymous SNPs
# col6: non-synonymous measure (pi/theta/d)
# col7: synonyomous measure (pi/theta/d)

col_syn <- c("gene", "N_length", "S_length", "N_SNPs", "S_SNPs", "PiN", "PiS")

## Read data in

F1.syn <- data.table::fread("F1.syn.tef.pi", na.strings = "na")
colnames(F1.syn) <- col_syn
F1.syn$pop <- c("F1")

F2.syn <- data.table::fread("F2.syn.tef.pi", na.strings = "na")
colnames(F2.syn) <- col_syn
F2.syn$pop <- c("F2")

N1.syn <- data.table::fread("N1.syn.tef.pi", na.strings = "na")
colnames(N1.syn) <- col_syn
N1.syn$pop <- c("N1")

N2.syn <- data.table::fread("N2.syn.tef.pi", na.strings = "na")
colnames(N2.syn) <- col_syn
N2.syn$pop <- c("N2")

pNS <- rbind(F1.syn,F2.syn,N1.syn,N2.syn)
pNS$expEvo <- substr(pNS$pop, 1,1)
pNS$Block <- substr(pNS$pop, 2,2)

## get some descriptives
pNS %>% 
  group_by(pop, expEvo, Block) %>% 
  summarise(n=n(),
            median_piN = median(PiN),
            median_piS = median(PiS))

## piN/piS
## contains diverged genes
xl <- readxl::read_xlsx("../../map/fet_genes_regulatory.xlsx")

## calculate piN/piS
pNS$pNS <- pNS$PiN / pNS$PiS

## get descriptives
pNS %>% 
  filter(!is.infinite(pNS)) %>% 
  #mutate(fet = ifelse(gene %in% xl$LOCid, 1, 0)) %>% 
  group_by(pop, expEvo, Block) %>% 
  summarise(n=n(),
            median_pNS = median(pNS,na.rm=T))

## test of piN/piS for all genes
kruskal.test(pNS ~ expEvo, data = subset(pNS, Block == "1"))
kruskal.test(pNS ~ expEvo, data = subset(pNS, Block == "2"))

## test of piN/piS for diverged genes
pNS %>% 
  filter(!is.infinite(pNS)) %>% 
  mutate(fet = ifelse(gene %in% xl$LOCid, 1, 0)) %>% 
  filter(fet == 1) %>% 
  group_by(pop, expEvo, Block) %>% 
  summarise(n=n(),
            median_pNS = median(pNS,na.rm=T))

pNS.fet <- pNS %>% 
  filter(!is.infinite(pNS)) %>% 
  mutate(fet = ifelse(gene %in% xl$LOCid, 1, 0)) %>% 
  filter(fet == 1)

kruskal.test(pNS ~ expEvo, data = subset(pNS.fet, Block == "1"))
kruskal.test(pNS ~ expEvo, data = subset(pNS.fet, Block == "2"))

## do FC harbour more variants at diverged genes?

library(boot)

# Set the number of bootstrap iterations
num_iterations <- 1000
sample_size <- 600

# Clean up df for bootstrapping
pNS.clean <- pNS %>% drop_na(pNS) %>% filter(!is.infinite(pNS))

# Get unique populations
populations <- unique(pNS.clean$pop)

# Create an empty data frame to store the bootstrap medians
bootstrap_medians_df <- data.frame(matrix(NA, nrow = num_iterations, ncol = length(populations)))
colnames(bootstrap_medians_df) <- populations

# Perform bootstrap resampling and calculate the median for each population
for (i in 1:num_iterations) {
  bootstrap_sample <- lapply(populations, function(pop) {
    subset_data <- subset(pNS.clean, pop == pop)$pNS
    sample_data <- sample(subset_data, replace = F, size = sample_size)
    median(sample_data)
  })
  bootstrap_medians_df[i, ] <- unlist(bootstrap_sample)
}

# Print the results
print(bootstrap_medians_df)


## get the values for selected genes
pop_medians <- pNS %>% 
  filter(!is.infinite(pNS)) %>% 
  mutate(fet = ifelse(gene %in% xl$LOCid, 1, 0)) %>% 
  filter(fet == 1) %>% 
  group_by(pop, expEvo, Block) %>% 
  summarise(n=n(),
            median_pNS = median(pNS,na.rm=T)) %>% 
  pull(median_pNS)

pop_medians
# [1] 0.02437635 0.02294070 0.00000000 0.00000000

## FC variation is greater than average sampling of genes?
length(which(bootstrap_medians_df$F1 > pop_medians[1])) / 1001 #0.01
length(which(bootstrap_medians_df$F2 > pop_medians[2])) / 1001 #0.0006

## NC is less than average sampling of genes?
length(which(bootstrap_medians_df$N1 < pop_medians[3])) / 1001 #0.00
length(which(bootstrap_medians_df$N2 < pop_medians[4])) / 1001 #0.00

## So what this suggests is that FC maintained more variation on average on those genes that were
## "selected", and likewise, NC reduced variation more though to different degrees between the NC replicates
## might reflect the slight difference in the strength of response in NC1 that appears throughout the measures

FC1 <- ggplot(bootstrap_medians_df, aes(x=F1)) +
  geom_histogram() +
  geom_vline(xintercept=pop_medians[1], color="red") +
  ggtitle("F1, p=0.145")

FC2 <- ggplot(bootstrap_medians_df, aes(x=F2)) +
  geom_histogram() +
  geom_vline(xintercept=pop_medians[2], color="red") +
  ggtitle("F2, p=0.211")

NC1 <- ggplot(bootstrap_medians_df, aes(x=N1)) +
  geom_histogram() +
  geom_vline(xintercept=pop_medians[3], color="red") +
  ggtitle("N1, p=0.0")

NC2 <- ggplot(bootstrap_medians_df, aes(x=N2)) +
  geom_histogram() +
  geom_vline(xintercept=pop_medians[4], color="red") +
  ggtitle("N2, p=0.0")

(FC1 / FC2 | NC1 / NC2)

### calculate bootstrap medians for diverged and not diverged genes

calculate_median <- function(data) {
  median_val <- median(data)
  boot_obj <- boot(data, statistic = function(d, i) median(d[i]), R = 1000)
  ci <- boot.ci(boot_obj, type = "basic", conf = 0.95)
  list(median = median_val, lower_ci = ci$basic[4], upper_ci = ci$basic[5])
}

ci_vals <- pNS.clean %>% 
  mutate(fet = ifelse(gene %in% xl$LOCid, 1, 0)) %>% 
  group_by(pop, fet) %>% 
  do(ci = calculate_median(.$pNS)) %>% 
  unnest_wider(ci) 

ci_vals$expEvo <- paste(substr(ci_vals$pop, 1, 1), "C", sep="")

ggplot(ci_vals, aes(x=pop, y=median, color=as.factor(fet))) +
  geom_point(geom_jitter=0.8) +
  geom_linerange(aes(ymin = lower_ci, ymax= upper_ci)) +
  theme_classic()

ggplot(ci_vals, aes(x=pop, y=median, color=expEvo)) +
  geom_point() +
  geom_linerange(aes(ymin = lower_ci, ymax= upper_ci)) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  facet_wrap(.~fet)