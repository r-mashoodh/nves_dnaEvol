## piN/piS

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
