### Genome wide Tajima's D

F1.td <- read.table("F1.500.ovlp.tef.td", stringsAsFactors = F, na.strings = "na")
F1.td$pop <-c("F1")
N1.td <- read.table("N1.500.ovlp.tef.td", stringsAsFactors = F, na.strings = "na")
N1.td$pop <-c("N1")
F2.td <- read.table("F2.500.ovlp.tef.td", stringsAsFactors = F, na.strings = "na")
F2.td$pop <-c("F2")
N2.td <- read.table("N2.500.ovlp.tef.td", stringsAsFactors = F, na.strings = "na")
N2.td$pop <-c("N2")

td.tmp <- rbind(F1.td, F2.td, N1.td, N2.td) 

td.df <- td.tmp %>% 
  filter(V4 >= 0.5) %>% 
  drop_na() %>% 
  mutate(region = paste(V1,V2,sep="_")) %>% 
  select(V5,pop,region) %>% 
  pivot_wider(., values_from=V5, names_from=pop) %>% 
  drop_na() %>% 
  pivot_longer(2:5, values_to = "td", names_to="Pop") %>% 
  mutate(expEvo = substr(Pop, 1,1),
         Block = substr(Pop, 2,2))

# kruskal.test(td ~ expEvo, data = subset(td.df, Block == "1"))
# kruskal.test(td ~ expEvo, data = subset(td.df, Block == "2"))

## looks pretty normally distributed ##
## run a t-test?
td.df %>% 
  ggplot(., aes(x=td, fill=Pop)) +
  geom_histogram() +
  facet_wrap(.~Pop)

## mean ###
calculate_mean <- function(data) {
  mean_val <- mean(data)
  boot_obj <- boot(data, statistic = function(d, i) mean(d[i]), R = 1000)
  ci <- boot.ci(boot_obj, type = "basic", conf = 0.95)
  list(mean = mean_val, lower_ci = ci$basic[4], upper_ci = ci$basic[5])
}

td.ci <- td.df %>% 
  group_by(Pop) %>% 
  do(ci = calculate_mean(.$td)) %>% 
  unnest_wider(ci) %>% 
  mutate(expEvo = substr(Pop, 1,1),
         Block = substr(Pop, 2,2))

td.ci %>% 
  ggplot(., aes(x=expEvo, y=mean, color=Block, group=Block)) +
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = lower_ci, ymax = upper_ci)) +
  ylab("tajimas D (mean)") +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  ylim(-0.2,0.03)

t.test(td ~ expEvo, data = subset(td.df, Block == "1"))
t.test(td ~ expEvo, data = subset(td.df, Block == "2"))

### What about diverged genes/loci
## make a bed file with window coordinates and scores for each pop
td.wide <- td.df %>%
  pivot_wider(., id_cols=region, values_from=td, names_from=c("Pop")) %>%
  separate(., col = region, into = c("contig", "pos"), sep = "_(?=\\w+$)", remove = T) %>% 
  mutate(F1xF2 = ifelse(F1 < 0 & F2 < 0, -F1 * F2, F1 * F2),
         N1xN2 = ifelse(N1 < 0 & N2 < 0, -N1 * N2, N1 * N2))

td.wide %>%
  mutate(start = as.numeric(pos),
         stop = as.numeric(pos) + 1) %>%
  select(contig, start, stop, F1, F2, N1, N2, F1xF2, N1xN2) %>%
  write.table("td_scores.bed", quote=F, row.names=F, col.names=F, sep="\t")

## create a file that has every gene +/- 5kb of Td
system("bedtools window -a genes.bed -b td_scores.bed -w 5000 > gene_td.txt")

gene_td <- data.table::fread("gene_td.txt") %>%
  mutate(LOCid = substr(V7, 1,12))

glimpse(gene_td)

colnames(gene_td)[1:16] <- c("contig_ref", "gene_start", "gene_stop", "ref", "score", "strand",
                             "gene_name", "contig", "pos", "pos_1", "F1", "F2", "N1", "N2", "F1xF2",
                             "N1xN2")


### all intervals close to selected genes

fet_td <- gene_td %>% 
  filter(LOCid %in% fet_genes$LOCid) %>% 
  select(LOCid, F1, F2, N1, N2) %>% 
  pivot_longer(., cols = 2:5, names_to="pop", values_to="Td") %>% 
  mutate(Block = substr(pop, 2, 2),
         expEvo = substr(pop, 1, 1))

t.test(Td ~ expEvo, data = subset(fet_td, Block == "1"))
t.test(Td ~ expEvo, data = subset(fet_td, Block == "2"))

fet_td %>% 
  ggplot(., aes(x=Td)) +
  geom_histogram() +
  facet_wrap(.~pop)
