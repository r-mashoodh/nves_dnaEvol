## Analysis of genome-wide differences in Pi

library(tidyverse)
library(boot)
library(ggridges)

## Read in Pi values
F1.pi <- read.table("F1.1000.tef.pi", stringsAsFactors = F, na.strings = "na")
F1.pi$pop <-c("F1")
N1.pi <- read.table("N1.1000.tef.pi", stringsAsFactors = F, na.strings = "na")
N1.pi$pop <-c("N1")
F2.pi <- read.table("F2.1000.tef.pi", stringsAsFactors = F, na.strings = "na")
F2.pi$pop <-c("F2")
N2.pi <- read.table("N2.1000.tef.pi", stringsAsFactors = F, na.strings = "na")
N2.pi$pop <-c("N2")

pi.df <- rbind(F1.pi, F2.pi, N1.pi, N2.pi) %>% 
  filter(V4 > 0.5) %>% 
  drop_na() %>% 
  mutate(region = paste(V1,V2,sep="_")) %>% 
  select(V5,pop,region) %>% 
  pivot_wider(., values_from=V5, names_from=pop) %>% 
  drop_na() %>% 
  pivot_longer(2:5, values_to = "pi", names_to="Pop") %>% 
  mutate(expEvo = substr(Pop, 1,1),
         Block = substr(Pop, 2,2)) %>% 
  mutate(BlockNum = ifelse(Block == 1, -1, 1),
         expEvoNum = ifelse(expEvo == "F", -1, 1))

### histogram of values
pi.df %>% 
  ggplot(., aes(x=pi, fill=Pop)) +
  geom_histogram() +
  facet_wrap(.~Pop)

## non parametric tests of significance
kruskal.test(pi ~ expEvo, data = subset(pi.df, Block == "1"))
kruskal.test(pi ~ expEvo, data = subset(pi.df, Block == "2"))

### calculate median and bootstrapped confidence intervals
calculate_ci <- function(data) {
  median_val <- median(data)
  boot_obj <- boot(data, statistic = function(d, i) median(d[i]), R = 1000)
  ci <- boot.ci(boot_obj, type = "basic", conf = 0.95)
  list(median = median_val, lower_ci = ci$basic[4], upper_ci = ci$basic[5])
}

pi.ci <- pi.df %>% 
  group_by(Pop) %>% 
  do(ci = calculate_ci(.$pi)) %>% 
  unnest_wider(ci) %>% 
  mutate(expEvo = substr(Pop, 1,1),
         Block = substr(Pop, 2,2))

pi.ci

### plot median and its confidence intervals 
x <- pi.ci %>% 
  ggplot(., aes(x=expEvo, y=median, color=expEvo, group=Block)) +
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = lower_ci, ymax = upper_ci)) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none") +
  ylim(0.0074, 0.0078)

(x|y)

## density plot
hist_p <- ggplot(pi.df, aes(y=expEvo, x=pi, fill=expEvo)) +
  ggridges::geom_density_ridges2(alpha=0.5) +
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.position = "none")

ggsave(hist_p,filename = "hist_p.eps", device="eps", width=1.75, height=1, units="in")