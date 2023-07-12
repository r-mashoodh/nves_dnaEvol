## Analysis of genome-wide differences in Theta

library(tidyverse)
library(boot)
library(ggridges)
library(patchwork)

## Read in data
F1.theta <- read.table("F1.1000.tef.theta", stringsAsFactors = F, na.strings = "na")
F1.theta$pop <-c("F1")
N1.theta <- read.table("N1.1000.tef.theta", stringsAsFactors = F, na.strings = "na")
N1.theta$pop <-c("N1")
F2.theta <- read.table("F2.1000.tef.theta", stringsAsFactors = F, na.strings = "na")
F2.theta$pop <-c("F2")
N2.theta <- read.table("N2.1000.tef.theta", stringsAsFactors = F, na.strings = "na")
N2.theta$pop <-c("N2")

## Clean up data
theta.df <- rbind(F1.theta, F2.theta, N1.theta, N2.theta) %>% 
  filter(V4 > 0.5) %>% 
  drop_na() %>% 
  mutate(region = paste(V1,V2,sep="_")) %>% 
  select(V5,pop,region) %>% 
  pivot_wider(., values_from=V5, names_from=pop) %>% 
  drop_na() %>% 
  pivot_longer(2:5, values_to = "theta", names_to="Pop") %>% 
  mutate(expEvo = substr(Pop, 1,1),
         Block = substr(Pop, 2,2)) %>% 
  mutate(BlockNum = ifelse(Block == 1, -1, 1),
         expEvoNum = ifelse(expEvo == "F", -1, 1))

## non-parametric tests of significance
kruskal.test(theta ~ expEvo, data = subset(theta.df, Block == "1"))
kruskal.test(theta ~ expEvo, data = subset(theta.df, Block == "2"))

## Calculate median and bootstrapped confidence intervals 

# calculate_ci <- function(data) {
#   median_val <- median(data)
#   boot_obj <- boot(data, statistic = function(d, i) median(d[i]), R = 1000)
#   ci <- boot.ci(boot_obj, type = "basic", conf = 0.95)
#   list(median = median_val, lower_ci = ci$basic[4], upper_ci = ci$basic[5])
# }

theta.ci <- theta.df %>% 
  group_by(Pop) %>% 
  do(ci = calculate_ci(.$theta)) %>% 
  unnest_wider(ci) %>% 
  mutate(expEvo = substr(Pop, 1,1),
         Block = substr(Pop, 2,2))


### plot median and its confidence intervals 
y <- theta.ci %>% 
  ggplot(., aes(x=expEvo, y=median, color=expEvo, group=Block)) +
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = lower_ci, ymax = upper_ci)) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  ylim(0.0072,0.0077) +
  theme(legend.position = "none")

plot <- (x|y)

ggsave(plot,filename = "diversity.eps", device="eps", width=5, height=2.2, units="in")


## density plot
hist_t <- ggplot(theta.df, aes(y=expEvo, x=theta, fill=expEvo)) +
  ggridges::geom_density_ridges2(alpha=0.5) +
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.position = "none")

ggsave(hist_t,filename = "hist_t.eps", device="eps", width=1.75, height=1, units="in")