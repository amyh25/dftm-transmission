require(tidyverse)
require(ggpubr)
require(binom)
theme_set(theme_pubr())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv(file="manual_avoid_effect.csv",
                 header=TRUE, sep=",", stringsAsFactors=FALSE)
data$morphotype <- factor(data$morphotype, levels = c("SNPV","MNPV"))
data$tree_sp <- factor(data$tree_sp, levels = c("GR","DO"))
data <- data %>% filter(control != "NA", treatment != "NA")
data$control[data$control < 0] <- 0
data$treatment[data$treatment < 0] <- 0
data <- data %>% cbind(D = data$control - data$treatment)

data <- data %>% cbind(list(D_hat = data$D))
data$D_hat[data$tree_sp=="DO"] <-
  data$D[data$tree_sp=="DO"] - mean(data$D[data$tree_sp=="DO" & data$isolate=="CTRL"])
data$D_hat[data$tree_sp=="GR"] <-
  data$D[data$tree_sp=="GR"] - mean(data$D[data$tree_sp=="GR" & data$isolate=="CTRL"])

## as a check, these should be zero:
# data %>% filter(isolate == "CTRL") %>% group_by(tree_sp) %>% summarise(mean(D_hat))



## avg_avoidance
## plot of average corrected avoidance metric for all morphotype-tree combos

data %>%
  filter(isolate != "CTRL") %>% 
  group_by(morphotype, tree_sp) %>% 
  summarise(CI = mean_cl_normal(D_hat)) %>% 
  ggplot() +
  geom_line(aes(x=morphotype,y=CI$y,color=tree_sp, group=tree_sp)) +
  geom_point(aes(x=morphotype,y=CI$y,color=tree_sp)) +
  geom_errorbar(width=.25,
                aes(x=morphotype,ymin=CI$ymin,ymax=CI$ymax,color=tree_sp,group=tree_sp)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  scale_y_continuous(breaks = seq(-.05,.2,.05)) +
  xlab("Morphotype") +
  ylab(expression(paste("Avoidance metric, ", widehat(italic(D)))))


data %>%
  filter(isolate != "CTRL") %>% 
  group_by(morphotype, tree_sp) %>% 
  summarise(mean_D = mean(D_hat), 
            se = sd(D_hat) / sqrt(n())) %>% 
  mutate(ymin = mean_D - se, ymax = mean_D + se) %>% 
  ggplot() +
  geom_line(aes(x=morphotype,y=mean_D,color=tree_sp, group=tree_sp)) +
  geom_point(aes(x=morphotype,y=mean_D,color=tree_sp)) +
  geom_errorbar(width=.25,
                aes(x=morphotype,ymin=ymin,ymax=ymax,color=tree_sp,group=tree_sp)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  scale_y_continuous(breaks = seq(-.05,.2,.05)) +
  xlab("Morphotype") +
  ylab(expression(paste("Avoidance metric, ", widehat(italic(D)))))
ggsave("../figures/avg_avoidance.pdf", height = 4, width = 5) 

data %>%
  group_by(isolate, tree_sp) %>% 
  summarise(CI = mean_cl_normal(D)) %>% 
  ggplot() +
  geom_point(aes(x=isolate,y=CI$y,color=tree_sp), position=position_dodge(width = .1)) +
  geom_errorbar(width=.25,
                aes(x=isolate,ymin=CI$ymin,ymax=CI$ymax,color=tree_sp,group=tree_sp), position=position_dodge(width = .1)) +
  geom_hline(yintercept=0, linetype="dashed") +
  xlab("Isolate") +
  ylab(expression(paste("Uncorrected avoidance metric, ", italic(D))))




## avoidance_model
## plot of best avoidance model
## see stan_diagnostics.R





## t-tests

data %>%
  filter(morphotype=="SNPV", tree_sp=="DO") %>% 
  select(D_hat) %>% 
  t.test(alternative="l")


data %>%
  filter(morphotype=="SNPV", tree_sp=="GR") %>% 
  select(D_hat) %>% 
  t.test(alternative="g")

data %>%
  filter(morphotype=="MNPV", tree_sp=="DO") %>% 
  select(D_hat) %>% 
  t.test(alternative="g")

data %>%
  filter(morphotype=="MNPV", tree_sp=="GR") %>% 
  select(D_hat) %>% 
  t.test(alternative="g")



  


