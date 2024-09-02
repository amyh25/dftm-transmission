require(tidyverse)
require(rstan) # had to delete contents of .Renviron and .R/Makevars to install rstan on R 4.0
require(pracma)
require(loo)
require(ggpubr)
theme_set(theme_pubr())

#### load fits if available
# morphotype_and_tree <- readRDS("stan_fits/morphotype_and_tree.rds")
# tree_only <- readRDS("stan_fits/tree_only.rds")
# morphotype_only <- readRDS("stan_fits/morphotype_only.rds")
# neither_morphotype_nor_tree <- readRDS("stan_fits/neither_morphotype_nor_tree.rds")


### load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
SOK_data <- read.csv(file="../data/SOK_data.csv",
                     header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
  filter(capsid != "none")

SOK_data$strain <- factor(SOK_data$strain,levels=c("COL","CUB","LOV","LST","DRY","KLP","TAM","TMB"))
SOK_data$capsid <- factor(SOK_data$capsid,levels=c("SNPV","MNPV"))
SOK_data$tree_sp <- factor(SOK_data$tree_sp,levels=c("GR","DO"))


### tidy the data
days <- colnames(SOK_data)[9:31]

day_num <- function(str) {
  return(substring(str,5) %>% strtoi())
}

tidy <- SOK_data %>%
  pivot_longer(days, names_to="day")

tidy <- tidy %>% cbind(numeric_day = unlist(lapply(tidy$day, day_num)))

tidy <- tidy[rep(1:nrow(tidy), times=tidy$value),]


### set up stan
set.seed(134523)

options(mc.cores = 1)
rstan_options(auto_write = FALSE)

factor_to_int <- function(factor, col_name) {
  factors <- t(unique(SOK_data[col_name]))
  dict <- as.list(1:length(factors))
  names(dict) <- factors
  return(dict[[factor]])
}

capsid_factor_to_int = function(i) factor_to_int(i,"capsid")
tree_factor_to_int = function(i) factor_to_int(i,"tree_sp")

N <- nrow(tidy) # number of virus-killed larvae = 421
H <- length(unique(tidy$capsid)) # number of capsids = 2
J <- length(unique(tidy$tree_sp)) # number of tree species = 2

cid <- sapply(tidy$capsid, capsid_factor_to_int) # SNPV = 1, MNPV = 2
tid <- sapply(tidy$tree_sp, tree_factor_to_int) # GR = 1, DO = 2



### stan models

morphotype_and_tree <- stan(file="morphotype_and_tree.stan",
                            data=list(N=N,H=H,J=J,
                                      cid=cid,tid=tid,
                                      y=tidy$numeric_day),
                            chains=4,
                            iter=4000,
                            init=1,
                            control = list(adapt_delta=0.99, max_treedepth=20))


tree_only <- stan(file="tree_only.stan",
                  data=list(N=N,J=J,
                            tid=tid,
                            y=tidy$numeric_day),
                  chains=4,
                  iter=4000,
                  init=1,
                  control = list(adapt_delta=0.99, max_treedepth=20))


morphotype_only <- stan(file="morphotype_only.stan",
                            data=list(N=N,H=H,
                                      cid=cid,
                                      y=tidy$numeric_day),
                        chains=4,
                        iter=4000,
                        init=1,
                        control = list(adapt_delta=0.99, max_treedepth=20))


neither_morphotype_nor_tree <- stan(file="neither_morphotype_nor_tree.stan",
                            data=list(N=N,
                                      y=tidy$numeric_day),
                            chains=4,
                            iter=4000,
                            init=1,
                            control = list(adapt_delta=0.99, max_treedepth=20))



### save fits
saveFit <- function(fit_obj, fit_name) {
  fit_obj@stanmodel@dso <- new("cxxdso")
  saveRDS(fit_obj, file=fit_name)
  fit <- fit_obj
  rm(fit_obj)
  return(fit)
}

fit <- saveFit(morphotype_and_tree, "stan_fits/morphotype_and_tree.rds")
fit <- saveFit(tree_only, "stan_fits/tree_only.rds")
fit <- saveFit(morphotype_only, "stan_fits/morphotype_only.rds")
fit <- saveFit(neither_morphotype_nor_tree, "stan_fits/neither_morphotype_nor_tree.rds")


#### diagnostics

loo1 <- loo(morphotype_and_tree)
loo2 <- loo(tree_only)
loo3 <- loo(morphotype_only)
loo4 <- loo(neither_morphotype_nor_tree)

loo_table <- loo_compare(loo1, loo2, loo3, loo4)
rownames(loo_table) <- c("M and T","T only","M only","Neither M nor T")[order(order(rownames(loo_table)))]
loo_table





## SOK_model
## plot of best SOK model

fit_hier <- morphotype_and_tree # best model

pars <- summary(fit_hier)$summary[1:8,"mean"]

model_output <- data.frame()
xs <- seq(4,27,.01)

for (x in xs) {
  model_output <-
    rbind(model_output,
          list(capsid="SNPV", tree_sp="GR", x=x, y=dgamma(x+.5,shape=pars["alpha[1,1]"],scale=pars["beta[1,1]"])),
          list(capsid="SNPV", tree_sp="DO", x=x, y=dgamma(x+.5,shape=pars["alpha[1,2]"],scale=pars["beta[1,2]"])),
          list(capsid="MNPV", tree_sp="GR", x=x, y=dgamma(x+.5,shape=pars["alpha[2,1]"],scale=pars["beta[2,1]"])),
          list(capsid="MNPV", tree_sp="DO", x=x, y=dgamma(x+.5,shape=pars["alpha[2,2]"],scale=pars["beta[2,2]"])))
}

model_output$capsid <- factor(model_output$capsid,levels=c("SNPV","MNPV"))
model_output$tree_sp <- factor(model_output$tree_sp,levels=c("GR","DO"))
tidy$capsid <- factor(tidy$capsid,levels=c("SNPV","MNPV"))
tidy$tree_sp <- factor(tidy$tree_sp,levels=c("GR","DO"))


tree.labs <- c("Grand fir","Douglas fir")
names(tree.labs) <- c("GR","DO")

ggplot() +
  geom_histogram(data=tidy,aes(x=numeric_day,y=stat(density),fill=tree_sp),
                 binwidth=1,boundary=0,color="black",size=.1) +
  geom_line(data=model_output,aes(x=x,y=y)) +
  facet_grid(capsid~tree_sp,labeller=labeller(tree_sp=tree.labs)) +
  coord_cartesian(xlim=c(4,27),ylim=c(0,.2)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_discrete(name = "Tree", labels = c("Grand fir", "Douglas fir")) +
  xlab("Speed of kill (days)") +
  ylab("Proportion killed") +
  theme(strip.background = element_blank(),legend.position = "none")



## avg_SOK
## plot of average speed of kill for morphotype-tree combos

avg_SOK %>% 
  ggplot() + 
  aes(x = capsid, y = y, color = tree_sp) + 
  geom_point() + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.25) + 
  geom_line(aes(group = tree_sp)) + 
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  xlab("Morphotype") +
  ylab("Average speed of kill (days)") +
  theme(legend.position="bottom")




## t-tests

t.test(tidy[tidy$capsid=="MNPV" & tidy$tree_sp=="DO","numeric_day"],
       tidy[tidy$capsid=="MNPV" & tidy$tree_sp=="GR","numeric_day"], alternative = "l")
std_err(tidy[tidy$capsid=="MNPV" & tidy$tree_sp=="DO","numeric_day"])
std_err(tidy[tidy$capsid=="MNPV" & tidy$tree_sp=="GR","numeric_day"])


t.test(tidy[tidy$capsid=="SNPV" & tidy$tree_sp=="DO","numeric_day"],
       tidy[tidy$capsid=="SNPV" & tidy$tree_sp=="GR","numeric_day"], alternative = "g")
std_err(tidy[tidy$capsid=="SNPV" & tidy$tree_sp=="DO","numeric_day"])
std_err(tidy[tidy$capsid=="SNPV" & tidy$tree_sp=="GR","numeric_day"])

