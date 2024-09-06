library(tidyverse)
library(rstan)
require(loo)
require(ggpubr)
theme_set(theme_pubr())

### settings
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(1000)

options(mc.cores = 1)
rstan_options(auto_write = TRUE)

input <- read.csv(file="../data/manual_avoid_effect.csv",
                 header=TRUE, sep=",", stringsAsFactors=FALSE)
input$morphotype <- factor(input$morphotype, levels = c("SNPV","MNPV"))
input$tree_sp <- factor(input$tree_sp, levels = c("GR","DO"))
input <- input %>% filter(control != "NA", treatment != "NA")
input$control[input$control < 0] <- 0
input$treatment[input$treatment < 0] <- 0
input <- input %>% cbind(D = input$control - input$treatment)
input <- input %>% cbind(list(D_hat = input$D))
input$D_hat[input$tree_sp=="DO"] <-
  input$D[input$tree_sp=="DO"] - mean(input$D[input$tree_sp=="DO" & input$isolate=="CTRL"])
input$D_hat[input$tree_sp=="GR"] <-
  input$D[input$tree_sp=="GR"] - mean(input$D[input$tree_sp=="GR" & input$isolate=="CTRL"])

input <- input %>% mutate(avoid_effect = D_hat, type = case_when(isolate == "CTRL" ~ "CTRL",
                                                                             TRUE ~ "TREAT"))
input$isolate <- factor(input$isolate, levels = c("CTRL", "COL", "LOV", "DRY", "TAM"))
input$tree_sp <- factor(input$tree_sp, levels=c("GR", "DO"))

data <- input


#### load fits if available
## might have to manually open the files first to decompress them
# fit_hier_tree_int <- readRDS("stan_fits/fit_hier_tree_int.rds")
# fit_hier_tree <- readRDS("stan_fits/fit_hier_tree.rds")
# fit_hier <- readRDS("stan_fits/fit_hier.rds")
# fit_tree <- readRDS("stan_fits/fit_tree.rds")
# fit_tree_only <- readRDS("stan_fits/fit_tree_only.rds")
# fit_treatment <- readRDS("stan_fits/fit_treatment.rds")
# fit_intercept <- readRDS("stan_fits/fit_intercept.rds")




###############
# FUNCTIONS
################

make_hier_tree_int_data <- function() {
  n_obs <- nrow(input)
  X <- model.matrix(~ tree_sp * isolate, data = input)
  data <- list(
    N = n_obs, 
    K = ncol(X), 
    X = X, 
    y = input$avoid_effect, 
    map = c(1, 2, 3, 3, 4, 4, 3, 3, 4, 4), 
    hier = 6
  )
  return (data)
}

make_hier_tree_data <- function() {
  n_obs <- nrow(input)
  X <- model.matrix(~ tree_sp + isolate, data = input)
  data <- list(
    N = n_obs, 
    K = ncol(X), 
    X = X, 
    y = input$avoid_effect, 
    map = c(1, 2, 3, 3, 4, 4), 
    hier = 2
  )
  return (data)
}

make_hier_data <- function() {
  n_obs <- nrow(input)
  X <- model.matrix(~ isolate, data = input)
  data <- list(
    N = n_obs, 
    K = ncol(X), 
    X = X, 
    y = input$avoid_effect, 
    map = c(1, 2, 2, 3, 3), 
    hier = 2
  )
  return (data)
}

make_tree_data <- function() {
  n_obs <- nrow(input)
  X <- model.matrix(~ type + tree_sp, data = input)
  data <- list(
    N = n_obs, 
    K = ncol(X),
    X = X, 
    y = input$avoid_effect, 
    map = c(1, 2, 3), 
    hier = 0
  )
  return (data)
}

make_tree_only_data <- function() {
  n_obs <- nrow(input)
  X <- model.matrix(~ tree_sp, data = input)
  data <- list(
    N = n_obs, 
    K = ncol(X), 
    X = X, 
    y = input$avoid_effect, 
    map = c(1, 2), 
    hier = 0
  )
  return (data)
}


make_treatment_data <- function() {
  n_obs <- nrow(input)
  X <- model.matrix(~ type, data = input)
  data <- list(
    N = n_obs, 
    K = ncol(X),
    X = X, 
    y = input$avoid_effect, 
    map = c(1, 2), 
    hier = 0
  )
  return (data)
}


make_intercept_data <- function() {
  n_obs <- nrow(input)
  X <- model.matrix(~ 1, data = input)
  data <- list(
    N = n_obs, 
    K = ncol(X), 
    X = X, 
    y = input$avoid_effect
  )
  return (data)
}

##################
# ALL PARAMETERS #
##################

### save fit
# fit@stanmodel@dso <- new("cxxdso")
saveFit <- function(fit_obj, fit_name) {
  saveRDS(fit_obj, file=fit_name)
  fit <- fit_obj
  rm(fit_obj)
  return(fit)
}

### make data and fit stan models

nchains <- 5
nwarmup <- 1000
niter <- 1e4
nthin <- 1


data_hier_tree_int <- make_hier_tree_int_data()
data_hier_tree <- make_hier_tree_data()
data_hier <- make_hier_data()
data_tree <- make_tree_data()
data_tree_only <- make_tree_only_data()
data_treatment <- make_treatment_data()
data_intercept <- make_intercept_data()



fit_hier_tree_int <- stan(file = "model.stan", data=data_hier_tree_int, 
                          chains=nchains, warmup=nwarmup, iter=niter, thin=nthin)
saveFit(fit_hier_tree_int, "stan_fits/fit_hier_tree_int.rds")

fit_hier_tree <- stan(file = "model.stan", data = data_hier_tree, 
                 chains = nchains, warmup = nwarmup, iter = niter,thin = nthin)
saveFit(fit_hier_tree, "stan_fits/fit_hier_tree.rds")

fit_hier <- stan(file = "model.stan", data = data_hier, 
                 chains = nchains, warmup = nwarmup, iter = niter,thin = nthin)
saveFit(fit_hier, "stan_fits/fit_hier.rds")

fit_tree <- stan(file = "model.stan", data = data_tree, 
                 chains = nchains, warmup = nwarmup, iter = niter,thin = nthin)
saveFit(fit_tree, "stan_fits/fit_tree.rds")

fit_tree_only <- stan(file = "model.stan", data=data_tree_only, 
                      chains=nchains, warmup=nwarmup, iter=niter, thin=nthin)
saveFit(fit_tree_only, "stan_fits/fit_tree_only.rds")

fit_treatment <- stan(file = "model.stan", data=data_treatment, 
                      chains=nchains, warmup=nwarmup, iter=niter, thin=nthin)
saveFit(fit_treatment, "stan_fits/fit_treatment.rds")

fit_intercept <- stan(file = "model_intercept.stan", data = data_intercept,
                      chains = nchains, warmup = nwarmup, iter = niter, thin = nthin)
saveFit(fit_intercept, "stan_fits/fit_intercept.rds")



#### diagnostics

loo1 <- loo(fit_hier_tree_int)
loo2 <- loo(fit_hier_tree)
loo3 <- loo(fit_hier)
loo4 <- loo(fit_tree)
loo5 <- loo(fit_tree_only)
loo6 <- loo(fit_treatment)
loo7 <- loo(fit_intercept)

loo_table <- loo_compare(loo1, loo2, loo3, loo4, loo5, loo6, loo7)
loo_table <- cbind(elpd=c(loo1$estimates[1,1],
                          loo2$estimates[1,1],
                          loo3$estimates[1,1],
                          loo4$estimates[1,1],
                          loo5$estimates[1,1],
                          loo6$estimates[1,1],
                          loo7$estimates[1,1])[order(order(rownames(loo_table)))],
                   loo_table[,1:2])
rownames(loo_table) <-
  c("treatment + morphotype + tree sp + interaction",
    "treatment + morphotype + tree sp",
    "treatment + morphotype",
    "treatment + tree sp",
    "tree sp",
    "treatment",
    "intercept")[order(order(rownames(loo_table)))]
loo_table


## avoidance_model
## plot of best avoidance model

fit <- fit_hier_tree_int
estimates <- summary(fit)$summary[,"mean"]
betas <- estimates[str_detect(names(estimates),pattern="^beta")]


data <- data %>% cbind(prediction = model.matrix(~ tree_sp * isolate, data = input) %*% betas)

data %>%
  filter(isolate != "CTRL") %>% 
  group_by(tree_sp, isolate) %>% 
  summarise(CI=mean_cl_normal(D_hat), prediction=mean(prediction)) %>%
  ungroup() %>% 
  ggplot() +
  geom_point(aes(x=tree_sp,y=CI$y,color=tree_sp), size=1.7) +
  geom_errorbar(width=.31,
                aes(x=tree_sp,ymin=CI$ymin,ymax=CI$ymax,color=tree_sp,group=tree_sp)) +
  geom_point(aes(x=tree_sp,y=prediction), shape=4, size=3.2) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~isolate, nrow=2, scales="free_x") +
  scale_x_discrete(labels = str_wrap(c("Grand fir", "Douglas fir"), width=7)) +
  scale_y_continuous(breaks = seq(-.1,.3,.1)) +
  xlab("Tree species") +
  ylab(expression(paste("Avoidance metric, ", widehat(italic(D))))) +
  theme(strip.background = element_blank(),
        panel.spacing= unit(1.5, "lines"),
        legend.position = "none",
        axis.text.x=element_text(size=10))





## avg_avoidance
## plot of average corrected avoidance metric for all morphotype-tree combos

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




## t-tests

t.test(data$treatment, data$control, alternative="l")
std_err(data$treatment)
std_err(data$control)

t.test(data[data$isolate=="CTRL","treatment"],data[data$isolate=="CTRL","control"])
std_err(data[data$isolate=="CTRL","treatment"])
std_err(data[data$isolate=="CTRL","control"])


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

