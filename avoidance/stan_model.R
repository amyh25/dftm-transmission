library(tidyverse)
library(rstan)
library(bayesplot)

### settings
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(1000)
rstan_options(auto_write = TRUE)

input <- read.csv(file="manual_avoid_effect.csv",
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

### make data and fit stan model

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




looic <- sapply(c(fit_hier_tree_int,
                  fit_hier_tree,
                  fit_hier,
                  fit_tree,
                  fit_tree_only,
                  fit_treatment,
                  fit_intercept),
               function(fit) loo(fit)$estimate[3])
weight <- exp(-looic/2) / sum(exp(-looic/2))




