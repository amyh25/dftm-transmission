require(tidyverse)
require(rstan) # had to delete contents of .Renviron and .R/Makevars to install rstan on R 4.0
require(pracma)
require(loo)

### load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
SOK_data <- read.csv(file="SOK_data.csv",
                               header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
            filter(capsid != "none")

data_MNPV_first <- SOK_data
data_SNPV_first <- SOK_data
data_SNPV_first$capsid <- factor(data_MNPV_first$capsid,levels=c("SNPV","MNPV"))


### set up stan
set.seed(683038)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

factor_to_int <- function(factor, col_name) {
  factors <- t(unique(SOK_data[col_name]))
  dict <- as.list(1:length(factors))
  names(dict) <- factors
  return(dict[[factor]])
}

strain_factor_to_int = function(i) factor_to_int(i,"strain")
tree_factor_to_int = function(i) factor_to_int(i,"tree_sp")

N <- nrow(SOK_data) # number of treatments = 48
H <- length(unique(SOK_data$capsid)) # number of capsids = 2
I <- length(unique(SOK_data$strain)) # number of strains = 8
J <- length(unique(SOK_data$tree_sp)) # number of tree species = 2

cid <- c(rep(1,4),rep(2,4)) # SNPV = 1, MNPV = 2
sid <- sapply(SOK_data$strain, strain_factor_to_int) # COL=1,CUB=2,LOV=3,LST=4,DRY=5,KLP=6,TAM=7,TMB=8
tid <- sapply(SOK_data$tree_sp, tree_factor_to_int) # GR = 1, DO = 2

x <- SOK_data$ob_count
y <- SOK_data$total_virus
total <- SOK_data$total_n


## morphotype_and_tree

glm_SNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_SNPV_first, weights=total_n)
glm_MNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_MNPV_first, weights=total_n)

prior_mu_alpha <- matrix(nrow = I, ncol = J)
prior_mu_beta <- matrix(nrow = I, ncol = J)
p <- as.vector(predict(glm_SNPV_first))
for (n in seq(1,N,3)) {
  prior_mu_beta[sid[n],tid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  if (prior_mu_beta[sid[n],tid[n]] < 0) {
    prior_mu_beta[sid[n],tid[n]] <- 0
    prior_mu_alpha[sid[n],tid[n]] <- mean(p[n:(n+2)])
  } else {
    prior_mu_alpha[sid[n],tid[n]] <- p[n] - prior_mu_beta[sid[n],tid[n]] * x[n]
  }
}

prior_sigma_alpha <- c(Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),
                       Norm(coef(summary(glm_MNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]))
prior_sigma_beta <- c(Norm(coef(summary(glm_SNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]),
                      Norm(coef(summary(glm_MNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]))

morphotype_and_tree <- stan(file="morphotype_and_tree.stan",
                        data=list(N=N,H=H,I=I,J=J,
                                  cid=cid,sid=sid,tid=tid,
                                  x=x,y=y,total=total,
                                  prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                                  prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                        chains=4,
                        iter=5000,
                        control = list(adapt_delta=0.99, max_treedepth=25))




## tree_only

glm_SNPV_first <- glm(dose_response ~ ob_count * tree_sp,
                      family = "binomial", data = data_SNPV_first, weights=total_n)

prior_mu_alpha <- rep(0,J)
prior_mu_beta <- rep(0,J)
p <- as.vector(predict(glm_SNPV_first))
for (n in c(1,4)) {
  prior_mu_beta[tid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  if (prior_mu_beta[tid[n]] < 0) {
    prior_mu_beta[tid[n]] <- 0
    prior_mu_alpha[tid[n]] <- mean(p[n:(n+2)])
  } else {
    prior_mu_alpha[tid[n]] <- p[n] - prior_mu_beta[tid[n]] * x[n]
  }
}

prior_sigma_alpha <- Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"])
prior_sigma_beta <- Norm(coef(summary(glm_SNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"])

tree_only <- stan(file="tree_only.stan",
                  data=list(N=N,I=I,J=J,
                            sid=sid,tid=tid,
                            x=x,y=y,total=total,
                            prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                            prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                  chains=4,
                  iter=5000,
                  control = list(adapt_delta=0.99, max_treedepth=25))



## morphotype_only

glm_SNPV_first <- glm(dose_response ~ ob_count * capsid,
                      family = "binomial", data = data_SNPV_first, weights=total_n)
glm_MNPV_first <- glm(dose_response ~ ob_count * capsid,
                      family = "binomial", data = data_MNPV_first, weights=total_n)

prior_mu_alpha <- rep(0,I)
prior_mu_beta <- rep(0,I)
p <- as.vector(predict(glm_SNPV_first))
for (n in seq(1,N,6)) {
  prior_mu_beta[sid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  if (prior_mu_beta[sid[n]] < 0) {
    prior_mu_beta[sid[n]] <- 0
    prior_mu_alpha[sid[n]] <- mean(p[n:(n+5)])
  } else {
    prior_mu_alpha[sid[n]] <- p[n] - prior_mu_beta[sid[n]] * x[n]
  }
}

prior_sigma_alpha <- c(coef(summary(glm_SNPV_first))["(Intercept)", "Std. Error"],
                       coef(summary(glm_MNPV_first))["(Intercept)", "Std. Error"])
prior_sigma_beta <- c(coef(summary(glm_SNPV_first))["ob_count", "Std. Error"],
                      coef(summary(glm_MNPV_first))["ob_count", "Std. Error"])

morphotype_only <- stan(file="morphotype_only.stan",
                    data=list(N=N,H=H,I=I,
                              cid=cid,sid=sid,
                              x=x,y=y,total=total,
                              prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                              prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                    chains=4,
                    iter=5000,
                    control = list(adapt_delta=0.99, max_treedepth=25))


## neither_morphotype_nor_tree

glm_SNPV_first <- glm(dose_response ~ ob_count,
                      family = "binomial", data = data_SNPV_first, weights=total_n)

prior_mu_alpha <- coef(summary(glm_SNPV_first))["(Intercept)", "Estimate"]
prior_mu_beta <- coef(summary(glm_SNPV_first))["ob_count", "Estimate"]

prior_sigma_alpha <- coef(summary(glm_SNPV_first))["(Intercept)", "Std. Error"]
prior_sigma_beta <- coef(summary(glm_SNPV_first))["ob_count", "Std. Error"]
if (prior_sigma_beta < 0) {
  prior_sigma_beta <- 0
}

neither_morphotype_nor_tree <- stan(file="neither_morphotype_nor_tree.stan",
                                    data=list(N=N,I=I,
                                              sid=sid,
                                              x=x,y=y,total=total,
                                              prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                                              prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                                    chains=4,
                                    iter=5000,
                                    control = list(adapt_delta=0.99, max_treedepth=25))


## complete_hierarchy

glm_SNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_SNPV_first, weights=total_n)
glm_MNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_MNPV_first, weights=total_n)

prior_mu_alpha <- matrix(nrow = H, ncol = J)
prior_mu_beta <- matrix(nrow = H, ncol = J)
p <- as.vector(predict(glm_SNPV_first))
for (n in c(1,4,25,28)) {
  prior_mu_beta[cid[sid[n]],tid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  prior_mu_alpha[cid[sid[n]],tid[n]] <- p[n] - prior_mu_beta[cid[sid[n]],tid[n]] * x[n]
}

prior_sigma_alpha <- c(Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),
                       Norm(coef(summary(glm_MNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]))
prior_sigma_beta <- c(Norm(coef(summary(glm_SNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]),
                      Norm(coef(summary(glm_MNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]))

complete_hierarchy <- stan(file="complete_hierarchy.stan",
                           data=list(N=N,H=H,I=I,J=J,
                                     cid=cid,sid=sid,tid=tid,
                                     x=x,y=y,total=total,
                                     prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                                     prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                           chains=4,
                           iter=5000,
                           control = list(adapt_delta=0.99, max_treedepth=25))


## no_hierarchy

glm_SNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_SNPV_first, weights=total_n)
glm_MNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_MNPV_first, weights=total_n)

prior_mu_alpha <- matrix(nrow = I, ncol = J)
prior_mu_beta <- matrix(nrow = I, ncol = J)
p <- as.vector(predict(glm_SNPV_first))
for (n in seq(1,N,3)) {
  prior_mu_beta[sid[n],tid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  prior_mu_alpha[sid[n],tid[n]] <- p[n] - prior_mu_beta[sid[n],tid[n]] * x[n]
}

prior_sigma_alpha <- c(rep(Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),4),
                       rep(Norm(coef(summary(glm_MNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),4))
prior_sigma_beta <- c(rep(Norm(coef(summary(glm_SNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]),4),
                      rep(Norm(coef(summary(glm_MNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]),4))

no_hierarchy <- stan(file="no_hierarchy.stan",
                        data=list(N=N,I=I,J=J,
                                  sid=sid,tid=tid,
                                  x=x,y=y,total=total,
                                  prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                                  prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                        chains=4,
                        iter=5000,
                        control = list(adapt_delta=0.99, max_treedepth=25))



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
fit <- saveFit(complete_hierarchy, "stan_fits/complete_hierarchy.rds")
fit <- saveFit(no_hierarchy, "stan_fits/no_hierarchy.rds")



