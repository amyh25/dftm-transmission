### run after SOK_models.R, or at least ensure the stan fits exist for the SOK models

require(deSolve)
require(tidyverse)
require(rstan)
require(ggpubr)
theme_set(theme_pubr())



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

SOK_complete <-
  tidy %>%
  count(capsid,tree_sp,numeric_day) %>%
  complete(capsid,tree_sp,numeric_day, fill=list(n=0))



### load fits from SOK_models.R
morphotype_and_tree <- readRDS("stan_fits/morphotype_and_tree.rds")
neither_morphotype_nor_tree <- readRDS("stan_fits/neither_morphotype_nor_tree.rds")

pars_diff_dists <- summary(morphotype_and_tree)$summary[1:8,"mean"]
pars_same_dist <- summary(neither_morphotype_nor_tree)$summary[1:2,"mean"]

alpha_SNPV_GR <- pars_diff_dists["alpha[1,1]"]
beta_SNPV_GR <- pars_diff_dists["beta[1,1]"]
alpha_SNPV_DO <- pars_diff_dists["alpha[1,2]"]
beta_SNPV_DO <- pars_diff_dists["beta[1,2]"]
alpha_MNPV_GR <- pars_diff_dists["alpha[2,1]"]
beta_MNPV_GR <- pars_diff_dists["beta[2,1]"]
alpha_MNPV_DO <- pars_diff_dists["alpha[2,2]"]
beta_MNPV_DO <- pars_diff_dists["beta[2,2]"]

alpha_same_dist <- pars_same_dist["alpha"]
beta_same_dist <- pars_same_dist["beta"]

SOK_var_model <- function(t,y,p){
  with(as.list(p), {
    S = y[1]
    Es = y[2:(m+1)]
    P = y[m+2]
    
    
    dS.dt = -nu * S * P * (S/S0)^(C^2)
    dEs.dt = c(nu * S * P * (S/S0)^(C^2) - m * delta * Es[1], -m * delta * diff(Es))
    dP.dt = m * delta * Es[m] - mu * P
    
    return(list(c(dS.dt, dEs.dt, dP.dt)))
  })
}

P0 <- .01
C <- .978
nu <- .024
mu <- 0
ts <- seq(0,50,.1) # a 50-day epizootic

data <- c()

for (diff_variances in c(FALSE, TRUE)) {
  for (diff_means in c(FALSE, TRUE)) {
    for (S0 in seq(.5,100,.5)) {
      print(S0)
      
      # MNPV_GR
      mean <- if (diff_means) alpha_MNPV_GR * beta_MNPV_GR else alpha_same_dist * beta_same_dist
      variance <- if (diff_variances) alpha_MNPV_GR * beta_MNPV_GR^2 else alpha_same_dist * beta_same_dist^2
      m <- mean^2 / variance
      m <- if (m %% 1 == .5) ceiling(m) else round(m) + (round(m)==0)
      p <- list(S0 = S0, C = C, nu = nu, mu = mu, m = m, delta = 1/mean)
      y0 <- c(S0, rep(0,m), P0)
      out_MNPV_GR = dede(y=y0, times=ts, func=SOK_var_model, parms=p)
      
      # MNPV_DO
      mean <- if (diff_means) alpha_MNPV_DO * beta_MNPV_DO else alpha_same_dist * beta_same_dist
      variance <- if (diff_variances) alpha_MNPV_DO * beta_MNPV_DO^2 else alpha_same_dist * beta_same_dist^2
      m <- mean^2 / variance
      m <- if (m %% 1 == .5) ceiling(m) else round(m) + (round(m)==0)
      p <- list(S0 = S0, C = C, nu = nu, mu = mu, m = m, delta = 1/mean)
      y0 <- c(S0, rep(0,m), P0)
      out_MNPV_DO = dede(y=y0, times=ts, func=SOK_var_model, parms=p)
      
      # SNPV_GR
      mean <- if (diff_means) alpha_SNPV_GR * beta_SNPV_GR else alpha_same_dist * beta_same_dist
      variance <- if (diff_variances) alpha_SNPV_GR * beta_SNPV_GR^2 else alpha_same_dist * beta_same_dist^2
      m <- mean^2 / variance
      m <- if (m %% 1 == .5) ceiling(m) else round(m) + (round(m)==0)
      p <- list(S0 = S0, C = C, nu = nu, mu = mu, m = m, delta = 1/mean)
      y0 <- c(S0, rep(0,m), P0)
      out_SNPV_GR = dede(y=y0, times=ts, func=SOK_var_model, parms=p)
      
      # SNPV_DO
      mean <- if (diff_means) alpha_SNPV_DO * beta_SNPV_DO else alpha_same_dist * beta_same_dist
      variance <- if (diff_variances) alpha_SNPV_DO * beta_SNPV_DO^2 else alpha_same_dist * beta_same_dist^2
      m <- mean^2 / variance
      m <- if (m %% 1 == .5) ceiling(m) else round(m) + (round(m)==0)
      p <- list(S0 = S0, C = C, nu = nu, mu = mu, m = m, delta = 1/mean)
      y0 <- c(S0, rep(0,m), P0)
      out_SNPV_DO = dede(y=y0, times=ts, func=SOK_var_model, parms=p)
      
      
      data <-
        rbind(data,
              data.frame(diff_variances = diff_variances,
                         diff_means = diff_means,
                         S0 = S0,
                         morphotype = factor(rep(c("MNPV","SNPV"), each=2), levels=c("SNPV","MNPV")),
                         tree_sp = factor(rep(c("GR","DO"), 2), levels=c("GR","DO")),
                         y = 1-c(out_MNPV_GR[length(ts),2],
                                 out_MNPV_DO[length(ts),2],
                                 out_SNPV_GR[length(ts),2],
                                 out_SNPV_DO[length(ts),2])/S0))
    }
  }
}




## SOK_var_same_means_diff_vars

data %>%
  filter(!diff_means, diff_variances) %>%
  ggplot() +
  geom_line(aes(x=S0, y=y, group=interaction(morphotype,tree_sp), color=tree_sp, lty=morphotype)) +
  scale_x_continuous(name=expression(paste("Initial density of susceptible larvae (per ",m^2,")")),limits=c(0,100),expand=expansion(c(0,0))) +
  scale_y_continuous(name="Total fraction of larvae infected",
                     expand=expansion(c(0,0)), limits=c(0,1)) +
  scale_color_discrete(name="Tree", labels=c("Grand fir","Douglas-fir")) +
  scale_linetype_manual(name="Morphotype", values=c(2,1)) +
  ggtitle("Same means, different variances")



## SOK_var_diff_means_diff_vars

data %>%
  filter(diff_means, diff_variances) %>%
  ggplot() +
  geom_line(aes(x=S0, y=y, group=interaction(morphotype,tree_sp), color=tree_sp, lty=morphotype)) +
  scale_x_continuous(name=expression(paste("Initial density of susceptible larvae (per ",m^2,")")),limits=c(0,100),expand=expansion(c(0,0))) +
  scale_y_continuous(name="Total fraction of larvae infected",
                     expand=expansion(c(0,0)), limits=c(0,1)) +
  scale_color_discrete(name="Tree", labels=c("Grand fir","Douglas-fir")) +
  scale_linetype_manual(name="Morphotype", values=c(2,1)) +
  ggtitle("Different means, different variances")



## SOK_var_diff_means_same_vars in supplement

data %>%
  filter(diff_means, !diff_variances) %>%
  ggplot() +
  geom_line(aes(x=S0, y=y, group=interaction(morphotype,tree_sp), color=tree_sp, lty=morphotype)) +
  scale_x_continuous(name=expression(paste("Initial density of susceptible larvae (per ",m^2,")")),limits=c(0,100),expand=expansion(c(0,0))) +
  scale_y_continuous(name="Total fraction of larvae infected",
                     expand=expansion(c(0,0)), limits=c(0,1)) +
  scale_color_discrete(name="Tree", labels=c("Grand fir","Douglas-fir")) +
  scale_linetype_manual(name="Morphotype", values=c(2,1)) +
  ggtitle("Different means, same variances")




