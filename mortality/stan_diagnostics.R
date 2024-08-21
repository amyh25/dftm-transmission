require(rstan)
require(tidyverse)
require(shinystan)
require(bayesplot)
require(coda)
require(loo)
require(pracma)
require(binom)
require(ggpubr)
require(boot)
theme_set(theme_pubr())

### load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
SOK_data <- read.csv(file="SOK_data.csv",
                     header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
            filter(capsid != "none")
SOK_data$strain <- factor(SOK_data$strain,levels=c("COL","CUB","LOV","LST","DRY","KLP","TAM","TMB"))
SOK_data$capsid <- factor(SOK_data$capsid,levels=c("SNPV","MNPV"))
SOK_data$tree_sp <- factor(SOK_data$tree_sp,levels=c("GR","DO"))
y <- SOK_data$dose_response

### load fits
# might have to manually open the files first to decompress them
morphotype_and_tree <- readRDS("stan_fits/morphotype_and_tree.rds")
tree_only <- readRDS("stan_fits/tree_only.rds")
morphotype_only <- readRDS("stan_fits/morphotype_only.rds")
neither_morphotype_nor_tree <- readRDS("stan_fits/neither_morphotype_nor_tree.rds")
complete_hierarchy <- readRDS("stan_fits/complete_hierarchy.rds")
no_hierarchy <- readRDS("stan_fits/no_hierarchy.rds")

# LOOIC table with supplementary models
bhms <- c(morphotype_and_tree,
          tree_only,morphotype_only,neither_morphotype_nor_tree,no_hierarchy,complete_hierarchy)
name <- c("M and T","T only","M only","Neither M nor T","No hierarchy","Complete hierarchy")
looic <- sapply(bhms,function(model) loo(model)$estimates[3])
delta_looic <- looic - min(looic)
weight <- exp(-looic/2) / sum(exp(-looic/2))
p_looic <- sapply(bhms,function(model) loo(model)$estimates[2])
looic_table <- data.frame(name,looic,delta_looic,weight,p_looic)

# LOOIC table without supplementary models
bhms <- c(morphotype_and_tree,
          tree_only,morphotype_only,neither_morphotype_nor_tree)
name <- c("M and T","T only","M only","Neither M nor T")
looic <- sapply(bhms,function(model) loo(model)$estimates[3])
delta_looic <- looic - min(looic)
weight <- exp(-looic/2) / sum(exp(-looic/2))
p_looic <- sapply(bhms,function(model) loo(model)$estimates[2])
looic_table <- data.frame(name,looic,delta_looic,weight,p_looic)




fit_hier <- morphotype_and_tree # best model



## mortality_model
## plot of best mortality model

factor_to_int <- function(factor, col_name) {
  factors <- t(unique(SOK_data[col_name]))
  dict <- as.list(1:length(factors))
  names(dict) <- factors
  return(dict[[factor]])
}

strain_factor_to_int = function(i) factor_to_int(i,"strain")
tree_factor_to_int = function(i) factor_to_int(i,"tree_sp")

estimates <- summary(fit_hier)$summary[,"mean"]
alphas <- estimates[str_detect(names(estimates),pattern="^alpha")]
betas <- estimates[str_detect(names(estimates),pattern="^beta")]

doses <- seq(0,6500,100)
predictions <- data.frame(dose=numeric(16*length(doses)),
                             strain=character(16*length(doses)),
                             tree_sp=character(16*length(doses)),
                             prediction=numeric(16*length(doses)),
                             stringsAsFactors = FALSE)
i <- 1
for (s in unique(SOK_data[,"strain"])) {
  s_int <- strain_factor_to_int(s)
  for (t in unique(SOK_data[,"tree_sp"])) {
    t_int <- tree_factor_to_int(t)
    for (d in doses) {
      predictions[i,] <- list(dose=d,
                                 strain=s,
                                 tree_sp=t,
                                 prediction = alphas[paste("alpha[",s_int,",",t_int,"]",sep="")] +
                                              betas[paste("beta[",s_int,",",t_int,"]",sep="")] * d)
      i <- i+1
    }
  }
}

predictions$strain <- factor(predictions$strain, levels = unique(SOK_data[,"strain"]))
predictions$tree_sp <- factor(predictions$tree_sp, levels = unique(SOK_data[,"tree_sp"]))


SOK_data_grouped_isolate <-
  SOK_data %>%
  group_by(strain,density,tree_sp) %>%
  summarise(total_virus=sum(total_virus) + .5 * (sum(total_virus)==0),
            total_n=sum(total_n),
            dose_response=logit(sum(total_virus)/sum(total_n)),
            dose_response_lower=logit(binom.confint(total_virus,total_n,method="wilson",conf.level=.95)$lower),
            dose_response_upper=logit(binom.confint(total_virus,total_n,method="wilson",conf.level=.95)$upper),
            dose_var=mean(ob_count))

ggplot(SOK_data_grouped_isolate) +
  geom_point(aes(x=dose_var/1000,y=dose_response,color=tree_sp)) +
  geom_errorbar(aes(x=dose_var/1000,ymin=dose_response_lower,ymax=dose_response_upper,
                    color=tree_sp),width=.30) +
  geom_line(data=predictions,aes(x=dose/1000,y=prediction,color=tree_sp)) +
  facet_wrap(~strain,nrow=2,scales="free_x") +
  geom_hline(yintercept=0, linetype="dashed",linewidth=.3) +
  scale_x_continuous(limits=c(0,6.5),expand = c(0,0)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6),limits=c(-6.5,6),expand = c(0,0)) +
  xlab("Dose (thousands of occlusion bodies)") + ylab("logit (Proportion virus-killed)") +
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  theme(strip.background = element_blank(),
        panel.spacing= unit(1.5, "lines"),
        plot.margin = margin(0,20,0,10))




# pulling out the fits
# NUTS (the No-U-Turn Sampler) optimizes HMC adaptively
posterior_fit_hier <- as.array(fit_hier)
nuts_params_fit_hier <- nuts_params(fit_hier)
log_posterior_fit_hier <- log_posterior(fit_hier)

#summary(fit_hier, pars = c("sigma_alpha","sigma_beta"))$summary
#head(log_posterior_fit_hier)

# trace plots for the runs
# to make sure nothing is being wonky
mcmc_trace(posterior_fit_hier, regex="alpha",np=nuts_params_fit_hier)
mcmc_trace(posterior_fit_hier, regex="beta",np=nuts_params_fit_hier)



mcmc_dens(posterior_fit_hier, regex="^alpha.{2,}1",facet_args=list(ncol=8)) +
  scale_x_continuous(limits=c(-6,2),expand = c(0,0)) +
  theme(text=element_text(size=11,family="Palatino"),
        panel.spacing= unit(1.2, "lines")) +
  facet_text(FALSE)
mcmc_dens(posterior_fit_hier, regex="^alpha.{2,}2",facet_args=list(ncol=8)) +
  scale_x_continuous(limits=c(-6,2),expand = c(0,0)) +
  theme(text=element_text(size=11,family="Palatino"),
        panel.spacing= unit(1.2, "lines")) +
  facet_text(FALSE)

mcmc_dens(posterior_fit_hier, regex="^beta.{2,}1",facet_args=list(ncol=8)) +
  scale_x_continuous(limits=c(-.0004,.0006),expand = c(0,0),
                     breaks=c(-.0004,-.0002,0,.0002,.0004,.0006),
                     labels=c("-.0004","","0","",".0004","")) +
  theme(text=element_text(size=11,family="Palatino"),
        strip.background = element_blank(),
        panel.spacing= unit(1.2, "lines")) +
  facet_text(FALSE)
mcmc_dens(posterior_fit_hier, regex="^beta.{2,}2",facet_args=list(ncol=8)) +
  scale_x_continuous(limits=c(-.0004,.0006),expand = c(0,0),
                     breaks=c(-.0004,-.0002,0,.0002,.0004,.0006),
                     labels=c("-.0004","","0","",".0004","")) +
  theme(text=element_text(size=11,family="Palatino"),
        strip.background = element_blank(),
        panel.spacing= unit(1.2, "lines")) +
  facet_text(FALSE)


mcmc_dens(posterior_fit_hier, regex="sigma_alpha") +
  scale_x_continuous(limits=c(0,3),expand = c(0,0)) +
  theme(text=element_text(size=11,family="Palatino"),
        panel.spacing= unit(1.2, "lines"),
        plot.margin=margin(0,20,0,10)) +
  facet_text(FALSE)
mcmc_dens(posterior_fit_hier, regex="sigma_beta") +
  scale_x_continuous(limits=c(0,.0006),breaks=c(0,.0002,.0004,.0006),
                     labels=c("0",".0002",".0004",".0006"),expand = c(0,0)) +
  theme(text=element_text(size=11,family="Palatino"),
        panel.spacing= unit(1.2, "lines"),
        plot.margin=margin(0,20,0,10)) +
  facet_text(FALSE)


# looks at the divergence in the plots
color_scheme_set("red")
mcmc_nuts_divergence(nuts_params_fit_hier, log_posterior_fit_hier)

# calculates the rhats for the fit
# if they are above 1.05 that is bad
# still working on getting the Gelman-Rubin statistic
# that is just a split r-hat (I believe)
# last line plots the rhats
rhats <- rhat(fit_hier)
# print(rhats_simple)
mcmc_rhat(rhats)

# plots the fit and deviation of the parameters
plot(fit_hier, pars = "alpha")
plot(fit_hier, pars = "beta")
plot(fit_hier, pars = c("theta"), fill_color = "red", show_density=TRUE) +
  geom_point(data = data.frame(x_coords=logit(y),y_coords=48:1),
             mapping = aes(x = x_coords, y = y_coords),
             colour = "blue",
             shape = 4,
             size = 3)

f <- function(i) return(paste("inv_logit_theta[",i,"]",sep=""))
plot(fit_hier, pars = sapply(48:1,f), fill_color = "red", show_density=TRUE) +
  geom_point(data = data.frame(x_coords=y,y_coords=1:48),
             mapping = aes(x = x_coords, y = y_coords),
             colour = "blue",
             shape = 4,
             size = 3) +
  geom_hline(yintercept=25) +
  geom_hline(yintercept=c(7,13,19,31,37,43),linetype="dashed")

plot(fit_hier, pars = "sigma_alpha", show_density=TRUE)
plot(fit_hier, pars = "sigma_beta", show_density=TRUE)


# comparing posterior distributions for the sigmas
sigma_alpha_SNPV <- c(sapply(1:4,function(i) fit_hier@sim$samples[[i]]$`sigma_alpha[1]`))
sigma_alpha_MNPV <- c(sapply(1:4,function(i) fit_hier@sim$samples[[i]]$`sigma_alpha[2]`))
sigma_beta_SNPV <- c(sapply(1:4,function(i) fit_hier@sim$samples[[i]]$`sigma_beta[1]`))
sigma_beta_MNPV <- c(sapply(1:4,function(i) fit_hier@sim$samples[[i]]$`sigma_beta[2]`))

ks.test(sigma_alpha_SNPV,sigma_alpha_MNPV,alternative="greater")
ks.test(sigma_beta_SNPV,sigma_beta_MNPV,alternative="greater")

t.test(sigma_alpha_SNPV,sigma_alpha_MNPV)
t.test(sigma_beta_SNPV,sigma_beta_MNPV)


## measuring the predicted difference in percent mortality over the dose range measured

min_dose <- min(SOK_data$ob_count)
min_dose_floor <- floor(min_dose/100)*100
min_dose_ceiling <- ceiling(min_dose/100)*100
min_predictions_floor <- filter(predictions, dose==min_dose_floor) %>% select(prediction)
min_predictions_ceiling <- filter(predictions, dose==min_dose_ceiling) %>% select(prediction)
min_predictions_actual <- min_predictions_floor + (min_predictions_ceiling-min_predictions_floor)*(min_dose-min_dose_floor)/100
min_predictions <- inv.logit(t(min_predictions_actual)[1,])

max_dose <- max(SOK_data$ob_count)
max_dose_floor <- floor(max_dose/100)*100
max_dose_ceiling <- ceiling(max_dose/100)*100
max_predictions_floor <- filter(predictions, dose==max_dose_floor) %>% select(prediction)
max_predictions_ceiling <- filter(predictions, dose==max_dose_ceiling) %>% select(prediction)
max_predictions_actual <- max_predictions_floor + (max_predictions_ceiling-max_predictions_floor)*(max_dose-max_dose_floor)/100
max_predictions <- inv.logit(t(max_predictions_actual)[1,])

SNPV_differences <- max_predictions[1:8] - min_predictions[1:8]
MNPV_differences <- max_predictions[9:16] - min_predictions[9:16]

summary(SNPV_differences)
summary(MNPV_differences)


