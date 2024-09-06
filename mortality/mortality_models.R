require(tidyverse)
require(rstan)
require(pracma)
require(loo)
require(binom)
require(boot)
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


#### load fits if available
morphotype_and_tree <- readRDS("stan_fits/morphotype_and_tree.rds")
tree_only <- readRDS("stan_fits/tree_only.rds")
morphotype_only <- readRDS("stan_fits/morphotype_only.rds")
neither_morphotype_nor_tree <- readRDS("stan_fits/neither_morphotype_nor_tree.rds")
no_hierarchy <- readRDS("stan_fits/no_hierarchy.rds")



### set up stan
set.seed(683038)

options(mc.cores = 1)
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


## stan models

morphotype_and_tree <-
  stan(file="morphotype_and_tree.stan",
       data=list(N=N,H=H,I=I,J=J,
                 cid=cid,sid=sid,tid=tid,
                 x=x,y=y,total=total),
       chains=4,
       iter=4000,
       init=1,
       control = list(adapt_delta=0.99, max_treedepth=20))


tree_only <-
  stan(file="tree_only.stan",
       data=list(N=N,I=I,J=J,
                  sid=sid,tid=tid,
                 x=x,y=y,total=total),
       chains=4,
       iter=4000,
       init=1,
       control = list(adapt_delta=0.99, max_treedepth=20))


morphotype_only <-
  stan(file="morphotype_only.stan",
       data=list(N=N,H=H,I=I,
                 cid=cid,sid=sid,
                 x=x,y=y,total=total),
       chains=4,
       iter=4000,
       init=1,
       control = list(adapt_delta=0.99, max_treedepth=20))


neither_morphotype_nor_tree <-
  stan(file="neither_morphotype_nor_tree.stan",
       data=list(N=N,I=I,
                 sid=sid,
                 x=x,y=y,total=total),
       chains=4,
       iter=4000,
       init=1,
       control = list(adapt_delta=0.99, max_treedepth=20))


no_hierarchy <-  # this model only presented in the supplememnt
  stan(file="no_hierarchy.stan",
       data=list(N=N,I=I,J=J,
                 sid=sid,tid=tid,
                 x=x,y=y,total=total),
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
fit <- saveFit(no_hierarchy, "stan_fits/no_hierarchy.rds")







#### diagnostics

loo1 <- loo(morphotype_and_tree)
loo2 <- loo(tree_only)
loo3 <- loo(morphotype_only)
loo4 <- loo(neither_morphotype_nor_tree)
loo5 <- loo(no_hierarchy) # this model only presented in the supplememnt

loo_table <- loo_compare(loo1, loo2, loo3, loo4, loo5)
loo_table <- cbind(elpd=c(loo1$estimates[1,1],
                          loo2$estimates[1,1],
                          loo3$estimates[1,1],
                          loo4$estimates[1,1],
                          loo5$estimates[1,1])[order(order(rownames(loo_table)))],
                   loo_table[,1:2])
rownames(loo_table) <- c("M and T","T only","M only","Neither M nor T","No hierarchy")[order(order(rownames(loo_table)))]
loo_table



## mortality_model
## plot of best mortality model

fit_hier <- morphotype_and_tree # best model

y <- SOK_data$dose_response

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
  summarise(total_virus=sum(total_virus) + .5 * (sum(total_virus)==0), # for visual purposes on the logit scale, assumes treatments with 0 mortality would have had one dead larvae if sample size was doubled
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
  xlab("Dose (thousands of occlusion bodies)") + ylab("Proportion virus-killed (logistic scale)") +
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  theme(strip.background = element_blank(),
        panel.spacing= unit(1.5, "lines"),
        plot.margin = margin(0,20,0,10))



## avg_mortality
## plot of average mortality rate for morphotype-tree combos

grouped <- SOK_data %>%
  group_by(capsid, tree_sp) %>%
  summarise(y = sum(total_virus)/sum(total_n),
            ymin=binom.confint(sum(total_virus),sum(total_n),method="wilson")$lower,
            ymax=binom.confint(sum(total_virus),sum(total_n),method="wilson")$upper) %>% 
  as.data.frame()

ggplot(data = grouped) +
  geom_line(aes(x=capsid,y=y,color=tree_sp,group=tree_sp)) +
  geom_point(aes(x=capsid,y=y,color=tree_sp)) +
  geom_errorbar(width=.25,
                aes(x=capsid,ymin=ymin,ymax=ymax,color=tree_sp))+
  scale_y_continuous(limits = c(0,1),expand = c(0,0)) +
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  xlab("Morphotype") +
  ylab("Proportion virus-killed")








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

