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
require(stringr)
theme_set(theme_pubr())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv(file="manual_avoid_effect.csv",
                 header=TRUE, sep=",", stringsAsFactors=FALSE)
data$morphotype <- factor(data$morphotype, levels = c("SNPV","MNPV"))
data$tree_sp <- factor(data$tree_sp, levels = c("GR","DO"))
data$isolate <- factor(data$isolate,level=c("COL","LOV","DRY","TAM","CTRL"))
data <- data %>% filter(control != "NA", treatment != "NA")
data$control[data$control < 0] <- 0
data$treatment[data$treatment < 0] <- 0
data <- data %>% cbind(D = data$control - data$treatment)

data <- data %>% cbind(list(D_hat = data$D))
data$D_hat[data$tree_sp=="DO"] <-
  data$D[data$tree_sp=="DO"] - mean(data$D[data$tree_sp=="DO" & data$isolate=="CTRL"])
data$D_hat[data$tree_sp=="GR"] <-
  data$D[data$tree_sp=="GR"] - mean(data$D[data$tree_sp=="GR" & data$isolate=="CTRL"])


### load fits
# might have to manually open the files first to decompress them
fit_hier_tree_int <- readRDS("stan_fits/fit_hier_tree_int.rds")
fit_hier_tree <- readRDS("stan_fits/fit_hier_tree.rds")
fit_hier <- readRDS("stan_fits/fit_hier.rds")
fit_tree <- readRDS("stan_fits/fit_tree.rds")
fit_tree_only <- readRDS("stan_fits/fit_tree_only.rds")
fit_treatment <- readRDS("stan_fits/fit_treatment.rds")
fit_intercept <- readRDS("stan_fits/fit_intercept.rds")

bhms <- c(fit_hier_tree_int,
          fit_hier_tree,fit_hier,fit_tree,fit_tree_only,fit_treatment,fit_intercept)
name <- c("Hier+tree+interaction","Hier+tree","Hier","Treatment+tree","Tree","Treatment","Intercept")
looic <- sapply(bhms,function(model) loo(model)$estimates[3])
delta_looic <- looic - min(looic)
weight <- exp(-looic/2) / sum(exp(-looic/2))
p_looic <- sapply(bhms,function(model) loo(model)$estimates[2])
looic_table <- data.frame(name,looic,delta_looic,weight,p_looic)



# best model
fit <- fit_hier_tree_int
estimates <- summary(fit)$summary[,"mean"]
betas <- estimates[str_detect(names(estimates),pattern="^beta")]

# input defined at the beginning of avoidance/stan_model.R
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
