require(deSolve)
require(tidyverse)
require(rstan)
require(ggpubr)
theme_set(theme_pubr())

# make nice-looking scientific notation for plots
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  l <- gsub("\\+","",l)
  parse(text=l)
}


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

tau_min = 6 # minimum number of days from infection to death
tau_max = 26 # maximum number of days from infection to death
K = tau_max-tau_min+1

f_SNPV_GR <- dgamma(tau_min:tau_max+.5,shape=alpha_SNPV_GR,scale=beta_SNPV_GR) /
  sum(dgamma(tau_min:tau_max+.5,shape=alpha_SNPV_GR,scale=beta_SNPV_GR))
f_SNPV_DO <- dgamma(tau_min:tau_max+.5,shape=alpha_SNPV_DO,scale=beta_SNPV_DO) /
  sum(dgamma(tau_min:tau_max+.5,shape=alpha_SNPV_DO,scale=beta_SNPV_DO))
f_MNPV_GR <- dgamma(tau_min:tau_max+.5,shape=alpha_MNPV_GR,scale=beta_MNPV_GR) /
  sum(dgamma(tau_min:tau_max+.5,shape=alpha_MNPV_GR,scale=beta_MNPV_GR))
f_MNPV_DO <- dgamma(tau_min:tau_max+.5,shape=alpha_MNPV_DO,scale=beta_MNPV_DO) /
  sum(dgamma(tau_min:tau_max+.5,shape=alpha_MNPV_DO,scale=beta_MNPV_DO))

f_same_dist <- dgamma(tau_min:tau_max+.5,shape=alpha_same_dist,scale=beta_same_dist) /
  sum(dgamma(tau_min:tau_max+.5,shape=alpha_same_dist,scale=beta_same_dist))



SOK_var_model <- function(t,y,p){
  S = y[1]
  R = y[2]
  V = y[3]
  with(as.list(p), {
    lags = matrix(0, nrow=K, ncol=3)
    for (tau in tau_min:tau_max) {
      if (t > tau + offset) {
        lags[tau-tau_min+1,] = lagvalue(t-tau-offset)
      }
    }
    lagS = lags[,1]
    lagV = lags[,3]
    
    
    dS.dt = -beta * S * V - d * S + lambda * R - d * S
    dR.dt = (1-m) * beta * sum(f * lagS * lagV) - lambda * R - d * R
    dV.dt = m * beta * sum(exp(-d * tau_min:tau_max) * f * ns * lagS * lagV) - delta * V
    return(list(c(dS.dt, dR.dt, dV.dt)))
  })
}

S0 <- 1e6 # starting population of susceptible larvae
V0 <- 1e8 # starting number of virions in environment
ns <- 1e6/(1+exp((16-tau_min:tau_max))) # within-host growth of virions during infection over time in days
beta <- 1e-11 # contact/transmission rate
delta <- 0 # assuming viral decay in environment is minimal over the time span being modeled
d = 1/35 # natural death/metamorphosis rate of larvae
lambda = 1/7 # rate at which larvae who have been unsuccessfully infected lose temporary immunity (roughly 1/(time between molts))

ts = seq(0,200,.1)
y0 = c(S0, 0, V0)


ms <- # mortality rates for different morphotype-tree combos
  SOK_data %>%
  group_by(capsid,tree_sp) %>%
  summarise(m=sum(total_virus)/sum(total_n)) %>%
  ungroup()


for (diff_variances in c(FALSE, TRUE)) {
  for (diff_means in c(FALSE, TRUE)) {
    
    # MNPV_GR
    f = if (diff_variances) f_MNPV_GR else f_same_dist
    m_MNPV_GR = ms %>% filter(capsid=="MNPV", tree_sp=="GR") %>% pull(m)
    p = list(beta=beta,ns=ns,f=f,delta=delta,d=d,m=m_MNPV_GR,lambda=lambda,
             offset=
               if (diff_means)
                 sum(f_MNPV_GR*tau_min:tau_max)-sum(f*tau_min:tau_max)
             else
               sum(f_same_dist*tau_min:tau_max)-sum(f*tau_min:tau_max))
    out_MNPV_GR = dede(y=y0, times=ts, func=SOK_var_model, parms=p)
    
    # MNPV_DO
    f = if (diff_variances) f_MNPV_DO else f_same_dist
    m_MNPV_DO = ms %>% filter(capsid=="MNPV", tree_sp=="DO") %>% pull(m)
    p = list(beta=beta,ns=ns,f=f,delta=delta,d=d,m=m_MNPV_DO,lambda=lambda,
             offset=
               if (diff_means)
                 sum(f_MNPV_DO*tau_min:tau_max)-sum(f*tau_min:tau_max)
             else
               sum(f_same_dist*tau_min:tau_max)-sum(f*tau_min:tau_max))
    out_MNPV_DO = dede(y=y0, times=ts, func=SOK_var_model, parms=p)
    
    # SNPV_GR
    f = if (diff_variances) f_SNPV_GR else f_same_dist
    m_SNPV_GR = ms %>% filter(capsid=="SNPV", tree_sp=="GR") %>% pull(m)
    p = list(beta=beta,ns=ns,f=f,delta=delta,d=d,m=m_SNPV_GR,lambda=lambda,
             offset=
               if (diff_means)
                 sum(f_SNPV_GR*tau_min:tau_max)-sum(f*tau_min:tau_max)
             else
               sum(f_same_dist*tau_min:tau_max)-sum(f*tau_min:tau_max))
    out_SNPV_GR = dede(y=y0, times=ts, func=SOK_var_model, parms=p)

    # SNPV_DO
    f = if (diff_variances) f_SNPV_DO else f_same_dist
    m_SNPV_DO = ms %>% filter(capsid=="SNPV", tree_sp=="DO") %>% pull(m)
    p = list(beta=beta,ns=ns,f=f,delta=delta,d=d,m=m_SNPV_DO,lambda=lambda,
             offset=
               if (diff_means)
                 sum(f_SNPV_DO*tau_min:tau_max)-sum(f*tau_min:tau_max)
             else
               sum(f_same_dist*tau_min:tau_max)-sum(f*tau_min:tau_max))
    out_SNPV_DO = dede(y=y0, times=ts, func=SOK_var_model, parms=p)
    
    
    data <-
      data.frame(t=rep(ts,4),
                 morphotype=rep(c("MNPV","SNPV"),each=2*length(ts)),
                 tree_sp=rep(factor(c("GR","DO","GR","DO"),levels=c("GR","DO")),each=length(ts)),
                 y=log10(c(out_MNPV_GR[,4],
                           out_MNPV_DO[,4],
                           out_SNPV_GR[,4],
                           out_SNPV_DO[,4])))
    
    ## SOK_var_same_dists, SOK_var_diff_dists in main text
    ## same_means_same_vars, diff_means_same_vars, same_means_diff_vars, diff_means_diff_vars in supplement
    print(
      ggplot(data) +
        geom_line(aes(x=t, y=y, group=interaction(morphotype,tree_sp), color=tree_sp, lty=morphotype)) +
        scale_x_continuous(name="Time since first infection (days)") +
        scale_y_continuous(name="Virus particles in environment", labels=function(l) parse(text=paste0("10^",l)),
                           limits=c(8,10.2), breaks=8:10, expand=expansion(c(.02,.02))) +
        scale_color_discrete(name="Tree", labels=c("Grand fir","Douglas fir")) +
        scale_linetype_discrete(name="Morphotype") +
        ggtitle(paste0(if (diff_means) "Different means" else "Same means",
                       if (diff_variances) ", different variances" else ", same variances")))
  }
}



### viral_yield in supplement
### plot of viral yield upon larvae death based on days from infection to death

ggplot(data=data.frame(days=seq(tau_min,tau_max,.1),yield=1e6/(1+exp((16-seq(tau_min,tau_max,.1)))))) +
  geom_line(aes(days,yield)) +
  xlab("Time from infection to death (days)") +
  scale_y_continuous(name="Number of virus particles released at death",labels=fancy_scientific)

