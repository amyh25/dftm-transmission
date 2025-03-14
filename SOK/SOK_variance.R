require(deSolve)
require(tidyverse)
require(rstan)
require(cowplot)
require(ggpubr)
theme_set(theme_pubr())


#### load simulation results if available
# load("sims_same_means_DO.RData")
# load("sims_diff_means_DO.RData")
# load("sims_same_means_GR.RData")
# load("sims_diff_means_GR.RData")


##### SNPV indexed by 1, MNPV indexed by 2 #####

### load fits from SOK_models.R
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
morphotype_and_tree <- readRDS("stan_fits/morphotype_and_tree.rds")
tree_only <- readRDS("stan_fits/tree_only.rds")

pars_diff_means <- summary(morphotype_and_tree)$summary[1:8,"mean"]
pars_same_mean <- summary(tree_only)$summary[1:4,"mean"]

alpha_SNPV_GR <- unname(pars_diff_means["alpha[1,1]"])
beta_SNPV_GR <- unname(pars_diff_means["beta[1,1]"])
delta_SNPV_GR <- 1 / beta_SNPV_GR

alpha_SNPV_DO <- unname(pars_diff_means["alpha[1,2]"])
beta_SNPV_DO <- unname(pars_diff_means["beta[1,2]"])
delta_SNPV_DO <- 1 / beta_SNPV_DO

alpha_MNPV_GR <- unname(pars_diff_means["alpha[2,1]"])
beta_MNPV_GR <- unname(pars_diff_means["beta[2,1]"])
delta_MNPV_GR <- 1 / beta_MNPV_GR

alpha_MNPV_DO <- unname(pars_diff_means["alpha[2,2]"])
beta_MNPV_DO <- unname(pars_diff_means["beta[2,2]"])
delta_MNPV_DO <- 1 / beta_MNPV_DO


alpha_GR <- unname(pars_same_mean["alpha[1]"])
beta_GR <- unname(pars_same_mean["beta[1]"])
delta_GR <- 1 / beta_GR

alpha_DO <- unname(pars_same_mean["alpha[2]"])
beta_DO <- unname(pars_same_mean["beta[2]"])
delta_DO <- 1 / beta_DO



within_season <- function(t, y, parms) {
  with(as.list(parms), {
    S = y[1]
    nu1 = y[2]
    nu2 = y[3]
    P1 = y[4]
    P2 = y[5]
    i1 = y[6]
    i2 = y[7]
    E1s = y[8:(k1+7)]
    E2s = y[(k1+8):(k1+k2+7)]
    
    dS = -P1 * nu1 * S - P2 * nu2 * S
    dnu1 = -P1 * nu1^2 * C1^2 - P2 * rho * C1 * C2 * nu1 * nu2
    dnu2 = -P2 * nu2^2 * C2^2 - P1 * rho * C1 * C2 * nu1 * nu2
    dP1 = k1 * delta1 * E1s[k1] - mu * P1
    dP2 = k2 * delta2 * E2s[k2] - mu * P2
    di1 = P1 * nu1 * S / N
    di2 = P2 * nu2 * S / N
    if (k1 > 1) {
      dE1s = c(P1 * nu1 * S, k1 * delta1 * E1s[1:(k1-1)]) - k1 * delta1 * E1s
    } else
      dE1s = P1 * nu1 * S - k1 * delta1 * E1s
    if (k2 > 1) {
      dE2s = c(P2 * nu2 * S, k2 * delta2 * E2s[1:(k2-1)]) - k2 * delta2 * E2s
    } else {
      dE2s = P2 * nu2 * S - k2 * delta2 * E2s
    }
    
    return(list(c(dS, dnu1, dnu2, dP1, dP2, di1, di2, dE1s, dE2s)))
  })
}

between_season <- function(y, parms) {
  with(as.list(parms), {
    nu1 = unname(y[2])
    nu2 = unname(y[3])
    i1 = unname(y[6])
    i2 = unname(y[7])
    
    return(c(
      N = N * (1 - i1 - i2) * r * (1 + s1 * nu1 + s2 * nu2),
      nu1 = (nu1 + s1 * (b1^2 * C1^2 +1) * nu1^2 + s2 * (rho * b1 * C1 * b2 * C2 + 1) * nu1 * nu2) / (1 + s1 * nu1 + s2 * nu2),
      nu2 = (nu2 + s2 * (b2^2 * C2^2 +1) * nu1^2 + s1 * (rho * b1 * C1 * b2 * C2 + 1) * nu1 * nu2) / (1 + s1 * nu1 + s2 * nu2),
      Z1 = phi1 * N * i1 + gamma * Z1 + 1e-5,
      Z2 = phi2 * N * i2 + gamma * Z2 + 1e-5
    ))
  })
}


ts <- seq(0,70,1)
parms_GR <- c(C1=3.7, C2=3.3, s1=1396, s2=60.7, mu=.55, r=.2, gamma=.5, b1=1, b2=1, phi1=10, phi2=10, rho=.101)
parms_DO <- c(C1=3.3, C2=2.5, s1=1.1e-5, s2=5.95, mu=.55, r=.2, gamma=.5, b1=1, b2=1, phi1=10, phi2=10, rho=.101)



## same means, on Douglas fir

sims_same_means_DO <- data.frame()

for (k1 in 1:30) {
  print(paste("k1 =",k1))
  
for (k2 in 1:30) {
  print(paste("  k2 =",k2))
  
  if (k2 == 1) {
    Ns <- 1
    nu1s <- .1
    nu2s <- .1
    Z1s <- 1
    Z2s <- 1
    i1s <- 0
    i2s <- 0
  } else {
    Ns <- unname(Ns[length(Ns)])
    nu1s <- unname(nu1s[length(nu1s)])
    nu2s <- unname(nu2s[length(nu2s)])
    Z1s <- unname(Z1s[length(Z1s)])
    Z2s <- unname(Z2s[length(Z2s)])
    i1s <- 0
    i2s <- 0
  }
  
  for (i in 2:ifelse(k2==1, 20000, 10000)) {
    y0 = c(Ns[i-1], nu1s[i-1], nu2s[i-1], Z1s[i-1], Z2s[i-1], 0, 0, rep(0, k1+k2))
    parms <- c(parms_DO, N=Ns[i-1], Z1=Z1s[i-1], Z2=Z2s[i-1], k1=k1, k2=k2, delta1=delta_DO, delta2=delta_DO)
    season <- ode(y0, ts, within_season, parms)
    
    i1s[i] <- season[length(ts),7]
    i2s[i] <- season[length(ts),8]
    
    beginning_of_next_season <- between_season(season[length(ts),-1], parms)
    Ns[i] <- beginning_of_next_season["N"]
    nu1s[i] <- beginning_of_next_season["nu1"]
    nu2s[i] <- beginning_of_next_season["nu2"]
    Z1s[i] <- beginning_of_next_season["Z1"]
    Z2s[i] <- beginning_of_next_season["Z2"]
  }
  
  sims_same_means_DO <- rbind(sims_same_means_DO,
                              data.frame(k1 = k1,
                                         k2 = k2,
                                         gen = 1:10000,
                                         N = Ns[(length(Ns)-9999):length(Ns)],
                                         nu1 = nu1s[(length(Ns)-9999):length(Ns)],
                                         nu2 = nu2s[(length(Ns)-9999):length(Ns)],
                                         Z1 = Z1s[(length(Ns)-9999):length(Ns)],
                                         Z2 = Z2s[(length(Ns)-9999):length(Ns)],
                                         i1 = i1s[(length(Ns)-9999):length(Ns)],
                                         i2 = i2s[(length(Ns)-9999):length(Ns)]))
}
  
  # saves data incrementally in case of R crashing
  if (k1 == 30) {
    save(sims_same_means_DO, file="sims_same_means_DO.RData")
  } else if (k1 %% 5 == 0) {
    save(sims_same_means_DO, file=paste0("sims_same_means_DO_",k1 %/% 5,".RData"))
  }
}



## different means, on Douglas fir

sims_diff_means_DO <- data.frame()

for (k1 in 1:30) {
  print(paste("k1 =",k1))
  
  for (k2 in 1:30) {
    print(paste("  k2 =",k2))
    
    if (k2 == 1) {
      Ns <- 1
      nu1s <- .1
      nu2s <- .1
      Z1s <- 1
      Z2s <- 1
      i1s <- 0
      i2s <- 0
    } else {
      Ns <- unname(Ns[length(Ns)])
      nu1s <- unname(nu1s[length(nu1s)])
      nu2s <- unname(nu2s[length(nu2s)])
      Z1s <- unname(Z1s[length(Z1s)])
      Z2s <- unname(Z2s[length(Z2s)])
      i1s <- 0
      i2s <- 0
    }
    
    for (i in 2:ifelse(k2==1, 20000, 10000)) {
      y0 = c(Ns[i-1], nu1s[i-1], nu2s[i-1], Z1s[i-1], Z2s[i-1], 0, 0, rep(0, k1+k2))
      parms <- c(parms_DO, N=Ns[i-1], Z1=Z1s[i-1], Z2=Z2s[i-1], k1=k1, k2=k2, delta1=delta_SNPV_DO, delta2=delta_MNPV_DO)
      season <- ode(y0, ts, within_season, parms)
      
      i1s[i] <- season[length(ts),7]
      i2s[i] <- season[length(ts),8]
      
      beginning_of_next_season <- between_season(season[length(ts),-1], parms)
      Ns[i] <- beginning_of_next_season["N"]
      nu1s[i] <- beginning_of_next_season["nu1"]
      nu2s[i] <- beginning_of_next_season["nu2"]
      Z1s[i] <- beginning_of_next_season["Z1"]
      Z2s[i] <- beginning_of_next_season["Z2"]
    }
    
    sims_diff_means_DO <- rbind(sims_diff_means_DO,
                                data.frame(k1 = k1,
                                           k2 = k2,
                                           gen = 1:10000,
                                           N = Ns[(length(Ns)-9999):length(Ns)],
                                           nu1 = nu1s[(length(Ns)-9999):length(Ns)],
                                           nu2 = nu2s[(length(Ns)-9999):length(Ns)],
                                           Z1 = Z1s[(length(Ns)-9999):length(Ns)],
                                           Z2 = Z2s[(length(Ns)-9999):length(Ns)],
                                           i1 = i1s[(length(Ns)-9999):length(Ns)],
                                           i2 = i2s[(length(Ns)-9999):length(Ns)]))
  }
  
  if (k1 == 30) {
    save(sims_diff_means_DO, file="sims_diff_means_DO.RData")
  } else if (k1 %% 5 == 0) {
    save(sims_diff_means_DO, file=paste0("sims_diff_means_DO_",k1 %/% 5,".RData"))
  }
}




## freq_SNPV_DO_same_means
## frequency of SNPV on Douglas fir for different k1 and k2 values,
## with SOK distributions coming from tree_only SOK model
sims_same_means_DO %>%
  filter(gen > 2000,
         i1 + i2 > .1) %>%
  group_by(k1, k2) %>%
  summarise(freq = mean(i1 / (i1+i2))) %>%
  ggplot() +
  geom_tile(aes(x=k1,y=k2,fill=freq,color=freq)) +
  geom_point(data=data.frame(k1=alpha_DO, k2=alpha_DO), aes(x=k1,y=k2), shape=4, size=3) + 
  coord_fixed() +
  scale_x_continuous(name="SNPV speed of kill CV",
                     expand=c(0,0),breaks=c(1,10,20,30),
                     labels=c(1, expression(1/sqrt(10)), expression(1/sqrt(20)), expression(1/sqrt(30)))) + 
  scale_y_continuous(name="MNPV speed of kill CV",
                     expand=c(0,0),breaks=c(1,10,20,30),
                     labels=c(1, expression(1/sqrt(10)), expression(1/sqrt(20)), expression(1/sqrt(30)))) +
  scale_fill_gradientn(colors=rev(rainbow(7))[-1], name="Average\nfrequency\nof SNPV",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_color_gradientn(colors=rev(rainbow(7))[-1], guide = "none") +
  ggtitle("Same mean speeds of kill") +
  theme(strip.background=element_blank(),
        legend.position="right",
        legend.title=element_text(size = 10),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        strip.text.x=element_text(size = 12),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(hjust = 0.5))


## freq_SNPV_DO_diff_means
## frequency of SNPV on Douglas fir for different k1 and k2 values,
## with SOK distributions coming from morphotype_and_tree SOK model
sims_diff_means_DO %>%
  filter(gen > 2000,
         i1 + i2 > .1) %>%
  group_by(k1, k2) %>%
  summarise(freq = mean(i1 / (i1+i2))) %>%
  ggplot() +
  geom_tile(aes(x=k1,y=k2,fill=freq,color=freq)) +
  geom_point(data=data.frame(k1=alpha_SNPV_DO, k2=alpha_MNPV_DO), aes(x=k1,y=k2), shape=4, size=3) + 
  coord_fixed() +
  scale_x_continuous(name="SNPV speed of kill CV",
                     expand=c(0,0),breaks=c(1,10,20,30),
                     labels=c(1, expression(1/sqrt(10)), expression(1/sqrt(20)), expression(1/sqrt(30)))) + 
  scale_y_continuous(name="MNPV speed of kill CV",
                     expand=c(0,0),breaks=c(1,10,20,30),
                     labels=c(1, expression(1/sqrt(10)), expression(1/sqrt(20)), expression(1/sqrt(30)))) +
  scale_fill_gradientn(colors=rev(rainbow(7))[-1], name="Average\nfrequency\nof SNPV",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_color_gradientn(colors=rev(rainbow(7))[-1], guide = "none") +
  ggtitle("Different mean speeds of kill") +
  theme(strip.background=element_blank(),
        legend.position="right",
        legend.title=element_text(size = 10),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        strip.text.x=element_text(size = 12),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(hjust = 0.5))












##### showing that k1, k2, and same/different means have no effect on grand fir

sims_same_means_GR <- data.frame()

for (k1 in 1:30) {
  print(paste("k1 =",k1))
  
  for (k2 in 1:30) {
    print(paste("  k2 =",k2))
    
    if (k2 == 1) {
      Ns <- 1
      nu1s <- .1
      nu2s <- .1
      Z1s <- 1
      Z2s <- 1
      i1s <- 0
      i2s <- 0
    } else {
      Ns <- unname(Ns[length(Ns)])
      nu1s <- unname(nu1s[length(nu1s)])
      nu2s <- unname(nu2s[length(nu2s)])
      Z1s <- unname(Z1s[length(Z1s)])
      Z2s <- unname(Z2s[length(Z2s)])
      i1s <- 0
      i2s <- 0
    }
    
    for (i in 2:ifelse(k2==1, 5000, 2500)) {
      y0 = c(Ns[i-1], nu1s[i-1], nu2s[i-1], Z1s[i-1], Z2s[i-1], 0, 0, rep(0, k1+k2))
      parms <- c(parms_GR, N=Ns[i-1], Z1=Z1s[i-1], Z2=Z2s[i-1], k1=k1, k2=k2, delta1=delta_GR, delta2=delta_GR)
      season <- ode(y0, ts, within_season, parms)
      
      i1s[i] <- season[length(ts),7]
      i2s[i] <- season[length(ts),8]
      
      beginning_of_next_season <- between_season(season[length(ts),-1], parms)
      Ns[i] <- beginning_of_next_season["N"]
      nu1s[i] <- beginning_of_next_season["nu1"]
      nu2s[i] <- beginning_of_next_season["nu2"]
      Z1s[i] <- beginning_of_next_season["Z1"]
      Z2s[i] <- beginning_of_next_season["Z2"]
    }
    
    sims_same_means_GR <- rbind(sims_same_means_GR,
                                data.frame(k1 = k1,
                                           k2 = k2,
                                           gen = 1:2500,
                                           N = Ns[(length(Ns)-2499):length(Ns)],
                                           nu1 = nu1s[(length(Ns)-2499):length(Ns)],
                                           nu2 = nu2s[(length(Ns)-2499):length(Ns)],
                                           Z1 = Z1s[(length(Ns)-2499):length(Ns)],
                                           Z2 = Z2s[(length(Ns)-2499):length(Ns)],
                                           i1 = i1s[(length(Ns)-2499):length(Ns)],
                                           i2 = i2s[(length(Ns)-2499):length(Ns)]))
  }
  
  if (k1 == 30) {
    save(sims_same_means_GR, file="sims_same_means_GR.RData")
  } else if (k1 %% 5 == 0) {
    save(sims_same_means_GR, file=paste0("sims_same_means_GR_",k1 %/% 5,".RData"))
  }
}



sims_diff_means_GR <- data.frame()

for (k1 in 1:30) {
  print(paste("k1 =",k1))
  
  for (k2 in 1:30) {
    print(paste("  k2 =",k2))
    
    if (k2 == 1) {
      Ns <- 1
      nu1s <- .1
      nu2s <- .1
      Z1s <- 1
      Z2s <- 1
      i1s <- 0
      i2s <- 0
    } else {
      Ns <- unname(Ns[length(Ns)])
      nu1s <- unname(nu1s[length(nu1s)])
      nu2s <- unname(nu2s[length(nu2s)])
      Z1s <- unname(Z1s[length(Z1s)])
      Z2s <- unname(Z2s[length(Z2s)])
      i1s <- 0
      i2s <- 0
    }
    
    for (i in 2:ifelse(k2==1, 5000, 2500)) {
      y0 = c(Ns[i-1], nu1s[i-1], nu2s[i-1], Z1s[i-1], Z2s[i-1], 0, 0, rep(0, k1+k2))
      parms <- c(parms_GR, N=Ns[i-1], Z1=Z1s[i-1], Z2=Z2s[i-1], k1=k1, k2=k2, delta1=delta_SNPV_GR, delta2=delta_MNPV_GR)
      season <- ode(y0, ts, within_season, parms)
      
      i1s[i] <- season[length(ts),7]
      i2s[i] <- season[length(ts),8]
      
      beginning_of_next_season <- between_season(season[length(ts),-1], parms)
      Ns[i] <- beginning_of_next_season["N"]
      nu1s[i] <- beginning_of_next_season["nu1"]
      nu2s[i] <- beginning_of_next_season["nu2"]
      Z1s[i] <- beginning_of_next_season["Z1"]
      Z2s[i] <- beginning_of_next_season["Z2"]
    }
    
    sims_diff_means_GR <- rbind(sims_diff_means_GR,
                                data.frame(k1 = k1,
                                           k2 = k2,
                                           gen = 1:2500,
                                           N = Ns[(length(Ns)-2499):length(Ns)],
                                           nu1 = nu1s[(length(Ns)-2499):length(Ns)],
                                           nu2 = nu2s[(length(Ns)-2499):length(Ns)],
                                           Z1 = Z1s[(length(Ns)-2499):length(Ns)],
                                           Z2 = Z2s[(length(Ns)-2499):length(Ns)],
                                           i1 = i1s[(length(Ns)-2499):length(Ns)],
                                           i2 = i2s[(length(Ns)-2499):length(Ns)]))
  }
  
  
  
  if (k1 == 30) {
    save(sims_diff_means_GR, file="sims_diff_means_GR.RData")
  } else if (k1 %% 5 == 0) {
    save(sims_diff_means_GR, file=paste0("sims_diff_means_GR_",k1 %/% 5,".RData"))
  }
}


## freq_SNPV_GR_same_means
## frequency of SNPV on grand fir for different k1 and k2 values,
## with SOK distributions coming from tree_only SOK model
sims_same_means_GR %>%
  filter(gen > 2000,
         i1 + i2 > .1) %>%
  group_by(k1, k2) %>%
  summarise(freq = mean(i1 / (i1+i2))) %>%
  ggplot() +
  geom_tile(aes(x=k1,y=k2,fill=freq,color=freq)) +
  geom_point(data=data.frame(k1=alpha_GR, k2=alpha_GR), aes(x=k1,y=k2), shape=4, size=3) + 
  coord_fixed() +
  scale_x_continuous(name="SNPV speed of kill CV",
                     expand=c(0,0),breaks=c(1,10,20,30),
                     labels=c(1, expression(1/sqrt(10)), expression(1/sqrt(20)), expression(1/sqrt(30)))) + 
  scale_y_continuous(name="MNPV speed of kill CV",
                     expand=c(0,0),breaks=c(1,10,20,30),
                     labels=c(1, expression(1/sqrt(10)), expression(1/sqrt(20)), expression(1/sqrt(30)))) +
  scale_fill_gradientn(colors=rev(rainbow(7))[-1], name="Average\nfrequency\nof SNPV",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_color_gradientn(colors=rev(rainbow(7))[-1], guide = "none") +
  ggtitle("Same mean speeds of kill") +
  theme(strip.background=element_blank(),
        legend.position="right",
        legend.title=element_text(size = 10),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        strip.text.x=element_text(size = 12),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(hjust = 0.5))


## freq_SNPV_GR_diff_means
## frequency of SNPV on grand fir for different k1 and k2 values,
## with SOK distributions coming from morphotype_and_tree SOK model
sims_diff_means_GR %>%
  filter(gen > 2000,
         i1 + i2 > .1) %>%
  group_by(k1, k2) %>%
  summarise(freq = mean(i1 / (i1+i2))) %>%
  ggplot() +
  geom_tile(aes(x=k1,y=k2,fill=freq,color=freq)) +
  geom_point(data=data.frame(k1=alpha_SNPV_GR, k2=alpha_MNPV_GR), aes(x=k1,y=k2), shape=4, size=3) + 
  coord_fixed() +
  scale_x_continuous(name="SNPV speed of kill CV",
                     expand=c(0,0),breaks=c(1,10,20,30),
                     labels=c(1, expression(1/sqrt(10)), expression(1/sqrt(20)), expression(1/sqrt(30)))) + 
  scale_y_continuous(name="MNPV speed of kill CV",
                     expand=c(0,0),breaks=c(1,10,20,30),
                     labels=c(1, expression(1/sqrt(10)), expression(1/sqrt(20)), expression(1/sqrt(30)))) +
  scale_fill_gradientn(colors=rev(rainbow(7))[-1], name="Average\nfrequency\nof SNPV",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_color_gradientn(colors=rev(rainbow(7))[-1], guide = "none") +
  ggtitle("Different mean speeds of kill") +
  theme(strip.background=element_blank(),
        legend.position="right",
        legend.title=element_text(size = 10),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        strip.text.x=element_text(size = 12),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(hjust = 0.5))



