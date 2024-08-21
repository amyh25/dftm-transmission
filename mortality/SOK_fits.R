require(tidyverse)
require(ggpubr)
theme_set(theme_pubr())

### load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
SOK_data <- read.csv(file="SOK_data.csv",header=TRUE,sep=",") %>%
            filter(capsid != "none")
SOK_data$capsid <- factor(SOK_data$capsid, levels = c("SNPV","MNPV"))
SOK_data$tree_sp <- factor(SOK_data$tree_sp, levels = c("GR","DO"))


### tidy the data
days <- colnames(SOK_data)[9:31]

day_num <- function(str) {
  return(substring(str,5) %>% strtoi())
}

tidy <- SOK_data %>%
  pivot_longer(days, names_to="day")

tidy <- tidy %>% cbind(numeric_day = unlist(lapply(tidy$day, day_num)))

#tidy %>% group_by(capsid, tree_sp) %>% summarise(n=sum(value)) # total virus-kill in each group

avg_SOK <- tidy %>% 
  group_by(capsid, tree_sp) %>% 
  summarise(value = rep(numeric_day, value)) %>% 
  summarise(avg_SOK = mean(value), 
            sd = sd(value), 
            se = sd(value) / sqrt(n())) %>% 
  ungroup() %>% 
  mutate(y = avg_SOK, ymin = y - se, ymax = y + se)


MNPV_DO <- filter(tidy, capsid == "MNPV", tree_sp == "DO")
MNPV_GR <- filter(tidy, capsid == "MNPV", tree_sp == "GR")
SNPV_DO <- filter(tidy, capsid == "SNPV", tree_sp == "DO")
SNPV_GR <- filter(tidy, capsid == "SNPV", tree_sp == "GR")

MNPV_DO_days <- rep(MNPV_DO$numeric_day,MNPV_DO$value)
MNPV_GR_days <- rep(MNPV_GR$numeric_day,MNPV_GR$value)
SNPV_DO_days <- rep(SNPV_DO$numeric_day,SNPV_DO$value)
SNPV_GR_days <- rep(SNPV_GR$numeric_day,SNPV_GR$value)


## avg_SOK
## plot of average speed of kill for morphotype-tree combos

tree.labs <- c("Grand fir","Douglas fir")
names(tree.labs) <- c("GR","DO")

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
ggsave("../figures/avg_SOK.pdf", height = 4, width = 5) 





#likelihood function

logLHoodHtg<-function(data, par){
  alpha = par[1];
  beta  = par[2];
  
  logLHood = dgamma(data, shape=alpha, scale=beta, log=TRUE); #2. Type in the log-likelihood function for the pure death model
  return(-sum(logLHood));  #We are looking for the negative sum of the log likelihoods across experimental units, because optim is a minimizer
  
}


#MNPV_DO

hist(MNPV_DO_days, freq=FALSE, breaks = 15)

MNPV_DO_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = MNPV_DO_days);

xs <- seq(4,26,0.1)
y_MNPV_DO <- dgamma(xs, shape = MNPV_DO_OptOut$par[1],
                    scale = MNPV_DO_OptOut$par[2])# * length(MNPV_DO_days)

lines(xs,y_MNPV_DO)

MNPV_DO_LHood <- MNPV_DO_OptOut$value


#MNPV_GR

hist(MNPV_GR_days, freq=FALSE, breaks = 15)

MNPV_GR_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = MNPV_GR_days);

y_MNPV_GR <- dgamma(xs, shape = MNPV_GR_OptOut$par[1],
                    scale = MNPV_GR_OptOut$par[2])# * length(MNPV_GR_days)

lines(xs,y_MNPV_GR)

MNPV_GR_LHood <- MNPV_GR_OptOut$value



#SNPV_DO

hist(SNPV_DO_days, freq=FALSE, breaks = 15)

SNPV_DO_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = SNPV_DO_days);

y_SNPV_DO <- dgamma(xs, shape = SNPV_DO_OptOut$par[1],
                    scale = SNPV_DO_OptOut$par[2])# * length(SNPV_DO_days)

lines(xs,y_SNPV_DO)

SNPV_DO_LHood <- SNPV_DO_OptOut$value



#SNPV_GR

hist(SNPV_GR_days, freq=FALSE, breaks = 15)

SNPV_GR_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = SNPV_GR_days);

y_SNPV_GR <- dgamma(xs, shape = SNPV_GR_OptOut$par[1],
                    scale = SNPV_GR_OptOut$par[2])# * length(SNPV_GR_days)

lines(xs,y_SNPV_GR)

SNPV_GR_LHood <- SNPV_GR_OptOut$value







##grouping by morphotype

MNPV <- filter(tidy, capsid == "MNPV")
SNPV <- filter(tidy, capsid == "SNPV")

MNPV_days <- rep(MNPV$numeric_day, MNPV$value)
SNPV_days <- rep(SNPV$numeric_day, SNPV$value)


#MNPV

hist(MNPV_days, freq=FALSE, breaks = 15)

MNPV_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = MNPV_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = MNPV_OptOut$par[1], scale = MNPV_OptOut$par[2])

lines(xnums,ynums)

MNPV_LHood <- MNPV_OptOut$value



#SNPV

hist(SNPV_days, freq=FALSE, breaks = 15)

SNPV_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = SNPV_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = SNPV_OptOut$par[1], scale = SNPV_OptOut$par[2])

lines(xnums,ynums)

SNPV_LHood <- SNPV_OptOut$value



##grouping by tree

DO <- filter(tidy, tree_sp == "DO")
GR <- filter(tidy, tree_sp == "GR")

DO_days <- rep(DO$numeric_day, DO$value)
GR_days <- rep(GR$numeric_day, GR$value)


#DO

hist(DO_days, freq=FALSE, breaks = 15)

DO_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = DO_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = DO_OptOut$par[1], scale = DO_OptOut$par[2])

lines(xnums,ynums)

DO_LHood <- DO_OptOut$value



#GR

hist(GR_days, freq=FALSE, breaks = 15)

GR_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = GR_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = GR_OptOut$par[1], scale = GR_OptOut$par[2])

lines(xnums,ynums)

GR_LHood <- GR_OptOut$value



##grouping by everything

SOK_days <- rep(tidy$numeric_day, tidy$value)

hist(SOK_days, freq=FALSE, breaks = 15)

SOK_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = SOK_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = SOK_OptOut$par[1], scale = SOK_OptOut$par[2])

lines(xnums,ynums)

SOK_LHood <- SOK_OptOut$value



## AICcs
AICcCalc <- function(negLogLHood, K, N) {
  2*negLogLHood + 2*K + 2*K*(K+1)/(N-K-1)
}


# neither morphotype nor tree
C1 <- AICcCalc(SOK_LHood, 2*1, length(SOK_days))

# morphotype only
C2 <- AICcCalc(MNPV_LHood + SNPV_LHood, 2*2, length(SOK_days))

# tree only
C3 <- AICcCalc(DO_LHood + GR_LHood, 2*2, length(SOK_days))

# morphotype and tree
C4 <- AICcCalc(MNPV_DO_LHood + MNPV_GR_LHood + SNPV_DO_LHood + SNPV_GR_LHood, 2*4, length(SOK_days))

name <- c("C4","C2","C3","C1")
K <- c(2*4,2*2,2*2,2*1)
SOK_AICc <- c(C4,C2,C3,C1)
delta_AICc <- SOK_AICc - min(SOK_AICc)
weight <- exp(-delta_AICc/2)/sum(exp(-delta_AICc/2))
LHood <- -c(MNPV_DO_LHood + MNPV_GR_LHood + SNPV_DO_LHood + SNPV_GR_LHood,
            MNPV_LHood + SNPV_LHood,
            DO_LHood + GR_LHood,
            SOK_LHood)
SOK_table <- data.frame(name,LHood,K,SOK_AICc,delta_AICc,weight)


## SOK_model
## plot of best speed of kill distribution model

C4 <- data.frame(capsid=c(rep("SNPV",2*length(xs)),rep("MNPV",2*length(xs))),
                 tree_sp=rep(c(rep("GR",length(xs)),rep("DO",length(xs))),2),
                 x=rep(xs,4),
                 y=c(y_SNPV_GR,y_SNPV_DO,y_MNPV_GR,y_MNPV_DO))
C4$tree_sp <- factor(C4$tree_sp,levels=c("GR","DO"))

days <- data.frame(capsid=c(rep("SNPV",length(SNPV_days)),rep("MNPV",length(MNPV_days))),
                   tree_sp=c(rep("GR",length(SNPV_GR_days)),rep("DO",length(SNPV_DO_days)),
                             rep("GR",length(MNPV_GR_days)),rep("DO",length(MNPV_DO_days))),
                   day=c(SNPV_GR_days,SNPV_DO_days,MNPV_GR_days,MNPV_DO_days))
days$tree_sp <- factor(days$tree_sp,levels=c("GR","DO"))

tree.labs <- c("Grand fir","Douglas fir")
names(tree.labs) <- c("GR","DO")

ggplot() +
  geom_histogram(data=days,aes(x=day,y=stat(density),fill=tree_sp),
                 binwidth=1,boundary=0,color="black",size=.1) +
  geom_line(data=C4,aes(x=x,y=y)) +
  facet_grid(capsid~tree_sp,labeller=labeller(tree_sp=tree.labs)) +
  coord_cartesian(xlim=c(4,27),ylim=c(0,.2)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_discrete(name = "Tree", labels = c("Grand fir", "Douglas fir")) +
  xlab("Speed of kill (days)") +
  ylab("Proportion killed") +
  theme(strip.background = element_blank(),legend.position = "none")


## t-tests

t.test(MNPV_DO_days, MNPV_GR_days, alternative = "l")
t.test(SNPV_DO_days, SNPV_GR_days, alternative = "g")

