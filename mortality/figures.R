require(tidyverse)
require(AICcmodavg)
require(boot)
require(binom)
require(gmodels)
require(ggstance)
require(scales)
require(ggpubr)
require(Hmisc)
theme_set(theme_pubr())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

SOK_data <- read.csv(file="SOK_data.csv",header=TRUE,sep=",") %>%
            filter(capsid != "none")
SOK_data$strain <- factor(SOK_data$strain,levels=c("COL","CUB","LOV","LST","DRY","KLP","TAM","TMB"))
SOK_data$capsid <- factor(SOK_data$capsid,levels=c("SNPV","MNPV"))
SOK_data$tree_sp <- factor(SOK_data$tree_sp,levels=c("GR","DO"))

dose_var <- SOK_data$ob_count



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



## avg_SOK
## plot of average speed of kill for morphotype-tree combos
## see SOK_fits.R



## mortality_model
## plot of best mortality model
## see stan_diagonistics.R



## SOK_model
## plot of best speed of kill distribution model
## see SOK_fits.R


