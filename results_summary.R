#### code to make results_summary figure

require(cowplot)

### load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# avoidance data
avoidance_data <- read.csv(file="data/manual_avoid_effect.csv",
                  header=TRUE, sep=",", stringsAsFactors=FALSE)
avoidance_data$morphotype <- factor(avoidance_data$morphotype, levels = c("SNPV","MNPV"))
avoidance_data <- avoidance_data %>% filter(control != "NA", treatment != "NA")
avoidance_data$control[avoidance_data$control < 0] <- 0
avoidance_data$treatment[avoidance_data$treatment < 0] <- 0
avoidance_data <- avoidance_data %>% cbind(D = avoidance_data$control - avoidance_data$treatment)
avoidance_data <- avoidance_data %>% cbind(list(D_hat = avoidance_data$D))
avoidance_data$D_hat[avoidance_data$tree_sp=="DO"] <-
  avoidance_data$D[avoidance_data$tree_sp=="DO"] - mean(avoidance_data$D[avoidance_data$tree_sp=="DO" & avoidance_data$isolate=="CTRL"])
avoidance_data$D_hat[avoidance_data$tree_sp=="GR"] <-
  avoidance_data$D[avoidance_data$tree_sp=="GR"] - mean(avoidance_data$D[avoidance_data$tree_sp=="GR" & avoidance_data$isolate=="CTRL"])
avoidance_data <- avoidance_data %>% mutate(avoid_effect = D_hat, type = case_when(isolate == "CTRL" ~ "CTRL",
                                                                 TRUE ~ "TREAT"))
avoidance_data$isolate <- factor(avoidance_data$isolate, levels = c("CTRL", "COL", "LOV", "DRY", "TAM"))
avoidance_data$tree_sp <- factor(avoidance_data$tree_sp, levels=c("GR", "DO"))


# mortality/speed of kill data
SOK_data <- read.csv(file="data/SOK_data.csv",
                     header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
  filter(capsid != "none")
SOK_data$strain <- factor(SOK_data$strain,levels=c("COL","CUB","LOV","LST","DRY","KLP","TAM","TMB"))
SOK_data$capsid <- factor(SOK_data$capsid,levels=c("SNPV","MNPV"))
SOK_data$tree_sp <- factor(SOK_data$tree_sp,levels=c("GR","DO"))
days <- colnames(SOK_data)[9:31]
day_num <- function(str) {
  return(substring(str,5) %>% strtoi())
}
tidy <- SOK_data %>%
  pivot_longer(days, names_to="day")
tidy <- tidy %>% cbind(numeric_day = unlist(lapply(tidy$day, day_num)))
tidy <- tidy[rep(1:nrow(tidy), times=tidy$value),]




D_bar <- data.frame(strain=rep(factor(c("SNPV","MPNV"),
                                      levels=c("SNPV","MPNV")),2),
                    tree_sp=rep(factor(c("GF","DF"),
                                       levels=c("GF","DF")),each=2),
                    y=avoidance_data %>% group_by(tree_sp,morphotype) %>% summarise(y=mean(D_hat)) %>% pull(y))

P_death <- data.frame(strain=rep(factor(c("SNPV","MPNV"),
                                        levels=c("SNPV","MPNV")),2),
                      tree_sp=rep(factor(c("GF","DF"),
                                         levels=c("GF","DF")),each=2),
                      y=SOK_data %>% group_by(tree_sp, capsid) %>% summarise(y = sum(total_virus)/sum(total_n)) %>% pull(y))

SOK_mean <- data.frame(strain=rep(factor(c("SNPV","MPNV"),
                                         levels=c("SNPV","MPNV")),2),
                      tree_sp=rep(factor(c("GF","DF"),
                                         levels=c("GF","DF")),each=2),
                      y=tidy %>% group_by(tree_sp, capsid) %>% summarise(y = mean(numeric_day)) %>% pull(y))

SOK_CV <- data.frame(strain=rep(factor(c("SNPV","MPNV"),
                                        levels=c("SNPV","MPNV")),2),
                       tree_sp=rep(factor(c("GF","DF"),
                                          levels=c("GF","DF")),each=2),
                       y=tidy %>% group_by(tree_sp, capsid) %>% summarise(y = sd(numeric_day) / mean(numeric_day)) %>% pull(y))





p1 <-
  D_bar %>%
  ggplot() +
  geom_point(aes(x=tree_sp,y=y,color=strain)) +
  geom_line(aes(x=tree_sp,y=y,color=strain,group=strain)) +
  scale_color_manual(name="Morphotype",values=c("goldenrod3","mediumpurple3")) +
  scale_x_discrete(name="Tree species") +
  ylab(expression("Avoidance metric, " * hat(D))) +
  theme(axis.text.y=element_text(size=10))

p2 <-
  P_death %>%
  ggplot() +
  geom_point(aes(x=tree_sp,y=y,color=strain)) +
  geom_line(aes(x=tree_sp,y=y,color=strain,group=strain)) +
  scale_color_manual(name="Morphotype",values=c("goldenrod3","mediumpurple3"),
                       labels=c("Grand fir (GF)","Douglas-fir (DF)")) +
  scale_x_discrete(name="Tree species") +
  ylab("Probability of death\ngiven infection") +
  theme(axis.text.y=element_text(size=10),
        legend.position="none")

p3 <-
  SOK_mean %>%
  ggplot() +
  geom_point(aes(x=tree_sp,y=y,color=strain)) +
  geom_line(aes(x=tree_sp,y=y,color=strain,group=strain)) +
  scale_color_manual(values=c("goldenrod3","mediumpurple3")) +
  scale_x_discrete(name="Tree species") +
  ylab("Mean speed of kill") +
  theme(axis.text.y=element_text(size=10),
        legend.position="none")

p4 <-
  SOK_CV %>%
  ggplot() +
  geom_point(aes(x=tree_sp,y=y,color=strain)) +
  geom_line(aes(x=tree_sp,y=y,color=strain,group=strain)) +
  scale_color_manual(values=c("goldenrod3","mediumpurple3")) +
  scale_x_discrete(name="Tree species") +
  ylab("Speed of kill\ncoefficient of variation") +
  theme(axis.text.y=element_text(size=10),
        legend.position="none")
  

grobs <- ggplotGrob(p1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

pg <- plot_grid(p1+theme(legend.position="none"),
                p2,p3,p4,align="v",axis="tbrl",ncol=4)

plot_grid(pg,legend,ncol=1,rel_heights=c(1,.1))


