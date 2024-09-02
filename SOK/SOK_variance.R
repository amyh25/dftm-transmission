require(deSolve)
require(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tau_min = 6
tau_max = 28
K = tau_max-tau_min+1

SOK_full <- function(t,y,p){
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
    
    
    dS.dt = -beta * S * V - d * S# + lambda * R - d * S
    dR.dt = 0#(1-m) * beta * sum(f * lagS * lagV) - lambda * R - d * R
    dV.dt = m * beta * sum(exp(-d * tau_min:tau_max) * f * ns * lagS * lagV) - delta * V
    return(list(c(dS.dt, dR.dt, dV.dt)))
  })
}

S0 = 1e6
V0 = 1e8
ns = exp(6:28*.25)#rep(5e8, K)#1e9/(1+exp((-(6:28)+20)/3))#6:tau_max/(10+6:tau_max)*1e9#seq(from=1e9*6/7,by=1e9/7,length.out=K)
ns = ns / mean(ns) * 1e6
beta = 1e-11
delta = 0#1/100
d = 0#1/50
m = .2
lambda = 1
ts = seq(0,365,.1)
y0 = c(S0, 0, V0)



for (concavity in c(0,.25)) {
  for (diff_variances in c(FALSE, TRUE)) {
    for (diff_means in c(FALSE, TRUE)) {
      
      
      ts = seq(0,365,.1)
      
      #ns = exp(6:28*concavity)#rep(5e8, K)#1e9/(1+exp((-(6:28)+20)/3))#6:tau_max/(10+6:tau_max)*1e9#seq(from=1e9*6/7,by=1e9/7,length.out=K)
      ns <- 1e6/(1+.1*exp(20-(6:28)))#10^(8.75 + (1/seq(1/.014,1/.006,(1/.006-1/.014)/(tau_max-tau_min))-.006) / (.014-.006) * (8.25-8.75))
      #ns <- 1:K/(1:K+1)
      #ns = ns / mean(ns) * 1e6
      
      ms <- SOK_data %>% group_by(capsid,tree_sp) %>% summarise(m=sum(total_virus)/sum(total_n))
      #mgf <- function(f, x) sum(f * exp(x * tau_min:tau_max))
      f_MNPV_GR = MNPV_GR %>% group_by(numeric_day) %>% summarise(n=sum(value)) %>% pull(n)
      f_MNPV_GR = f_MNPV_GR/sum(f_MNPV_GR)
      f = f_MNPV_GR
      m_MNPV_GR = ms %>% filter(capsid=="MNPV", tree_sp=="GR") %>% pull(m)
      #gamma = -uniroot(function(x) mgf(f_MNPV_GR,x)-m_MNPV_GR, c(-2,0))$root - d
      p = list(beta=beta,ns=ns,f=f,delta=delta,d=d,m=m_MNPV_GR,lambda=lambda,
               offset=0)
      #y0 = c(S0, rep(0,K), 0, f_MNPV_GR*P0)
      out_MNPV_GR = dede(y=y0, times=ts, func=SOK_full, parms=p)
      
      f_MNPV_DO = MNPV_DO %>% group_by(numeric_day) %>% summarise(n=sum(value)) %>% pull(n)
      f_MNPV_DO = f_MNPV_DO/sum(f_MNPV_DO)
      f = if (diff_variances) f_MNPV_DO else f_MNPV_GR
      m_MNPV_DO = ms %>% filter(capsid=="MNPV", tree_sp=="DO") %>% pull(m)
      #gamma = -uniroot(function(x) mgf(f_MNPV_DO,x)-m_MNPV_DO, c(-2,0))$root - d
      p = list(beta=beta,ns=ns,f=f,delta=delta,d=d,m=m_MNPV_DO,lambda=lambda,
               offset=
                 if (diff_means)
                   sum(f_MNPV_DO*tau_min:tau_max)-sum(f*tau_min:tau_max)
               else
                 sum(f_MNPV_GR*tau_min:tau_max)-sum(f*tau_min:tau_max))
      #y0 = c(S0, rep(0,K), 0, f_MNPV_DO*P0)
      out_MNPV_DO = dede(y=y0, times=ts, func=SOK_full, parms=p)
      
      f_SNPV_GR = SNPV_GR %>% group_by(numeric_day) %>% summarise(n=sum(value)) %>% pull(n)
      f_SNPV_GR = f_SNPV_GR/sum(f_SNPV_GR)
      f = if (diff_variances) f_SNPV_GR else f_MNPV_GR
      m_SNPV_GR = ms %>% filter(capsid=="SNPV", tree_sp=="GR") %>% pull(m)
      #gamma = -uniroot(function(x) mgf(f_SNPV_GR,x)-m_SNPV_GR, c(-2,0))$root - d
      p = list(beta=beta,ns=ns,f=f,delta=delta,d=d,m=m_SNPV_GR,lambda=lambda,
               offset=
                 if (diff_means)
                   sum(f_SNPV_GR*tau_min:tau_max)-sum(f*tau_min:tau_max)
               else
                 sum(f_MNPV_GR*tau_min:tau_max)-sum(f*tau_min:tau_max))
      #y0 = c(S0, rep(0,K), 0, f_SNPV_GR*P0)
      out_SNPV_GR = dede(y=y0, times=ts, func=SOK_full, parms=p)
      
      f_SNPV_DO = SNPV_DO %>% group_by(numeric_day) %>% summarise(n=sum(value)) %>% pull(n)
      f_SNPV_DO = f_SNPV_DO/sum(f_SNPV_DO)
      f = if (diff_variances) f_SNPV_DO else f_MNPV_GR
      m_SNPV_DO = ms %>% filter(capsid=="SNPV", tree_sp=="DO") %>% pull(m)
      #gamma = -uniroot(function(x) mgf(f_SNPV_DO,x)-m_SNPV_DO, c(-2,0))$root - d
      p = list(beta=beta,ns=ns,f=f,delta=delta,d=d,m=m_SNPV_DO,lambda=lambda,
               offset=
                 if (diff_means)
                   sum(f_SNPV_DO*tau_min:tau_max)-sum(f*tau_min:tau_max)
               else
                 sum(f_MNPV_GR*tau_min:tau_max)-sum(f*tau_min:tau_max))
      #y0 = c(S0, rep(0,K), 0, f_SNPV_DO*P0)
      out_SNPV_DO = dede(y=y0, times=ts, func=SOK_full, parms=p)
      
      
      data <-
        data.frame(t=rep(ts,4),
                   morphotype=rep(c("MNPV","SNPV"),each=2*length(ts)),
                   tree_sp=rep(factor(c("GR","DO","GR","DO"),levels=c("GR","DO")),each=length(ts)),
                   y=log10(c(out_MNPV_GR[,4],
                             out_MNPV_DO[,4],
                             out_SNPV_GR[,4],
                             out_SNPV_DO[,4])))
      print(
        ggplot(data) +
          geom_line(aes(x=t, y=y, group=interaction(morphotype,tree_sp), color=tree_sp, lty=morphotype)) +
          scale_x_continuous(name="Time since first infection (days)", expand=expansion(c(0,0))) +
          scale_y_continuous(name="Virus particles in environment", labels=function(l) parse(text=paste0("10^",l)),
                             limits=if (concavity==0) c(8,13) else c(8,13), expand=expansion(c(.02,.02))) +
          scale_color_discrete(name="Tree", labels=c("Grand fir","Douglas fir")) +
          scale_linetype_discrete(name="Morphotype") +
          ggtitle(paste0("concavity = ", concavity,
                         if (diff_means) ", different means" else ", same means",
                         if (diff_variances) ", different variances" else ", same variances")))
    }
  }
}
