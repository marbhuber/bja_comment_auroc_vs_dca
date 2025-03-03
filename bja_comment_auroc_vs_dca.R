############################################################################################################
# Beyond the area under the curve:
# the benefit of a decision curve analysis for evaluating the added value of a biomarker for decision-making
# Comment on Roshanov et al. 2025 (https://doi.org/10.1016/j.bja.2024.10.039)
#
# Markus Huber / 03 March 2025
#############################################################################################################

#######
# setup

rm(list=ls())
library(tidyverse)
library(simstudy)
library(data.table)
library(dcurves)
library(givitiR)
library(CalibrationCurves)
library(cowplot)

# define data frames

df.dca.treat     <- c()
df.dca           <- c()
df.relative      <- c()
df.absolute      <- c()
df.auroc         <- c()

N = 35815 # number of patients

# N=1e3
# one random seed should be set initially prior to the generation of data
set.seed(1234)

###########
# bootstrap

n.boot = 1e3
for (myboot in 1:n.boot){

  print(myboot)
  
  #################################################################################################
  # Patient and surgery characteristics with eGFR and non-linear eGFR
  # create one dataset, then fit three logistic regression models with varying degree of complexity
  
  coefs1 <- c(0.2,-0.08,0.08)
  myprev = 0.13
  myauc  = 0.77
  
  rm(d1,d2,d3,C1,d1a,dd)
  d1  <- defData(varname = "x1", formula = 0, variance = 1)
  d2  <- defData(d1,varname = "x2", formula = 0, variance = 1)
  d3  <- defData(d2,varname = "x3", formula = "x2*x2")
  C1  <- logisticCoefs(defCovar = d3, coefs = coefs1, popPrev = myprev, auc=myauc)
  d1a <- defData(d3, varname = "y",formula = "t(..C1) %*% c(1, x1, x2,x3)",dist = "binary", link = "logit")
  dd  <- genData(N, d1a) %>% as.data.frame()
  
  ############################################
  # model: patient and surgery characteristics

  rm(fit,df.fit)
  fit    <- glm(y ~ x1, data = dd,family=binomial())
  df.fit <- data.frame(y=fit$y,prob = fit$fitted.values)

  df.auroc <- rbind(df.auroc,data.frame(
    boot = myboot,
    type = "x1",
    auc  = as.numeric(pROC::auc(pROC::roc(df.fit$y,df.fit$prob)))
    ))
  
  mydca0       <- dca(y ~ prob, df.fit, thresholds = seq(0, 1, by = 0.02))
  df.dca.treat <- rbind(df.dca.treat,mutate(filter(mydca0$dca,variable!="prob")))

  ######################################################
  # model: patient and surgery characteristics with eGFR
  
  rm(fit,df.fit)
  fit    <- glm(y ~ x1+x2, data = dd,family=binomial())
  df.fit <- data.frame(y=fit$y,prob = fit$fitted.values)
  
  mydca1 <- dca(y ~ prob, df.fit, thresholds = seq(0, 1, by = 0.02))
  df.auroc <- rbind(df.auroc,data.frame(
    boot = myboot,
    type = "x1+x2",
    auc  = as.numeric(pROC::auc(pROC::roc(df.fit$y,df.fit$prob)))))
  
  ##########################################################################
  # model: patient and surgery characteristics with eGFR and non-linear eGFR
  
  rm(fit,df.fit)
  fit <- glm(y ~ x1+x2+I(x2^2), data = dd,family=binomial())
  df.fit <- data.frame(y=fit$y,prob = fit$fitted.values)
  mydca2 <- dca(y ~ prob, df.fit, thresholds = seq(0, 1, by = 0.02))

  df.auroc <- rbind(df.auroc,data.frame(
  boot = myboot,
  type = "x1+x2+I(x2)",
  auc  = as.numeric(pROC::auc(pROC::roc(df.fit$y,df.fit$prob)))))

  # data frame with net benefits
  df.dca.dummy <- rbind(
    mutate(filter(mydca0$dca,variable=="prob"),type="x1"),
    mutate(filter(mydca1$dca,variable=="prob"),type="x1+x2"),
    mutate(filter(mydca2$dca,variable=="prob"),type="x1+x2+I(x2)")) %>%
    mutate(boot=myboot)

  df.dca <- rbind(df.dca,df.dca.dummy)

  # gains per 1'000 patients
  absolute <-  rbind(
  mutate(filter(mydca0$dca,variable=="prob"),type="x1"),
  mutate(filter(mydca1$dca,variable=="prob"),type="x2"),
  mutate(filter(mydca2$dca,variable=="prob"),type="x3"))  %>% 
  mutate(net_benefit=1000*net_benefit) %>% 
  pivot_wider(
    names_from=type,values_from=net_benefit
  ) %>% 
  group_by(
    threshold
  ) %>% 
  fill(x1,.direction = "updown") %>% 
  fill(x2,.direction = "updown") %>% 
  fill(x3,.direction = "updown") %>% 
  select(threshold,x1,x2,x3) %>% 
  unique() %>% 
  mutate(
    boot = myboot
  )

  df.absolute <- rbind(df.absolute,absolute)

  # relative gains per 1'000 patients
  relative <-  rbind(
  mutate(filter(mydca0$dca,variable=="prob"),type="x1"),
  mutate(filter(mydca1$dca,variable=="prob"),type="x2"),
  mutate(filter(mydca2$dca,variable=="prob"),type="x3")) %>% 
  mutate(net_benefit=1000*net_benefit) %>% 
  pivot_wider(
    names_from=type,values_from=net_benefit
  ) %>% 
  group_by(
    threshold
  ) %>% 
  fill(x1,.direction = "updown") %>% 
  fill(x2,.direction = "updown") %>% 
  fill(x3,.direction = "updown") %>% 
  select(threshold,x1,x2,x3) %>% 
  unique() %>% 
  mutate(
    delta.x2 = (x2-x1),
    delta.x3 = (x3-x1),
    boot = myboot
  )

  df.relative <- rbind(df.relative,relative)

}

########
# FIGURE
########

# AUROC - could you group_by, but for following the derivations more easily without group_by

summary.auroc <- rbind(
  data.frame(
    type       = "x1",
    auroc.mean = filter(df.auroc,type=="x1") %>% na.omit() %>% pull(auc) %>% mean(),
    auroc.se   = sd(filter(df.auroc,type=="x1") %>% na.omit() %>% pull(auc))/sqrt(n.boot)
  ),
  data.frame(
    type       = "x1+x2",
    auroc.mean = filter(df.auroc,type=="x1+x2") %>% na.omit() %>% pull(auc) %>% mean(),
    auroc.se   = sd(filter(df.auroc,type=="x1+x2") %>% na.omit() %>% pull(auc))/sqrt(n.boot)
  ),
  data.frame(
    type       = "x1+x2+I(x2)",
    auroc.mean = filter(df.auroc,type=="x1+x2+I(x2)") %>% na.omit() %>% pull(auc) %>% mean(),
    auroc.se   = sd(filter(df.auroc,type=="x1+x2+I(x2)") %>% na.omit() %>% pull(auc))/sqrt(n.boot)
  )) %>% 
  mutate(
    auroc.lower = auroc.mean-1.96*auroc.se,
    auroc.upper = auroc.mean+1.96*auroc.se
  ) %>% 
  mutate(mystr = paste0(
    format(round(auroc.mean,2),nsmall = 2)))

# DCA

df.dca.scaled = df.dca %>% mutate(net_benefit = 1000*net_benefit)
summary.dca <- data.frame()

for (mythres in unique(df.dca.scaled$threshold)){
  summary.dca <- rbind(summary.dca,rbind(
  data.frame(
    threshold  = mythres,
    type       = "x1",
    mymean     = filter(df.dca.scaled,type=="x1" & threshold==mythres) %>% na.omit() %>% pull(net_benefit) %>% mean(),
    dca.se     = sd(filter(df.dca.scaled,type=="x1" & threshold==mythres) %>% na.omit() %>% pull(net_benefit))/sqrt(n.boot)
  ),
  data.frame(
    threshold  = mythres,
    type       = "x1+x2",
    mymean     = filter(df.dca.scaled,type=="x1+x2" & threshold==mythres) %>% na.omit() %>% pull(net_benefit) %>% mean(),
    dca.se     = sd(filter(df.dca.scaled,type=="x1+x2" & threshold==mythres) %>% na.omit() %>% pull(net_benefit))/sqrt(n.boot)
  ),
  data.frame(
    threshold  = mythres,
    type       = "x1+x2+I(x2)",
    mymean     = filter(df.dca.scaled,type=="x1+x2+I(x2)" & threshold==mythres) %>% na.omit() %>% pull(net_benefit) %>% mean(),
    dca.se     = sd(filter(df.dca.scaled,type=="x1+x2+I(x2)" & threshold==mythres) %>% na.omit() %>% pull(net_benefit))/sqrt(n.boot)
  )) %>% 
  mutate(
    mylower = mymean-1.96*dca.se,
    myupper = mymean+1.96*dca.se
  ) %>% select(-dca.se)) 
}

fig.dca.model <- summary.dca %>% 
  mutate(type = case_when(
    type=="x1"~paste0("Patient and surgery characteristics\nAUROC: ",summary.auroc$mystr[1]),
    type=="x1+x2"~paste0("Patient and surgery characteristics + eGFR\nAUROC: ",summary.auroc$mystr[2]),
    type=="x1+x2+I(x2)"~paste0("Patient and surgery characteristics + eGFR + eGFR^2\nAUROC: ",summary.auroc$mystr[3]))) %>% 
  select(threshold,mymean,mylower,myupper,type)
fig.dca.treat <- df.dca.treat %>% mutate(net_benefit=1000*net_benefit) %>% 
  group_by(threshold,label) %>% 
  summarise(
    mymean  = mean(net_benefit,na.rm=T),
    mylower = mean(net_benefit,na.rm=T),
    myupper = mean(net_benefit,na.rm=T)
  ) %>% 
  mutate(type = label) %>% 
  select(-label) %>% 
  select(threshold,mymean,mylower,myupper,type)

figA <- rbind(fig.dca.model,fig.dca.treat) %>% 
  mutate(
    type = factor(type,levels = c("Treat All","Treat None",
                                  paste0("Patient and surgery characteristics\nAUROC: ",summary.auroc$mystr[1]),
                                  paste0("Patient and surgery characteristics + eGFR\nAUROC: ",summary.auroc$mystr[2]),
                                  paste0("Patient and surgery characteristics + eGFR + eGFR^2\nAUROC: ",summary.auroc$mystr[3])))) %>% 
  ggplot(aes(x=threshold,y=mymean,ymin=mylower,ymax=myupper,color=type,fill=type))+
  geom_ribbon(alpha=0.1,colour=NA)+
  geom_line(size=1)+
  scale_x_continuous(labels=scales::percent,limits = c(0,0.6),n.breaks = 10)+
  scale_y_continuous(n.breaks = 10,breaks = seq(0,150,by=20),limits = c(-12,150))+
  xlab("Threshold probability")+
  theme_classic()+
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size=14)
  )+
  ggsci::scale_color_jama()+
  ggsci::scale_fill_jama()+
  ylab("Net benefit (per 1000 patients)")+
  guides(color = guide_legend(nrow = 2))+
  labs(title="A")

#######
# inlet
#######

mycol = ggsci::pal_jama()(5)[c(4,5)]

figB <- df.relative %>% 
  pivot_longer(contains("delta")) %>% 
  group_by(
    threshold,name
  ) %>% 
  summarise(
    mymean = mean(value,na.rm=T),
    mylower = mean(value,na.rm=T)-sd(value,na.rm=T),
    myupper = mean(value,na.rm=T)+sd(value,na.rm=T)
  ) %>% 
  ggplot(aes(x=threshold,y=mymean,ymin=mylower,ymax=myupper,color=name,fill=name))+
  geom_ribbon(alpha=0.1,colour=NA)+
  geom_line(size=1.2)+
  scale_x_continuous(labels=scales::percent,limits = c(0,0.6),n.breaks = 10)+
  scale_y_continuous(n.breaks = 10)+
  xlab("Threshold probability")+
  theme_classic()+
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size=13)
  )+
  scale_color_manual(values=mycol)+
  scale_fill_manual(values=mycol)+
  ylab("Increase in net benefit\n(per 1000 patients)")+
  labs(title="B")

# plot.with.inset <-
ggdraw() +
  draw_plot(figA) +
  draw_plot(figB, x = 0.4, y = .54, width = .53, height = .37)

ggsave("C:/My Programs KAS/KAS Research/M. Huber/AUROC vs DCA/BJA comment/revisions2/figure_revisions2_02mar25.jpeg",dpi=1000,width=10,height=6.5)


#########################################
# values for paper
# use monte carlo s.e. and derived 95%-CI
#########################################

df.dca %>%
  mutate(net_benefit=1000*net_benefit) %>% 
  group_by(threshold,type) %>% 
  summarise(
    mymean  = mean(net_benefit,na.rm=T),
    myse    = sd(net_benefit)/sqrt(n.boot)) %>% 
  mutate(
    mylower = mymean-1.96*myse,
    myupper = mymean+1.96*myse
  ) %>% 
  filter(type=="x1" & threshold==0.2)

df.relative %>% 
  pivot_longer(contains("delta")) %>% 
  group_by(
    threshold,name
  ) %>% 
  summarise(
    mymean = mean(value,na.rm=T),
    myse   = sd(value)/sqrt(n.boot)) %>% 
  mutate(
    mylower = mymean-1.96*myse,
    myupper = mymean+1.96*myse
  ) %>% 
  filter(name=="delta.x2" & threshold==0.2)


df.relative %>% 
  pivot_longer(contains("delta")) %>% 
  group_by(
    threshold,name
  ) %>% 
  summarise(
    mymean = mean(value,na.rm=T),
    myse   = sd(value)/sqrt(n.boot)) %>% 
  mutate(
    mylower = mymean-1.96*myse,
    myupper = mymean+1.96*myse
  ) %>% 
  filter(name=="delta.x3" & threshold==0.2)
