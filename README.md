# parametric-pooled-data
---
title: "pooled_data_parametric"
output: html_document
author: Thao Le
date: 21 March 2016
---

Fitting parametric models for pooled TBM data for HIV uninfected patients based on baseline data.
```{r global_options, include=FALSE}
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE)
```

```{r load data}
rm(list = ls())
setwd("W:/Project/TBMpooled/Thao/Analysis")
set.seed(123)

load(file = "W:/Project/TBMpooled/Thao/Analysis/Data/dt.Rdata")
load(file = "W:/Project/TBMpooled/Thao/Analysis/Data/pool_imp.Rdata")
load(file = "W:/Project/TBMpooled/Thao/Analysis/Data/pool_imp_long.Rdata")

source("W:/Project/TBMpooled/Thao/Summary/R functions/miscFun.R")
source("W:/Project/TBMpooled/Thao/Analysis/4_R_Function_Par.R")
library(survival)
library(Hmisc)
library(knitr)
library(mgcv)
library(plyr)
library(splines)
library(mice)
library(rms)
library(glmnet)
library(c060)
library(peperr)
library(pec)
library(timereg)
library(timeROC)
library(party)

dt$CSF_Lymp_Count <- dt$CSF_WCC*dt$CSF_Lymp/100

data.imp.long$logHist <- log2(data.imp.long$HistDay)
data.imp.long$logWCC  <- log2(data.imp.long$CSF_WCC)
data.imp.long$logPro  <- log2(data.imp.long$CSF_Pro)
data.imp.long$logLymp <- log2(data.imp.long$CSF_Lymp_Count)
summary(data.imp.long$logLymp)

d1 <- subset(data.imp.long, .imp == 1)
dt$logHist <- log2(dt$HistDay)
dt$logWCC <- log2(dt$CSF_WCC)
dt$logPro <- log2(dt$CSF_Pro)
dt$logLymp <- log2(dt$CSF_Lymp_Count)

fml.full <- Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Sex + Weight +
  Dex + TBMGrade + Prev_TB + Temp + FocalSign + Convulsion + CNP + Hemi + Para +
  Sodium + Glu + CSF_BGlu_Ratio + CSF_Lact + Resistance + CXR_Miliary + logHist + logWCC+
  logPro + logLymp
```


# Accelerated Failure Time models
Distributions for survival time being under consideration are: Weibull, lognormal and loglogistic.
Automatic model selection methods: Backward selection, Lasso, Stab
Explore functional form of covariates: assume all are linear
"The key assumption for an AFT model is that survival time accelerates (or decelerates)
by a constant factor when comparing different levels of covariates"
Should check whether parametric models fit the data well by exploring qq plot or residual plot
## Base on complete case analysis: pre-defined model

```{r Univariate}
Sdt <- with(dt, Surv(pmin(Ttdead,271), ifelse(Ttdead<= 271, Evdead, 0)))
sapply(dt[,c(4:27,31:36)], function(x){
   t <- summary(survreg(Sdt ~ x, data = dt, dist = "wei"))
   round(1- pchisq(t$chi, (t$df - t$idf)),4)})
sapply(dt[,c(4:27,31:36)], function(x){
   t <- summary(survreg(Sdt ~ x, data = dt, dist = "exp"))
   round(1- pchisq(t$chi, (t$df - t$idf)),4)})
sapply(dt[,c(4:27,31:36)], function(x){
   t <- summary(survreg(Sdt ~ x, data = dt, dist = "gau"))
   round(1- pchisq(t$chi, (t$df - t$idf)),4)})
sapply(dt[,c(4:27,31:36)], function(x){
   t <- summary(survreg(Sdt ~ x, data = dt, dist = "logistic"))
   round(1- pchisq(t$chi, (t$df - t$idf)),4)})
sapply(dt[,c(4:27,31:36)], function(x){
   t <- summary(survreg(Sdt ~ x, data = dt, dist = "lognormal"))
   round(1- pchisq(t$chi, (t$df - t$idf)),4)})
```

All AFT models select the same variables in univariate analysis: Age, Weight, Dex, (Hospital), TBMGrade, Prev_TB, GCS, FocalSign, Hemi, Resitance, logHist, logWCC, logLymp.
Fit models with predefined covariates based on completed case: 
```{r}
(wm <- survreg(Sdt ~ Age + Weight + Dex + TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
                Resistance + logWCC + logLymp + logHist + CSF_Lact,data = dt, dist = "w" ))
(em <- survreg(Sdt ~ Age + Weight + Dex + TBMGrade + Prev_TB + FocalSign + Sodium  + CSF_BGlu_Ratio + 
                Resistance + logWCC + logLymp + logHist + CSF_Lact,data = dt, dist = "exp" ))
(gaum <- survreg(Sdt ~ Age + Weight + Dex + TBMGrade + Prev_TB + FocalSign + Sodium  + CSF_BGlu_Ratio + 
                Resistance + logWCC + logLymp + logHist + CSF_Lact,data = dt, dist = "gau" ))
(logisticm <- survreg(Sdt ~ Age + Weight + Dex + TBMGrade + Prev_TB + FocalSign + Sodium  + CSF_BGlu_Ratio + 
                Resistance + logWCC + logLymp + logHist + CSF_Lact,data = dt, dist = "logistic" ))
(lognm <- survreg(Sdt ~ Age + Weight + Dex + TBMGrade + Prev_TB + FocalSign + Sodium  + CSF_BGlu_Ratio + 
                Resistance + logWCC + logLymp + logHist + CSF_Lact,data = dt, dist = "lognormal" ))

c(AIC(wm),AIC(em),AIC(gaum),AIC(logisticm),AIC(lognm))
```

Weibull and lognormal models gave smallest AIC values. Those models only use 399 patients, 552 patients deleted due to missingness.

## Base on imputed dataset - pre-defined model
Check if higher degree of continuous variable can improve the fit

```{r}
form.linear <- Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0)) ~  Age + 
  Weight + Dex + TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + Resistance + 
  I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + CSF_Lact

# df = 3
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 3, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "wei")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 3, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "exp")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 3, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "gau")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 3, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "logistic")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 3, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "lognormal")

# df = 4
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 4, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "wei")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 4, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "exp")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 4, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "gau")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 4, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "logistic")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 4, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "lognormal")

# df = 5
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 5, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "wei")
impAFTlinear(vars = labels(terms(form.linear)),
               outcome = "Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
               Data = data.imp, df = 5, 
               testvars = c("I(log2(CSF_WCC))", "I(log2(HistDay))", "CSF_Lact", "I(log2(CSF_Lymp_Count))", "Age", "Weight", "Sodium"), dist = "lognormal")

```

Check if there are any interactions between Age and Sex, Age and others, Sex and others
```{r check interaction}
# Compare with model with interaction between Age*sex

## Age * Sex
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Weight + Dex + 
    TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact + Sex + Age * Sex, dist = "wei")), fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Weight + 
    Dex + TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact, dist = "w")))$pvalue # not significant
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Weight + Dex + 
    TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact + Sex + Age * Sex, dist = "lognormal")), fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Weight + 
    Dex + TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact, dist = "lognormal")))$pvalue # not significant


# Age*[other]
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ Age *( Weight + Dex + 
    TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact), dist = "wei")), fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Weight + 
    Dex + TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact, dist = "wei")))$pvalue
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ Age *( Weight + Dex + 
    TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact), dist = "lognormal")), fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Weight + 
    Dex + TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact, dist = "lognormal")))$pvalue

# Sex*[other]
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ (Age + Weight + Dex + 
    TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact) * Sex, dist = "wei")), fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Weight + 
    Dex + TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact, dist = "wei")))$pvalue
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ (Age + Weight + Dex + 
    TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact) * Sex, dist = "lognormal")), fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Weight + 
    Dex + TBMGrade + Prev_TB + FocalSign + Sodium + CSF_BGlu_Ratio + 
    Resistance + I(log2(CSF_WCC)) + I(log2(HistDay)) + I(log2(CSF_Lymp_Count)) + 
    CSF_Lact, dist = "lognormal")))$pvalue

```

Final pre-defined models with all linear term
```{r}
# Question: how to pool scale ??
fit.linear.w <- with(data.imp, eval(parse(
  text = paste("expression(survreg(",deparse(form.linear, width.cutoff = 300L), "))"))))
t <- fit.linear.w$analyses[[1]]
sapply(fit.linear.w$analyses, function(x){x$scale})

```

## Based on imputed data - Automate selection 
### Stepwise slection

```{r backward selection}
# Weibull
bw.wei <- by(data = data.imp.long, 
             INDICES = data.imp.long$.imp,
             FUN = function(x){tmp <- selectAFT(formula = fml.full, data = x, rule = "aic", dist = "weibull")$In})
table(unlist(bw.wei))
fml.bw.wei <- Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Dex + Hemi + 
  I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact
impAFTDrop1(vars = labels(terms(fml.bw.wei)),
            outcome ="Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
            Data = data.imp, dist = "weibull")
## Interaction with Sex
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ (Age + Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact) * Sex, dist = "wei")), 
    fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact, dist = "wei")))$pvalue
## Interaction with Age
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ Age*(Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact), dist = "wei")), 
    fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact, dist = "wei")))$pvalue

# lognormal
bw.lognormal <- by(data = data.imp.long, 
             INDICES = data.imp.long$.imp,
             FUN = function(x){tmp <- selectAFT(formula = fml.full, data = x, rule = "aic", dist = "lognormal")$In})
table(unlist(bw.lognormal))
fml.bw.lognormal <- Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Dex + Hemi + 
  I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact
impAFTDrop1(vars = labels(terms(fml.bw.lognormal)),
            outcome ="Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0))",
            Data = data.imp, dist = "lognormal")
## Interaction with Sex
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ (Age + Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact) * Sex, dist = "lognormal")), 
    fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact, dist = "lognormal")))$pvalue
## Interaction with Age
pool.compare(fit1 = with(data.imp, survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ Age*(Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact), dist = "lognormal")), 
    fit0 = with(data.imp, survreg(Surv(pmin(Ttdead, 
    271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact, dist = "lognormal")))$pvalue

AIC(survreg(Surv(pmin(Ttdead, 271), 
    ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact, dist = "wei", data = d1))
AIC(survreg(Surv(pmin(Ttdead, 271), ifelse(Ttdead <= 271, Evdead, 0)) ~ Age + Dex + Hemi + 
    I(log2(CSF_Lymp_Count)) + Resistance + TBMGrade + CSF_Lact, dist = "lognormal", data = d1))


```


 
```{r}
# w.model <- survreg(Surv(pmin(Ttdead,271), ifelse(Ttdead<= 271, Evdead, 0)) ~ 
#                      Age + TBMGrade + Dex, data = d1, dist = "w")
# w.model
# y <- survfit(Surv(pmin(Ttdead,271), ifelse(Ttdead<= 271, Evdead, 0)) ~ TBMGrade, d1)
# plot(y)
# plot(y, fun = "cloglog",
#      xlab = "time in days using logarithmic scale",
#      ylab = "log-log survival", 
#      main = "log-log curves")
# yw <- survreg(Surv(pmin(Ttdead,271), ifelse(Ttdead<= 271, Evdead, 0)) ~ 
#                      TBMGrade , data = d1, dist = "w")
# 
# library(flexsurv)
# S <- with(d1, Surv(pmin(Ttdead,271), ifelse(Ttdead<= 271, Evdead, 0)))
# 
# sWei  <- flexsurvreg(S ~ TBMGrade, dist = 'weibull', data = d1)
# sLno  <- flexsurvreg(S ~ TBMGrade, dist = 'lnorm', data = d1)   
# plot(sWei)
# lines(sLno, col="blue")
# 
# sWei  <- flexsurvreg(S ~ 1, dist = 'weibull', data = d1)
# sLno  <- flexsurvreg(S ~ 1, dist = 'lnorm', data = d1)   
# plot(sWei, ci = F)
# lines(sLno, col="blue")























```



