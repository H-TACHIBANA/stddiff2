---
title: "stddiff"
author: "H.Tachibana"
date: "Sunday, January 04, 2015"
output: html_document
---
Calculate the standardized difference and make the table one for propensity score analysis (simulate SAS %stddiff)

```{r}
library(stddiff2)

treat <- rep(c(0,1),500)
age <- rnorm(1000,50,10)
crp <- rnorm(1000,10,5) + treat * rnorm(1000, 3, 3)
severe <- sample(c(0,1), 1000, replace = TRUE) * treat + sample(c(0,0,1), 1000, replace = TRUE) * (1 - treat)
score <- sample(c(0,1,2,3,4,5), 1000, replace = TRUE) * treat + sample(c(0,1,2,3,3,4,4,5,5), 1000, replace = TRUE) * (1 - treat)

fit <- glm(treat ~ age + crp + severe + score)
summary(fit)

ps <- predict(fit, type = "response")
invwt <- treat/ps + (1-treat)/(1-ps)
wt <- invwt * treat * sum(treat) / sum(invwt * treat) + (1 - treat) * invwt * sum(treat == 0) / sum(invwt * (1 - treat))

data <- data.frame(wt = wt, treat = treat, age = age, crp = crp, severe = severe, score = score)

res <- stddiff.sas(dat = data, 
                   wt.var = "wt",       
                   treat.var = "treat",  #treat var must be 0/1
                   val.cont = c("age", "crp"),  #continuous var
                   val.cat = c("severe"), #categolical var
                   val.ord = c("score")  #ranked var
)

tbl.1 <- stdiff.table.cont(res)
tbl.2 <- stdiff.table.cat(res)
tbl.3 <- stdiff.table.rank(res)

tbl.1
tbl.2
tbl.3

# > tbl.1
# $table.cont.uwt
# variable treat.mean treat.sd control.mean control.sd p.value stddiff
# 1      age      50.02    9.756        50.09      9.902   0.913  -0.007
# 2      crp      13.45    5.453        10.09      4.995   <.001   0.643
# 
# $table.cont.wt
# variable treat.mean treat.sd control.mean control.sd p.value stddiff
# 1      age      50.12    9.716        50.09      9.855   0.966   0.003
# 2      crp      11.72    5.536        11.77      5.493   0.917  -0.009
# 
# > tbl.2
# $table.cat.uwt
# variable categoly number.treat1 percent.treat1 number.treat0 percent.treat0 p.value stddiff
# 1   severe                                                                      <.001   0.299
# 2                 0           252           50.4           325             65                
# 3                 1           248           49.6           175             35                
# 
# $table.cat.wt
# variable categoly number.treat1 percent.treat1 number.treat0 percent.treat0 p.value stddiff
# 1   severe                                                                      0.746   0.023
# 2                 0         289.8             58           284           56.8                
# 3                 1         210.2             42           216           43.2                
# 
# > tbl.3
# $table.rank.uwt
# variable treat.mean treat.sd control.mean control.sd p.value stddiff
# 1    score      446.2    287.6        554.8      271.1   <.001  -0.389
# 
# $table.rank.wt
# variable treat.mean treat.sd control.mean control.sd p.value stddiff
# 1    score      501.9    287.9        505.8      275.2   0.797  -0.014


```
