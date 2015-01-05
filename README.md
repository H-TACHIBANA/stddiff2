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

fit <- glm(treat ~ age + crp + severe + score, family = binomial(logit))
summary(fit)

ps <- predict(fit, type = "response")
invwt <- treat/ps + (1-treat)/(1-ps)
wt <- invwt * treat * sum(treat) / sum(invwt * treat) + (1 - treat) * invwt * sum(treat == 0) / sum(invwt * (1 - treat))

data <- data.frame(wt = wt, treat = treat, age = age, crp = crp, severe = severe, score = score)

res <- stddiff.sas(dat = data, 
                   wt.var = "wt",       
                   treat.var = "treat",  #treat var must be 0/1
                   val.cont = c("age", "crp"),  #continuous var
                   val.cat = c("severe", "score"), #categolical var
                   val.ord = c("score")  #ranked var 
)

tbl.1 <- stdiff.table.cont(res)
tbl.2 <- stdiff.table.cat(res)
tbl.3 <- stdiff.table.rank(res)

tbl.1
tbl.2
tbl.3

##result

# > tbl.1
# $table.cont.uwt
# variable treat.mean treat.sd control.mean control.sd p.value stddiff
# 1      age       49.8    10.01         50.1       9.96   0.636   -0.03
# 2      crp      12.81    5.865        10.17      4.899   <.001   0.489
# 
# $table.cont.wt
# variable treat.mean treat.sd control.mean control.sd p.value stddiff
# 1      age       50.2    10.02           50      9.921    0.77    0.02
# 2      crp      11.24    5.957        11.14      4.909   0.794   0.018
# 
# > tbl.2
# $table.cat.uwt
# variable categoly number.treat1 percent.treat1 number.treat0 percent.treat0 p.value stddiff
# 1    severe                                                                      <.001   0.289
# 2                  0           265             53           335             67                
# 3                  1           235             47           165             33                
# 4     score                                                                      <.001   0.397
# 5                  0            83           16.6            45              9                
# 6                  1            99           19.8            55             11                
# 7                  2            73           14.6            67           13.4                
# 8                  3            90             18           123           24.6                
# 9                  4            81           16.2           103           20.6                
# 10                 5            74           14.8           107           21.4                
# 
# $table.cat.wt
# variable categoly number.treat1 percent.treat1 number.treat0 percent.treat0 p.value stddiff
# 1    severe                                                                      0.872   0.011
# 2                  0         304.5           60.9         307.2           61.4                
# 3                  1         195.5           39.1         192.8           38.6                
# 4     score                                                                      0.315   0.165
# 5                  0          62.8           12.6          60.5           12.1                
# 6                  1          86.1           17.2          71.2           14.2                
# 7                  2          70.5           14.1          73.8           14.8                
# 8                  3            88           17.6         117.9           23.6                
# 9                  4          93.7           18.7          90.3           18.1                
# 10                 5          98.9           19.8          86.4           17.3                
# 
# > tbl.3
# $table.rank.uwt
# variable treat.mean treat.sd control.mean control.sd p.value stddiff
# 1    score      450.6    288.1        550.4        272   <.001  -0.356
# 
# $table.rank.wt
# variable treat.mean treat.sd control.mean control.sd p.value stddiff
# 1    score      502.6    291.7        502.4      277.9   0.927   0.001


```
