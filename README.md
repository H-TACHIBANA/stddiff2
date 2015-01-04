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
```

