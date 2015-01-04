stddiff.sas <- function(   dat = data,        #data.frame
                           wt.var = "w",      #weight(already calculated to ATE or ATT or others)
                           treat.var = "treat",    #treat must be 1 / 0　
                           val.cont = "cont",      #continuous var
                           val.cat = "cat", #categolical var
                           val.ord = "O"     #ranked var
){
      library(survey)
      wt.mean <- function (x, wt){  #this is copy of package "SDMTools"
            s = which(is.finite(x * wt))
            wt = wt[s]
            x = x[s]
            return(sum(wt * x)/sum(wt))
      }
      wt.sd.sas  <- function (x, wt){ #this is from SAS proc means with weight
            s = which(is.finite(x + wt))  #remove NA,INF
            wt = wt[s]
            x = x[s]
            xbar = wt.mean(x, wt)
            return(sqrt(1 / (length(x) - 1) * sum(wt * (x - xbar)^2)))
      }
      stddiff.cat <- function(p1, p0){ #standardized difference for categorical variables
            S2 <- matrix(0, nrow = (length(p1) - 1), ncol = (length(p1) - 1))
            
            for (i in seq_len(length(p0) - 1)){
                  for (j in seq_len(length(p0) - 1)){
                        if (i == j){
                              S2[i,j] <- (p1[i]*(1-p1[i]) + p0[i] * (1-p0[i])) / 2
                        }else{
                              S2[i,j] <- - (p1[i]*p1[j] + p0[i]*p0[j]) / 2  #this line is different from original paper, but same as SAS proc %stddiff
                        }}
            }
            
            sqrt(t(p1[-length(p1)] - p0[-length(p0)]) %*% solve(S2) %*% (p1[-length(p1)] - p0[-length(p0)]))
      }
      stddiff.cont <- function(mu1, sd1, mu0, sd0){
            (mu1 - mu0)/sqrt((sd1^2 + sd0^2)/2)
      }
      
      res <- list()          #result
      
      #stabilized PS, treat must be 1/0
      dat$wt <- dat[,wt.var] * dat[,treat.var] * sum(dat[,treat.var]) / sum(dat[,wt.var] * dat[,treat.var]) + (1 - dat[,treat.var]) * dat[,wt.var] * sum(dat[,treat.var] == 0) / sum(dat[,wt.var] * (1 - dat[,treat.var]))
      dat$treat <- dat[,treat.var]
      
      treat0 <- subset(dat, treat == 0)
      treat1 <- subset(dat, treat == 1)
      
      for (j in c(val.cont, val.ord, val.cat)){
            res[[j]] <- list()
      }
      
      for (tmp1 in val.cont){   #continuous
            sd0 <- wt.sd.sas(treat0[,tmp1], treat0$wt)
            sd1 <- wt.sd.sas(treat1[,tmp1], treat1$wt)
            mu0 <- wt.mean(treat0[,tmp1], treat0$wt)
            mu1 <- wt.mean(treat1[,tmp1], treat1$wt)
            
            res[[tmp1]]["mean.treat0.cont.wt"] <- mu0
            res[[tmp1]]["sd.treat0.cont.wt"] <- sd0
            res[[tmp1]]["mean.treat1.cont.wt"] <- mu1
            res[[tmp1]]["sd.treat1.cont.wt"] <- sd1
            res[[tmp1]]["stddiff.cont.wt"] <- stddiff.cont(mu1, sd1, mu0, sd0)
            
            dat$tmp6 <- dat[,tmp1]
            design.temp <- svydesign(ids=~1, weights=~dat$wt, data=dat)
            res[[tmp1]]["ttest.p.cont.wt"] <- svyttest(tmp6 ~ treat, design = design.temp)$p.value
            
            sd0u <- sd(treat0[,tmp1])
            sd1u <- sd(treat1[,tmp1])
            mu0u <- mean(treat0[,tmp1])
            mu1u <- mean(treat1[,tmp1])
            
            res[[tmp1]]["mean.treat0.cont.uwt"] <- mu0u
            res[[tmp1]]["sd.treat0.cont.uwt"] <- sd0u
            res[[tmp1]]["mean.treat1.cont.uwt"] <- mu1u
            res[[tmp1]]["sd.treat1.cont.uwt"] <- sd1u
            res[[tmp1]]["stddiff.cont.unw"] <- stddiff.cont(mu1u, sd1u, mu0u, sd0u)
            res[[tmp1]]["ttest.p.cont.uwt"] <- t.test(dat$tmp6 ~ dat$treat)$p.value
      }
      
      for (tmp2 in val.cat){  #categoly
            dat$tmp3 <- dat[,tmp2]
            design.temp <- svydesign(ids=~1, weights=~dat$wt, data=dat)
            res[[tmp2]][["chisq.p.value.wt"]] <- svychisq(~tmp3 + treat, design = design.temp)$p.value
            res[[tmp2]][["table.wt"]] <- tbl.wt <- svytable(~tmp3 + treat, design.temp)
            res[[tmp2]][["treat0.No.wt"]] <- tbl.wt[,1]  
            res[[tmp2]][["treat0.rate.wt"]] <- p0 <- tbl.wt[,1] / sum(tbl.wt[,1])
            res[[tmp2]][["treat1.No.wt"]] <- tbl.wt[,2]  
            res[[tmp2]][["treat1.rate.wt"]] <- p1 <- tbl.wt[,2] / sum(tbl.wt[,2])
            res[[tmp2]][["stddiff.cat.wt"]]  <- stddiff.cat(p1, p0)
            
            res[[tmp2]][["table.unw"]] <- tbl.uwt <- table(dat$tmp3 , dat$treat)
            res[[tmp2]][["chisq.p.value.unw"]] <- chisq.test(tbl.uwt)$p.value
            res[[tmp2]][["treat0.No.uwt"]] <- tbl.uwt[,1]  
            res[[tmp2]][["treat0.rate.uwt"]] <- p0u <- tbl.uwt[,1] / sum(tbl.uwt[,1])
            res[[tmp2]][["treat1.No.uwt"]] <- tbl.uwt[,2]  
            res[[tmp2]][["treat1.rate.uwt"]] <- p1u <- tbl.uwt[,2] / sum(tbl.uwt[,2])
            res[[tmp2]][["stddiff.cat.uwt"]]  <- stddiff.cat(p1u, p0u)
      }
      
      
      for (tmp4 in val.ord){   #rank/ordinal
            dat$temp5 <- rank(dat[,tmp4], ties.method = "average")
            treat0 <- subset(dat, treat == 0)
            treat1 <- subset(dat, treat == 1)
            
            sd0 <- wt.sd.sas(treat0$temp5, treat0$wt)
            sd1 <- wt.sd.sas(treat1$temp5, treat1$wt)
            mu0 <- wt.mean(treat0$temp5, treat0$wt)
            mu1 <- wt.mean(treat1$temp5, treat1$wt)
            
            res[[tmp4]]["mean.treat0.rank.wt"] <- mu0
            res[[tmp4]]["sd.treat0.rank.wt"] <- sd0
            res[[tmp4]]["mean.treat1.rank.wt"] <- mu1
            res[[tmp4]]["sd.treat1.rank.wt"] <- sd1
            res[[tmp4]]["stddiff.rank.wt"] <- stddiff.cont(mu1, sd1, mu0, sd0)
            
            dat$tmp6 <- dat[,tmp4]
            design.temp <- svydesign(ids=~1, weights=~dat$wt, data=dat)
            res[[tmp4]]["ttest.p.rank.wt"] <- svyttest(tmp6 ~ treat, design = design.temp)$p.value
            
            sd0u <- sd(treat0$temp5)
            sd1u <- sd(treat1$temp5)
            mu0u <- mean(treat0$temp5)
            mu1u <- mean(treat1$temp5)
            
            res[[tmp4]]["mean.treat0.rank.uwt"] <- mu0u
            res[[tmp4]]["sd.treat0.rank.uwt"] <- sd0u
            res[[tmp4]]["mean.treat1.rank.uwt"] <- mu1u
            res[[tmp4]]["sd.treat1.rank.uwt"] <- sd1u
            res[[tmp4]]["stddiff.rank.uwt"] <- stddiff.cont(mu1u, sd1u, mu0u, sd0u)
            res[[tmp4]]["ttest.p.rank.uwt"] <- t.test(dat$tmp6 ~ dat$treat)$p.value
            
      }
      
      res[["wt.var"]] = wt.var
      res[["treat.var"]] = treat.var 
      res[["val.cont"]] = val.cont
      res[["val.cat"]] = val.cat
      res[["val.ord"]] = val.ord
      return(res)
}







stdiff.table.cont <- function(res = res2){
      table1 <- list()
      val.cont <- res[["val.cont"]]
      val.cat <- res[["val.cat"]]
      val.ord <- res[["val.ord"]]
      
      table.uwt <- data.frame(NULL)
      for (i1 in val.cont){  
            table.uwt <- rbind(table.uwt, data.frame(variable = i1,
                                                     treat.mean = as.character(signif(as.numeric(res[[i1]]["mean.treat1.cont.uwt"]), digits = 4)),
                                                     treat.sd = as.character(signif(as.numeric(res[[i1]]["sd.treat1.cont.uwt"]), digits =4)),
                                                     control.mean = as.character(signif(as.numeric(res[[i1]]["mean.treat0.cont.uwt"]), digits=4)),
                                                     control.sd = as.character(signif(as.numeric(res[[i1]]["sd.treat0.cont.uwt"]), digits=4)),
                                                     p.value = as.character(round(as.numeric(res[[i1]]["ttest.p.cont.uwt"]), digits = 3)),
                                                     stddiff = as.character(round(as.numeric(res[[i1]]["stddiff.cont.unw"]), digits = 3))
            ))
      }
      rownames(table.uwt) <- NULL
      for (i2 in colnames(table.uwt)){
            table.uwt[,i2] <- as.character(table.uwt[,i2]) 
      }
      table.uwt$p.value[table.uwt$p.value == "0"] <- "<.001"
      table.uwt$stddiff[table.uwt$stddiff == "0"] <- "<.001"
      
      table.wt <- data.frame(NULL)
      for (i3 in val.cont){  
            table.wt <- rbind(table.wt, data.frame(variable = i3,
                                                   treat.mean = as.character(signif(as.numeric(res[[i3]][["mean.treat1.cont.wt"]]), digits = 4)),
                                                   treat.sd = as.character(signif(as.numeric(res[[i3]]["sd.treat1.cont.wt"]), digits =4)),
                                                   control.mean = as.character(signif(as.numeric(res[[i3]]["mean.treat0.cont.wt"]), digits=4)),
                                                   control.sd = as.character(signif(as.numeric(res[[i3]]["sd.treat0.cont.wt"]), digits=4)),
                                                   p.value = as.character(round(as.numeric(res[[i3]]["ttest.p.cont.wt"]), digits = 3)),
                                                   stddiff = as.character(round(as.numeric(res[[i3]]["stddiff.cont.wt"]), digits = 3))
            ))
      }
      rownames(table.wt) <- NULL
      for (i4 in colnames(table.wt)){
            table.wt[,i4] <- as.character(table.wt[,i4]) 
      }
      table.wt$p.value[table.wt$p.value == "0"] <- "<.001"
      table.wt$stddiff[table.wt$stddiff == "0"] <- "<.001"
      
      table1[["table.cont.uwt"]] <- table.uwt
      table1[["table.cont.wt"]] <- table.wt
      return(table1)
}





stdiff.table.cat <- function(res = res2){
      table1 <- list()
      val.cont <- res[["val.cont"]]
      val.cat <- res[["val.cat"]]
      val.ord <- res[["val.ord"]]
      
      table.cat.uwt <- data.frame(NULL)
      for (i5 in val.cat){  
            table.cat.uwt <- rbind(table.cat.uwt, data.frame(variable = i5, categoly = "",
                                                             number.treat1 = "", percent.treat1 = "",
                                                             number.treat0 = "", percent.treat0 = "",
                                                             p.value = as.character(round(res[[i5]][["chisq.p.value.unw"]], digits = 3)),
                                                             stddiff = as.character(round(as.numeric(res[[i5]]["stddiff.cat.uwt"]), digits = 3))))
            table.cat.uwt <- rbind(table.cat.uwt, 
                                   data.frame(variable = "",
                                              categoly = rownames(res[[i5]][["table.unw"]]),
                                              number.treat1 = as.character(round(res[[i5]][["treat1.No.uwt"]], digits = 1)),
                                              percent.treat1 = as.character(round(res[[i5]][["treat1.rate.uwt"]]*100, digits =1)),
                                              number.treat0 = as.character(round(res[[i5]][["treat0.No.uwt"]], digits = 1)),
                                              percent.treat0 = as.character(round(res[[i5]][["treat0.rate.uwt"]]*100, digits = 1)),
                                              p.value = "", stddiff = ""))
      }
      rownames(table.cat.uwt) <- NULL
      for (i6 in colnames(table.cat.uwt)){
            table.cat.uwt[,i6] <- as.character(table.cat.uwt[,i6]) 
      }
      
      table.cat.uwt$p.value[table.cat.uwt$p.value == "0"] <- "<.001"
      table.cat.uwt$stddiff[table.cat.uwt$stddiff == "0"] <- "<.001"
      
      table.cat.wt <- data.frame(NULL)
      for (i7 in val.cat){  
            table.cat.wt <- rbind(table.cat.wt, data.frame(variable = i7, categoly = "",
                                                           number.treat1 = "", percent.treat1 = "",
                                                           number.treat0 = "", percent.treat0 = "",
                                                           p.value = as.character(round(res[[i7]][["chisq.p.value.wt"]], digits = 3)),
                                                           stddiff = as.character(round(as.numeric(res[[i7]]["stddiff.cat.wt"]), digits = 3))))
            table.cat.wt <- rbind(table.cat.wt, 
                                  data.frame(variable = "",
                                             categoly = rownames(res[[i7]][["table.wt"]]),
                                             number.treat1 = as.character(round(res[[i7]][["treat1.No.wt"]], digits = 1)),
                                             percent.treat1 = as.character(round(res[[i7]][["treat1.rate.wt"]]*100, digits =1)),
                                             number.treat0 = as.character(round(res[[i7]][["treat0.No.wt"]], digits = 1)),
                                             percent.treat0 = as.character(round(res[[i7]][["treat0.rate.wt"]]*100, digits = 1)),
                                             p.value = "", stddiff = ""))
      }
      rownames(table.cat.wt) <- NULL
      for (i8 in colnames(table.cat.wt)){
            table.cat.wt[,i8] <- as.character(table.cat.wt[,i8]) 
      }
      
      table.cat.wt$p.value[table.cat.wt$p.value == "0"] <- "<.001"
      table.cat.wt$stddiff[table.cat.wt$stddiff == "0"] <- "<.001"
      
      table1[["table.cat.uwt"]] <- table.cat.uwt
      table1[["table.cat.wt"]] <- table.cat.wt
      return(table1)
}



stdiff.table.rank <- function(res = res2){
      table1 <- list()
      val.cont <- res[["val.cont"]]
      val.cat <- res[["val.cat"]]
      val.ord <- res[["val.ord"]]
      
      table.rank.uwt <- data.frame(NULL)
      
      stdiff.table.rank <- function(res = res2){
            table1 <- list()
            val.cont <- res[["val.cont"]]
            val.cat <- res[["val.cat"]]
            val.ord <- res[["val.ord"]]
            
            table.uwt <- data.frame(NULL)
            for (i1 in val.ord){  
                  table.uwt <- rbind(table.uwt, data.frame(variable = i1,
                                                           treat.mean = as.character(signif(as.numeric(res[[i1]]["mean.treat1.rank.uwt"]), digits = 4)),
                                                           treat.sd = as.character(signif(as.numeric(res[[i1]]["sd.treat1.rank.uwt"]), digits =4)),
                                                           control.mean = as.character(signif(as.numeric(res[[i1]]["mean.treat0.rank.uwt"]), digits=4)),
                                                           control.sd = as.character(signif(as.numeric(res[[i1]]["sd.treat0.rank.uwt"]), digits=4)),
                                                           p.value = as.character(round(as.numeric(res[[i1]]["ttest.p.rank.uwt"]), digits = 3)),
                                                           stddiff = as.character(round(as.numeric(res[[i1]]["stddiff.rank.uwt"]), digits = 3))
                  ))
            }
            rownames(table.uwt) <- NULL
            for (i2 in colnames(table.uwt)){
                  table.uwt[,i2] <- as.character(table.uwt[,i2]) 
            }
            table.uwt$p.value[table.uwt$p.value == "0"] <- "<.001"
            table.uwt$stddiff[table.uwt$stddiff == "0"] <- "<.001"
            
            table.wt <- data.frame(NULL)
            for (i3 in val.ord){  
                  table.wt <- rbind(table.wt, data.frame(variable = i3,
                                                         treat.mean = as.character(signif(as.numeric(res[[i3]][["mean.treat1.rank.wt"]]), digits = 4)),
                                                         treat.sd = as.character(signif(as.numeric(res[[i3]]["sd.treat1.rank.wt"]), digits =4)),
                                                         control.mean = as.character(signif(as.numeric(res[[i3]]["mean.treat0.rank.wt"]), digits=4)),
                                                         control.sd = as.character(signif(as.numeric(res[[i3]]["sd.treat0.rank.wt"]), digits=4)),
                                                         p.value = as.character(round(as.numeric(res[[i3]]["ttest.p.rank.wt"]), digits = 3)),
                                                         stddiff = as.character(round(as.numeric(res[[i3]]["stddiff.rank.wt"]), digits = 3))
                  ))
            }
            rownames(table.wt) <- NULL
            for (i4 in colnames(table.wt)){
                  table.wt[,i4] <- as.character(table.wt[,i4]) 
            }
            table.wt$p.value[table.wt$p.value == "0"] <- "<.001"
            table.wt$stddiff[table.wt$stddiff == "0"] <- "<.001"
            
            table1[["table.rank.uwt"]] <- table.uwt
            table1[["table.rank.wt"]] <- table.wt
      }
      return(table1)
}
