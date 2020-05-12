if(!is.bull(dev.list())) dev.off()
cat("\014") 
rm(list=ls())


library(tidyverse)
set.seed(123)

# prevs <- (read.csv("~/Desktop/pneumo/prevalences y0to6 final NEGwVTs jewish cat time 3.csv"))[,-1]
# prevs <- prevs[-55,]

ds <- (read.csv(file = "carriage_final_Jewish.csv"))[,2:9]
new.steff <- read.csv("new_st_eff.csv")
new2 <- as.data.frame(new.steff[,2])
names(new2) <- "Serotype"
df1 <- new2


foo <- read.csv("~/Desktop/pneumo/multinomialALLPCV13final.csv")[,3]
foo2 <- as.data.frame((as.character(foo)))
colnames(foo2) <- "st"

# prevRsamps <- read.csv("~/Desktop/pneumo/jsamps.jewish.prevR7.csv")
# prevsamps  <- read.csv("~/Desktop/pneumo/jsamps.jewish.prev7.csv")


#IPD.dist=rpois(IPDN)

#IPD <- (read.csv("~/Desktop/pneumo/ageIPD noVT plus .5 v3.csv"))
IPD <- (read.csv("isr_ipd_maile.csv"))
IPD_1 <- IPD[which(IPD$agec==1),] 
IPD_2 <- IPD[which(IPD$agec==2),] 
IPD_3 <- IPD[which(IPD$agec==3),] 
IPD_4 <- IPD[which(IPD$agec==4),] 
IPD_5 <- IPD[which(IPD$agec==5),] 

IPD_1m <- matrix(as.numeric(as.character(unlist(IPD_1[,3:9]))), nrow=dim(IPD_1[,3:9])[1], ncol=7, byrow = F)
IPD_2m <- matrix(as.numeric(as.character(unlist(IPD_2[,3:9]))), nrow=dim(IPD_2[,3:9])[1], ncol=7, byrow = F)
IPD_3m <- matrix(as.numeric(as.character(unlist(IPD_3[,3:9]))), nrow=dim(IPD_3[,3:9])[1], ncol=7, byrow = F)
IPD_4m <- matrix(as.numeric(as.character(unlist(IPD_4[,3:9]))), nrow=dim(IPD_4[,3:9])[1], ncol=7, byrow = F)
IPD_5m <- matrix(as.numeric(as.character(unlist(IPD_5[,3:9]))), nrow=dim(IPD_5[,3:9])[1], ncol=7, byrow = F)
rownames(IPD_1m) <- as.character((IPD_1$serotype))
rownames(IPD_2m) <- as.character((IPD_2$serotype))
rownames(IPD_3m) <- as.character((IPD_3$serotype))
rownames(IPD_4m) <- as.character((IPD_4$serotype))
rownames(IPD_5m) <- as.character((IPD_5$serotype))
colnames(IPD_1m) <- colnames(IPD_2m) <- colnames(IPD_3m) <- colnames(IPD_4m) <- colnames(IPD_5m) <- c("t0","t1","t2","t3","t4","t5","t6")

IPD_1m[9,] <- IPD_1m[9,]+IPD_1m[10,]+IPD_1m[11,]
rownames(IPD_1m)[9] <- "15BC"
IPD_1m <- IPD_1m[-c(10,11,46),]
rownames(IPD_1m)[dim(IPD_1m)[1]] <- "NEG"

rownames(IPD_2m)[8] <- "15BC"
IPD_2m <- IPD_2m[-(dim(IPD_2m)[1]-1),]
rownames(IPD_2m)[dim(IPD_2m)[1]] <- "NEG"

IPD_3m[9,] <- IPD_3m[9,]+IPD_3m[10,]
rownames(IPD_3m)[9] <- "15BC"
IPD_3m <- IPD_3m[-c(10,(dim(IPD_3m)[1]-1)),]
rownames(IPD_3m)[dim(IPD_3m)[1]] <- "NEG"

IPD_4m[9,] <- IPD_4m[9,]+IPD_4m[10,]
rownames(IPD_4m)[9] <- "15BC"
IPD_4m <- IPD_4m[-c(10,(dim(IPD_4m)[1]-1)),]
rownames(IPD_4m)[dim(IPD_4m)[1]] <- "NEG"

IPD_5m[11,] <- IPD_5m[11,]+IPD_5m[12,]+IPD_5m[13,]
rownames(IPD_5m)[11] <- "15BC"
IPD_5m <- IPD_5m[-c(12, 13, (dim(IPD_5m)[1]-1)),]
rownames(IPD_5m)[dim(IPD_5m)[1]] <- "NEG"

IPD_1m <- as.data.frame(cbind(rownames(IPD_1m), IPD_1m))
IPD_2m <- as.data.frame(cbind(rownames(IPD_2m), IPD_2m))
IPD_3m <- as.data.frame(cbind(rownames(IPD_3m), IPD_3m))
IPD_4m <- as.data.frame(cbind(rownames(IPD_4m), IPD_4m))
IPD_5m <- as.data.frame(cbind(rownames(IPD_5m), IPD_5m))
colnames(IPD_1m) <- c("st", "t0_g1","t1_g1","t2_g1","t3_g1","t4_g1","t5_g1","t6_g1")
colnames(IPD_2m) <- c("st", "t0_g2","t1_g2","t2_g2","t3_g2","t4_g2","t5_g2","t6_g2")
colnames(IPD_3m) <- c("st", "t0_g3","t1_g3","t2_g3","t3_g3","t4_g3","t5_g3","t6_g3")
colnames(IPD_4m) <- c("st", "t0_g4","t1_g4","t2_g4","t3_g4","t4_g4","t5_g4","t6_g4")
colnames(IPD_5m) <- c("st", "t0_g5","t1_g5","t2_g5","t3_g5","t4_g5","t5_g5","t6_g5")

foo3 <- full_join(foo2, IPD_1m, by="st")
foo3 <- full_join(foo3, IPD_2m, by="st")
foo3 <- full_join(foo3, IPD_3m, by="st")
foo3 <- full_join(foo3, IPD_4m, by="st")
foo3 <- full_join(foo3, IPD_5m, by="st")
foo3[is.na(foo3)] <- 0
foo4 <- foo3[order(foo3$st),]
IPD_st <- foo4
IPD.stm     <- matrix(as.numeric(as.character(unlist(IPD_st[,2:36]))), nrow=75, ncol=35, byrow = F)
IPD_st_all <- rowSums(IPD.stm)
IPD_st_all <- as.data.frame(cbind(IPD_st$st, IPD_st_all))
IPD_st_NVTs <- IPD_st_all[-c(1,12,21,22,25,33,43,57,60,61,62,69,74),]


foo5 <- data.frame(foo4$st)
colnames(foo5) <- "st"
IPD_u5     <- full_join(foo5, IPD_1m, by="st")
IPD_5_17   <- full_join(foo5, IPD_2m, by="st")
IPD_18_39  <- full_join(foo5, IPD_3m, by="st")
IPD_40_64  <- full_join(foo5, IPD_4m, by="st")
IPD_65plus <- full_join(foo5, IPD_5m, by="st")

IPD_u5[is.na(IPD_u5)]         <- 0
IPD_5_17[is.na(IPD_5_17)]     <- 0
IPD_18_39[is.na(IPD_18_39)]   <- 0
IPD_40_64[is.na(IPD_40_64)]   <- 0
IPD_65plus[is.na(IPD_65plus)] <- 0

IPD_u5     <- IPD_u5[-75,]
IPD_5_17   <- IPD_5_17[-75,]
IPD_18_39  <- IPD_18_39[-75,]
IPD_40_64  <- IPD_40_64[-75,]
IPD_65plus <- IPD_65plus[-75,]

IPD_all <- cbind(IPD_u5[,1],
                 c(as.numeric(as.character(IPD_u5[,2]))+as.numeric(as.character(IPD_5_17[,2]))+as.numeric(as.character(IPD_18_39[,2]))+as.numeric(as.character(IPD_40_64[,2]))+as.numeric(as.character(IPD_65plus[,2]))),
                 c(as.numeric(as.character(IPD_u5[,3]))+as.numeric(as.character(IPD_5_17[,3]))+as.numeric(as.character(IPD_18_39[,3]))+as.numeric(as.character(IPD_40_64[,3]))+as.numeric(as.character(IPD_65plus[,3]))),
                 c(as.numeric(as.character(IPD_u5[,4]))+as.numeric(as.character(IPD_5_17[,4]))+as.numeric(as.character(IPD_18_39[,4]))+as.numeric(as.character(IPD_40_64[,4]))+as.numeric(as.character(IPD_65plus[,4]))),
                 c(as.numeric(as.character(IPD_u5[,5]))+as.numeric(as.character(IPD_5_17[,5]))+as.numeric(as.character(IPD_18_39[,5]))+as.numeric(as.character(IPD_40_64[,5]))+as.numeric(as.character(IPD_65plus[,5]))),
                 c(as.numeric(as.character(IPD_u5[,6]))+as.numeric(as.character(IPD_5_17[,6]))+as.numeric(as.character(IPD_18_39[,6]))+as.numeric(as.character(IPD_40_64[,6]))+as.numeric(as.character(IPD_65plus[,6]))),
                 c(as.numeric(as.character(IPD_u5[,7]))+as.numeric(as.character(IPD_5_17[,7]))+as.numeric(as.character(IPD_18_39[,7]))+as.numeric(as.character(IPD_40_64[,7]))+as.numeric(as.character(IPD_65plus[,7]))),
                 c(as.numeric(as.character(IPD_u5[,8]))+as.numeric(as.character(IPD_5_17[,8]))+as.numeric(as.character(IPD_18_39[,8]))+as.numeric(as.character(IPD_40_64[,8]))+as.numeric(as.character(IPD_65plus[,8])))
                 )

#1,12,21,22,25,33,43,57,60,61,62,69,74
IPD.NVTs_u5     <- IPD_u5[-c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.NVTs_5_17   <- IPD_5_17[-c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.NVTs_18_39  <- IPD_18_39[-c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.NVTs_40_64  <- IPD_40_64[-c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.NVTs_65plus <- IPD_65plus[-c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.NVTs_all    <- IPD_all[-c(1,12,21,22,25,33,43,57,60,61,62,69,74),] 
#IPD.VTs <- IPD[c(1,12,19,20,23,30,38,51,54,55,56,62,67),]
IPD.VTs_u5     <- IPD_u5[c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.VTs_5_17   <- IPD_5_17[c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.VTs_18_39  <- IPD_18_39[c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.VTs_40_64  <- IPD_40_64[c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.VTs_65plus <- IPD_65plus[c(1,12,21,22,25,33,43,57,60,61,62,69,74),]
IPD.VTs_all    <- IPD_all[c(1,12,21,22,25,33,43,57,60,61,62,69,74),]

IPD.NVTs_u5     <- matrix(as.numeric(as.character(unlist(IPD.NVTs_u5[,2:8]))), nrow=61, ncol=7, byrow = F)
IPD.NVTs_5_17   <- matrix(as.numeric(as.character(unlist(IPD.NVTs_5_17[,2:8]))), nrow=61, ncol=7, byrow = F)
IPD.NVTs_18_39  <- matrix(as.numeric(as.character(unlist(IPD.NVTs_18_39[,2:8]))), nrow=61, ncol=7, byrow = F)
IPD.NVTs_40_64  <- matrix(as.numeric(as.character(unlist(IPD.NVTs_40_64[,2:8]))), nrow=61, ncol=7, byrow = F)
IPD.NVTs_65plus <- matrix(as.numeric(as.character(unlist(IPD.NVTs_65plus[,2:8]))), nrow=61, ncol=7, byrow = F)
IPD.NVTs_all    <- matrix(as.numeric(as.character(unlist(IPD.NVTs_all[,2:8]))), nrow=61, ncol=7, byrow = F)

IPD_st_NVTs_vec <- as.numeric(as.character(IPD_st_NVTs$IPD_st_all))


prevratiomatx <- readRDS("~/Desktop/pneumo/prevalence ratio array.rds")


# prevmatxlist <- NULL
# for (i in 1:30000) {
#   prevmatx <- matrix(prevsamps[i,], ncol = 7, byrow = TRUE)
#   prevmatxlist <- array(c(prevmatxlist, prevmatx), dim = c(55,7,i))
#   print(i)
# }
# 
# 
# curr.sts<- df1
# names(curr.sts) <- "st"

foo6 <- data.frame((as.character(foo5[-c(1,12,21,22,25,33,43,57,60,61,62,69,74,75),])))
names(foo6) <- "st"
# 
# 
# prevratiomatx <- NULL
# for (i in 1:30000){
#   prevs <- matrix(as.numeric(as.character(unlist(prevmatxlist[,,i][-55,]))), nrow=54, ncol=7, byrow = F)
#   prevs.df <- data.frame(cbind(curr.sts, prevs))
#   all.sts <- full_join(foo6, prevs.df, "st")
#   all.sts <- all.sts[order(all.sts$st),]
#   all.sts[,2][is.na(all.sts[,2])] <- mean(all.sts[,2], na.rm = T)
#   all.sts[,3][is.na(all.sts[,3])] <- mean(all.sts[,3], na.rm = T)
#   all.sts[,4][is.na(all.sts[,4])] <- mean(all.sts[,4], na.rm = T)
#   all.sts[,5][is.na(all.sts[,5])] <- mean(all.sts[,5], na.rm = T)
#   all.sts[,6][is.na(all.sts[,6])] <- mean(all.sts[,6], na.rm = T)
#   all.sts[,7][is.na(all.sts[,7])] <- mean(all.sts[,7], na.rm = T)
#   all.sts[,8][is.na(all.sts[,8])] <- mean(all.sts[,8], na.rm = T)
#   prevs.all <- all.sts[,-1]
#   prevratio <- matrix(NA, nrow=61, ncol=7)
#   for(s in 1:61){
#     for (p in 1:7){
#       prevratio[s,p] <- prevs.all[s,p]/prevs.all[s,1]
#     }
#   }
#   prevratiomatx <- array(c(prevratiomatx, prevratio), dim = c(61,7,i))
#   print(i)
# }

# saveRDS(prevratiomatx,"~/Desktop/pneumo/prevalence ratio array.rds")

IPD.VTs_u5     <- matrix(as.numeric(as.character(unlist(IPD.VTs_u5[,-1]))), nrow=13, ncol=7, byrow = F)
IPD.VTs_5_17   <- matrix(as.numeric(as.character(unlist(IPD.VTs_5_17[,-1]))), nrow=13, ncol=7, byrow = F)
IPD.VTs_18_39  <- matrix(as.numeric(as.character(unlist(IPD.VTs_18_39[,-1]))), nrow=13, ncol=7, byrow = F)
IPD.VTs_40_64  <- matrix(as.numeric(as.character(unlist(IPD.VTs_40_64[,-1]))), nrow=13, ncol=7, byrow = F)
IPD.VTs_65plus <- matrix(as.numeric(as.character(unlist(IPD.VTs_65plus[,-1]))), nrow=13, ncol=7, byrow = F)
IPD.VTs_all    <- matrix(as.numeric(as.character(unlist(IPD.VTs_all[,-1]))), nrow=13, ncol=7, byrow = F)


obs.IPD_VT_u5 <- IPD.VTs_u5
  obs.IPD_VT_5_17 <- IPD.VTs_5_17
  obs.IPD_VT_18_39 <- IPD.VTs_18_39
  obs.IPD_VT_40_64 <- IPD.VTs_40_64
  obs.IPD_VT_65plus <- IPD.VTs_65plus
  obs.IPD_VT_all <- IPD.VTs_all
obs.IPD_NVT_u5 <- IPD.NVTs_u5
  obs.IPD_NVT_5_17 <- IPD.NVTs_5_17
  obs.IPD_NVT_18_39 <- IPD.NVTs_18_39
  obs.IPD_NVT_40_64 <- IPD.NVTs_40_64
  obs.IPD_NVT_65plus <- IPD.NVTs_65plus
  obs.IPD_NVT_all <- IPD.NVTs_all

exp.cases_u5     <- array(NA,dim = c(61,7,30000))
exp.cases_5_17   <- array(NA,dim = c(61,7,30000))
exp.cases_18_39  <- array(NA,dim = c(61,7,30000))
exp.cases_40_64  <- array(NA,dim = c(61,7,30000))
exp.cases_65plus <- array(NA,dim = c(61,7,30000))
exp.cases_all    <- array(NA,dim = c(61,7,30000))
for (i in 1:30000){
  for (st in 1:61) {
    for (t in 1:7) {
      exp.cases_u5[st,t,i]     <- obs.IPD_NVT_u5[st,t]    /prevratiomatx[st,t,i]
      exp.cases_5_17[st,t,i]   <- obs.IPD_NVT_5_17[st,t]  /prevratiomatx[st,t,i]
      exp.cases_18_39[st,t,i]  <- obs.IPD_NVT_18_39[st,t] /prevratiomatx[st,t,i]
      exp.cases_40_64[st,t,i]  <- obs.IPD_NVT_40_64[st,t] /prevratiomatx[st,t,i]
      exp.cases_65plus[st,t,i] <- obs.IPD_NVT_65plus[st,t]/prevratiomatx[st,t,i]
      exp.cases_all[st,t,i]    <- obs.IPD_NVT_all[st,t]   /prevratiomatx[st,t,i]
    
    }
  }
}


tot.exp.cases_u5 <- tot.exp.cases_5_17 <- tot.exp.cases_18_39 <- tot.exp.cases_40_64 <- tot.exp.cases_65plus <- tot.exp.cases_all <- matrix(NA, ncol=7, nrow=30000)

tot.obs.IPD_VT_u5     <- colSums(obs.IPD_VT_u5)
tot.obs.IPD_VT_5_17   <- colSums(obs.IPD_VT_5_17)
tot.obs.IPD_VT_18_39  <- colSums(obs.IPD_VT_18_39)
tot.obs.IPD_VT_40_64  <- colSums(obs.IPD_VT_40_64)
tot.obs.IPD_VT_65plus <- colSums(obs.IPD_VT_65plus)
tot.obs.IPD_VT_all    <- colSums(obs.IPD_VT_all)

tot.obs.IPD_NVT_u5     <- colSums(obs.IPD_NVT_u5)
tot.obs.IPD_NVT_5_17   <- colSums(obs.IPD_NVT_5_17)
tot.obs.IPD_NVT_18_39  <- colSums(obs.IPD_NVT_18_39)
tot.obs.IPD_NVT_40_64  <- colSums(obs.IPD_NVT_40_64)
tot.obs.IPD_NVT_65plus <- colSums(obs.IPD_NVT_65plus)
tot.obs.IPD_NVT_all    <- colSums(obs.IPD_NVT_all)

for (i in 1:30000){
    tot.exp.cases_u5[i,]     <- colSums(exp.cases_u5[,,i])
    tot.exp.cases_5_17[i,]   <- colSums(exp.cases_5_17[,,i])
    tot.exp.cases_18_39[i,]  <- colSums(exp.cases_18_39[,,i])
    tot.exp.cases_40_64[i,]  <- colSums(exp.cases_40_64[,,i])
    tot.exp.cases_65plus[i,] <- colSums(exp.cases_65plus[,,i])
    tot.exp.cases_all[i,]    <- colSums(exp.cases_all[,,i])
}



#################################All ages, times, iterations (NO ST)############################
tot.obs.IPD_NVT_u5.2 <- matrix(rep(tot.obs.IPD_NVT_u5,30000), ncol=7, nrow=30000, byrow = T)
tot.obs.IPD_NVT_5_17.2 <- matrix(rep(tot.obs.IPD_NVT_5_17,30000), ncol=7, nrow=30000, byrow = T)
tot.obs.IPD_NVT_18_39.2 <- matrix(rep(tot.obs.IPD_NVT_18_39,30000), ncol=7, nrow=30000, byrow = T)
tot.obs.IPD_NVT_40_64.2 <- matrix(rep(tot.obs.IPD_NVT_40_64,30000), ncol=7, nrow=30000, byrow = T)
tot.obs.IPD_NVT_65plus.2 <- matrix(rep(tot.obs.IPD_NVT_65plus,30000), ncol=7, nrow=30000, byrow = T)
tot.obs.IPD_NVT_all.2 <- matrix(rep(tot.obs.IPD_NVT_all,30000), ncol=7, nrow=30000, byrow = T)

cases.gained_u5     <- tot.obs.IPD_NVT_u5.2   - tot.exp.cases_u5
cases.gained_5_17   <- tot.obs.IPD_NVT_5_17.2  - tot.exp.cases_5_17
cases.gained_18_39  <- tot.obs.IPD_NVT_18_39.2 - tot.exp.cases_18_39
cases.gained_40_64  <- tot.obs.IPD_NVT_40_64.2 - tot.exp.cases_40_64
cases.gained_65plus <- tot.obs.IPD_NVT_65plus.2 - tot.exp.cases_65plus
cases.gained_all    <- tot.obs.IPD_NVT_all.2 - tot.exp.cases_all
tot.exp.cases_u5[cases.gained_u5<0] <- tot.obs.IPD_NVT_u5.2[cases.gained_u5<0]
tot.exp.cases_5_17[cases.gained_5_17<0] <- tot.obs.IPD_NVT_5_17.2[cases.gained_5_17<0]
tot.exp.cases_18_39[cases.gained_18_39<0] <- tot.obs.IPD_NVT_18_39.2[cases.gained_18_39<0]
tot.exp.cases_40_64[cases.gained_40_64<0] <- tot.obs.IPD_NVT_40_64.2[cases.gained_40_64<0]
tot.exp.cases_65plus[cases.gained_65plus<0] <- tot.obs.IPD_NVT_65plus.2[cases.gained_65plus<0]
tot.exp.cases_all[cases.gained_all<0] <- tot.obs.IPD_NVT_all.2[cases.gained_all<0]
cases.gained_u5[cases.gained_u5<0] <- 0
cases.gained_5_17[cases.gained_5_17<0] <- 0
cases.gained_18_39[cases.gained_18_39<0] <- 0
cases.gained_40_64[cases.gained_40_64<0] <- 0
cases.gained_65plus[cases.gained_65plus<0] <- 0
cases.gained_all[cases.gained_all<0] <- 0


#################################All ages, times, serotypes, iterations############################
cases.gained_u5_st     <- array(rep(c(obs.IPD_NVT_u5),30000),dim=c(61,7,30000)) - exp.cases_u5
cases.gained_5_17_st   <- array(rep(c(obs.IPD_NVT_5_17),30000),dim=c(61,7,30000)) - exp.cases_5_17
cases.gained_18_39_st  <- array(rep(c(obs.IPD_NVT_18_39),30000),dim=c(61,7,30000)) - exp.cases_18_39
cases.gained_40_64_st  <- array(rep(c(obs.IPD_NVT_40_64),30000),dim=c(61,7,30000)) - exp.cases_40_64
cases.gained_65plus_st <- array(rep(c(obs.IPD_NVT_65plus),30000),dim=c(61,7,30000)) - exp.cases_65plus
cases.gained_all_st    <- array(rep(c(obs.IPD_NVT_all),30000),dim=c(61,7,30000)) - exp.cases_all
exp.cases_u5[cases.gained_u5_st<0] <- (array(rep(c(obs.IPD_NVT_u5),30000),dim=c(61,7,30000)))[cases.gained_u5_st<0]
exp.cases_5_17[cases.gained_5_17_st<0] <- (array(rep(c(obs.IPD_NVT_5_17),30000),dim=c(61,7,30000)))[cases.gained_5_17_st<0] 
exp.cases_18_39[cases.gained_18_39_st<0] <- (array(rep(c(obs.IPD_NVT_18_39),30000),dim=c(61,7,30000)))[cases.gained_18_39_st<0]
exp.cases_40_64[cases.gained_40_64_st<0] <- (array(rep(c(obs.IPD_NVT_40_64),30000),dim=c(61,7,30000)))[cases.gained_40_64_st<0]
exp.cases_65plus[cases.gained_65plus_st<0] <- (array(rep(c(obs.IPD_NVT_65plus),30000),dim=c(61,7,30000)))[cases.gained_65plus_st<0]
exp.cases_all[cases.gained_all_st<0] <- (array(rep(c(obs.IPD_NVT_all),30000),dim=c(61,7,30000)))[cases.gained_all_st<0]
cases.gained_u5_st[cases.gained_u5_st<0] <- 0
cases.gained_5_17_st[cases.gained_5_17_st<0] <- 0
cases.gained_18_39_st[cases.gained_18_39_st<0] <- 0
cases.gained_40_64_st[cases.gained_40_64_st<0] <- 0
cases.gained_65plus_st[cases.gained_65plus_st<0] <- 0
cases.gained_all_st[cases.gained_all_st<0] <- 0



#################################All ages, serotypes, iterations (NO TIME)############################
cases.gained_u5_st.all <- cases.gained_5_17_st.all <- cases.gained_18_39_st.all <- cases.gained_40_64_st.all <- cases.gained_65plus_st.all <- cases.gained_all_st.all <- exp.cases_st.all <-  array(NA,dim = c(1,61,30000))
for (i in 1:30000) {
  cases.gained_u5_st.all[,,i] <- rowSums(cases.gained_u5_st[,,i])
  cases.gained_5_17_st.all[,,i] <- rowSums(cases.gained_5_17_st[,,i])
  cases.gained_18_39_st.all[,,i] <- rowSums(cases.gained_18_39_st[,,i])
  cases.gained_40_64_st.all[,,i] <- rowSums(cases.gained_40_64_st[,,i])
  cases.gained_65plus_st.all[,,i] <- rowSums(cases.gained_65plus_st[,,i])
  cases.gained_all_st.all[,,i] <- rowSums(cases.gained_all_st[,,i])
  exp.cases_st.all[,,i] <- rowSums(exp.cases_all[,,i])
  
}

cases.gained_time_st.all <-  array(NA,dim = c(1,7,30000))
cases.gained_u5time_st.all <-  array(NA,dim = c(1,7,30000))
cases.gained_5_17time_st.all <-  array(NA,dim = c(1,7,30000))
cases.gained_18_39time_st.all <-  array(NA,dim = c(1,7,30000))
cases.gained_40_64time_st.all <-  array(NA,dim = c(1,7,30000))
cases.gained_65plustime_st.all <-  array(NA,dim = c(1,7,30000))
cases.gained_alltime_st.all <-  array(NA,dim = c(1,7,30000))
for (i in 1:30000) {
  cases.gained_time_st.all[,,i] <- colSums(cases.gained_all_st[,,i])
  cases.gained_u5time_st.all[,,i] <- colSums(cases.gained_u5_st[,,i])
  cases.gained_5_17time_st.all[,,i] <- colSums(cases.gained_5_17_st[,,i])
  cases.gained_18_39time_st.all[,,i] <- colSums(cases.gained_18_39_st[,,i])
  cases.gained_40_64time_st.all[,,i] <- colSums(cases.gained_40_64_st[,,i])
  cases.gained_65plustime_st.all[,,i] <- colSums(cases.gained_65plus_st[,,i])
  cases.gained_alltime_st.all[,,i] <- colSums(cases.gained_all_st[,,i])
  
}


cases.gained.u5.st.quants <- cases.gained.5_17.st.quants <- cases.gained.18_39.st.quants <- cases.gained.40_64.st.quants <- cases.gained.65plus.st.quants <- cases.gained.all.st.quants <- exp.cases_st.all.st.quants <- matrix(NA, ncol = 3, nrow=61)
for (st in 1:61){
  cases.gained.u5.st.quants[st,] <- quantile(cases.gained_u5_st.all[1,st,], probs = c(.025,.5,.975))
  cases.gained.5_17.st.quants[st,] <- quantile(cases.gained_5_17_st.all[1,st,], probs = c(.025,.5,.975))
  cases.gained.18_39.st.quants[st,] <- quantile(cases.gained_18_39_st.all[1,st,], probs = c(.025,.5,.975))
  cases.gained.40_64.st.quants[st,] <- quantile(cases.gained_40_64_st.all[1,st,], probs = c(.025,.5,.975))
  cases.gained.65plus.st.quants[st,] <- quantile(cases.gained_65plus_st.all[1,st,], probs = c(.025,.5,.975))
  cases.gained.all.st.quants[st,] <- quantile(cases.gained_all_st.all[1,st,], probs = c(.025,.5,.975))
  exp.cases_st.all.st.quants[st,] <- quantile(exp.cases_st.all[1,st,], probs = c(.025,.5,.975))
}

rownames(cases.gained.u5.st.quants) <- rownames(cases.gained.5_17.st.quants) <- rownames(cases.gained.18_39.st.quants) <- rownames(cases.gained.40_64.st.quants) <- rownames(cases.gained.65plus.st.quants) <- rownames(cases.gained.all.st.quants) <- rownames(exp.cases_st.all.st.quants) <-  unlist(foo6)
colnames(cases.gained.u5.st.quants) <- colnames(cases.gained.5_17.st.quants) <- colnames(cases.gained.18_39.st.quants) <- colnames(cases.gained.40_64.st.quants) <- colnames(cases.gained.65plus.st.quants) <- colnames(cases.gained.all.st.quants) <- colnames(exp.cases_st.all.st.quants)<- c("lCI","median","uCI")
cases.gained.u5.st.quants[order(cases.gained.u5.st.quants[,2], decreasing = T),]
cases.gained.5_17.st.quants[order(cases.gained.5_17.st.quants[,2], decreasing = T),]
cases.gained.18_39.st.quants[order(cases.gained.18_39.st.quants[,2], decreasing = T),]
cases.gained.40_64.st.quants[order(cases.gained.40_64.st.quants[,2], decreasing = T),]
cases.gained.65plus.st.quants[order(cases.gained.65plus.st.quants[,2], decreasing = T),]
cases.gained.all.st.quants[order(cases.gained.all.st.quants[,2], decreasing = T),]


#write into formatted data for AJE
cases.gained.u5.st.quants <- as.data.frame(cases.gained.u5.st.quants)
cases.gained.u5.st.quants$CI <- with(cases.gained.u5.st.quants, sprintf("%.0f, %.0f",lCI,uCI))
cases.gained.u5.st.quants$Serotype <- unlist(foo6)
colnames(cases.gained.u5.st.quants) <- c("lCI","median","uCI", "CI","Serotype")
cases.gained.u5.st.quants2 <- as.data.frame(left_join(as.data.frame(df1), as.data.frame(cases.gained.u5.st.quants), "Serotype"))
write.csv(cases.gained.u5.st.quants2,"~/Desktop/pneumo.cases.gained.u5.st.quants.AJE.csv")

cases.gained.5_17.st.quants <- as.data.frame(cases.gained.5_17.st.quants)
cases.gained.5_17.st.quants$CI <- with(cases.gained.5_17.st.quants, sprintf("%.0f, %.0f",lCI,uCI))
cases.gained.5_17.st.quants$Serotype <- unlist(foo6)
colnames(cases.gained.5_17.st.quants) <- c("lCI","median","uCI", "CI","Serotype")
cases.gained.5_17.st.quants2 <- as.data.frame(left_join(as.data.frame(df1), as.data.frame(cases.gained.5_17.st.quants), "Serotype"))
write.csv(cases.gained.5_17.st.quants2,"~/Desktop/pneumo.cases.gained.5_17.st.quants.AJE.csv")

cases.gained.18_39.st.quants <- as.data.frame(cases.gained.18_39.st.quants)
cases.gained.18_39.st.quants$CI <- with(cases.gained.18_39.st.quants, sprintf("%.0f, %.0f",lCI,uCI))
cases.gained.18_39.st.quants$Serotype <- unlist(foo6)
colnames(cases.gained.18_39.st.quants) <- c("lCI","median","uCI", "CI","Serotype")
cases.gained.18_39.st.quants2 <- as.data.frame(left_join(as.data.frame(df1), as.data.frame(cases.gained.18_39.st.quants), "Serotype"))
write.csv(cases.gained.18_39.st.quants2,"~/Desktop/pneumo.cases.gained.18_39.st.quants.AJE.csv")

cases.gained.40_64.st.quants <- as.data.frame(cases.gained.40_64.st.quants)
cases.gained.40_64.st.quants$CI <- with(cases.gained.40_64.st.quants, sprintf("%.0f, %.0f",lCI,uCI))
cases.gained.40_64.st.quants$Serotype <- unlist(foo6)
colnames(cases.gained.40_64.st.quants) <- c("lCI","median","uCI", "CI","Serotype")
cases.gained.40_64.st.quants2 <- as.data.frame(left_join(as.data.frame(df1), as.data.frame(cases.gained.40_64.st.quants), "Serotype"))
write.csv(cases.gained.40_64.st.quants2,"~/Desktop/pneumo.cases.gained.40_64.st.quants.AJE.csv")

cases.gained.65plus.st.quants <- as.data.frame(cases.gained.65plus.st.quants)
cases.gained.65plus.st.quants$CI <- with(cases.gained.65plus.st.quants, sprintf("%.0f, %.0f",lCI,uCI))
cases.gained.65plus.st.quants$Serotype <- unlist(foo6)
colnames(cases.gained.65plus.st.quants) <- c("lCI","median","uCI", "CI","Serotype")
cases.gained.65plus.st.quants2 <- as.data.frame(left_join(as.data.frame(df1), as.data.frame(cases.gained.65plus.st.quants), "Serotype"))
write.csv(cases.gained.65plus.st.quants2,"~/Desktop/pneumo.cases.gained.65plus.st.quants.AJE.csv")

cases.gained.all.st.quants <- as.data.frame(cases.gained.all.st.quants)
cases.gained.all.st.quants$CI <- with(cases.gained.all.st.quants, sprintf("%.0f, %.0f",lCI,uCI))
cases.gained.all.st.quants$Serotype <- unlist(foo6)
colnames(cases.gained.all.st.quants) <- c("lCI","median","uCI", "CI","Serotype")
cases.gained.all.st.quants2 <- as.data.frame(left_join(as.data.frame(df1), as.data.frame(cases.gained.all.st.quants), "Serotype"))
write.csv(cases.gained.all.st.quants2,"~/Desktop/pneumo.cases.gained.all.st.quants.AJE.csv")


#################################OVERALL############################
tot.cases.gained <- NULL
for (i in 1:30000){
  tot.cases.gained[i] <- sum(cases.gained_all_st.all[,,i])
  
}
quantile(tot.cases.gained, probs = c(.025,.5,.975))
df2 <- df1

#################################by time overall############################
cases.gained.y0.quants <- quantile(cases.gained_time_st.all[1,1,], probs = c(.025,.5,.975))
cases.gained.y1.quants <- quantile(cases.gained_time_st.all[1,2,], probs = c(.025,.5,.975))
cases.gained.y2.quants <- quantile(cases.gained_time_st.all[1,3,], probs = c(.025,.5,.975))
cases.gained.y3.quants <- quantile(cases.gained_time_st.all[1,4,], probs = c(.025,.5,.975))
cases.gained.y4.quants <- quantile(cases.gained_time_st.all[1,5,], probs = c(.025,.5,.975))
cases.gained.y5.quants <- quantile(cases.gained_time_st.all[1,6,], probs = c(.025,.5,.975))
cases.gained.y6.quants <- quantile(cases.gained_time_st.all[1,7,], probs = c(.025,.5,.975))

#################################by age overall############################
  cases.gained.u5.sums <- cases.gained.5_17.sums <- cases.gained.18_39.sums <- cases.gained.40_64.sums <- cases.gained.65plus.sums <- NULL
  for (i in 1:30000){
  cases.gained.u5.sums[i] <- sum(cases.gained_u5_st.all[,,i])
  cases.gained.5_17.sums[i] <- sum(cases.gained_5_17_st.all[,,i])
  cases.gained.18_39.sums[i] <- sum(cases.gained_18_39_st.all[,,i])
  cases.gained.40_64.sums[i] <- sum(cases.gained_40_64_st.all[,,i])
  cases.gained.65plus.sums[i] <- sum(cases.gained_65plus_st.all[,,i])
  }
  
  cases.gained.u5.quants <- quantile(cases.gained.u5.sums, probs = c(.025,.5,.975))
  cases.gained.5_17.quants <- quantile(cases.gained.5_17.sums, probs = c(.025,.5,.975))
  cases.gained.18_39.quants <- quantile(cases.gained.18_39.sums, probs = c(.025,.5,.975))
  cases.gained.40_64.quants <- quantile(cases.gained.40_64.sums, probs = c(.025,.5,.975))
  cases.gained.65plus.quants <- quantile(cases.gained.65plus.sums, probs = c(.025,.5,.975))
  
  cases.gained_u5.mcmc <- cases.gained_5_17.mcmc <- cases.gained_18_39.mcmc <- cases.gained_40_64.mcmc <- cases.gained_65plus.mcmc <- cases.gained_all.mcmc <- NULL
  tot.obs.IPD_VT_u5.mcmc <- tot.obs.IPD_VT_5_17.mcmc <- tot.obs.IPD_VT_18_39.mcmc <- tot.obs.IPD_VT_40_64.mcmc <- tot.obs.IPD_VT_65plus.mcmc <- tot.obs.IPD_VT_all.mcmc <- NULL
  tot.exp.cases_u5.mcmc <- tot.exp.cases_5_17.mcmc <- tot.exp.cases_18_39.mcmc <- tot.exp.cases_40_64.mcmc <- tot.exp.cases_65plus.mcmc <- tot.exp.cases_all.mcmc <- NULL
  for (t in 1:7){
    cases.gained_u5.mcmc[t] <- median(cases.gained_u5time_st.all[1,t,])
    cases.gained_5_17.mcmc[t] <- median(cases.gained_5_17time_st.all[1,t,])
    cases.gained_18_39.mcmc[t] <- median(cases.gained_18_39time_st.all[1,t,])
    cases.gained_40_64.mcmc[t] <- median(cases.gained_40_64time_st.all[1,t,])
    cases.gained_65plus.mcmc[t] <- median(cases.gained_65plustime_st.all[1,t,])
    cases.gained_all.mcmc[t] <- median(cases.gained_alltime_st.all[1,t,])
    
    tot.obs.IPD_VT_u5.mcmc[t] <-tot.obs.IPD_VT_u5[t]
    tot.obs.IPD_VT_5_17.mcmc[t] <-tot.obs.IPD_VT_5_17[t]
    tot.obs.IPD_VT_18_39.mcmc[t] <-tot.obs.IPD_VT_18_39[t]
    tot.obs.IPD_VT_40_64.mcmc[t] <-tot.obs.IPD_VT_40_64[t]
    tot.obs.IPD_VT_65plus.mcmc[t] <-tot.obs.IPD_VT_65plus[t]
    tot.obs.IPD_VT_all.mcmc[t] <-tot.obs.IPD_VT_all[t]
    
    tot.exp.cases_u5.mcmc[t] <- median(tot.exp.cases_u5[,t])
    tot.exp.cases_5_17.mcmc[t] <- median(tot.exp.cases_5_17[,t])
    tot.exp.cases_18_39.mcmc[t] <- median(tot.exp.cases_18_39[,t])
    tot.exp.cases_40_64.mcmc[t] <- median(tot.exp.cases_40_64[,t])
    tot.exp.cases_65plus.mcmc[t] <- median(tot.exp.cases_65plus[,t])
    tot.exp.cases_all.mcmc[t] <- median(tot.exp.cases_all[,t])
  }


exp.cases_u5 <- rep(tot.exp.cases_u5.mcmc[1],7)
unexp.cases_u5 <- tot.exp.cases_u5.mcmc-exp.cases_u5

exp.cases_5_17 <- rep(tot.exp.cases_5_17.mcmc[1],7)
unexp.cases_5_17 <- tot.exp.cases_5_17.mcmc-exp.cases_5_17
unexp.cases_5_17[which(unexp.cases_5_17<0)] <- 0

exp.cases_18_39 <- rep(tot.exp.cases_18_39.mcmc[1],7)
unexp.cases_18_39 <- tot.exp.cases_18_39.mcmc-exp.cases_18_39
unexp.cases_18_39[which(unexp.cases_18_39<0)] <- 0

exp.cases_40_64 <- rep(tot.exp.cases_40_64.mcmc[1],7)
unexp.cases_40_64 <- tot.exp.cases_40_64.mcmc-exp.cases_40_64
unexp.cases_40_64[which(unexp.cases_40_64<0)] <- 0

exp.cases_65plus <- rep(tot.exp.cases_65plus.mcmc[1],7)
unexp.cases_65plus <- tot.exp.cases_65plus.mcmc-exp.cases_65plus
unexp.cases_65plus[which(unexp.cases_65plus<0)] <- 0

exp.cases_all <- rep(tot.exp.cases_all.mcmc[1],7)
unexp.cases_all <- tot.exp.cases_all.mcmc-exp.cases_all
unexp.cases_all[which(unexp.cases_all<0)] <- 0

# saveRDS(obs.IPD_VT_all,"~/Desktop/pneumo/obs VT IPD bytime,st all ages.rds")
# saveRDS(obs.IPD_VT_u5,"~/Desktop/pneumo/obs VT IPD bytime,st u5.rds")
# saveRDS(obs.IPD_VT_5_17,"~/Desktop/pneumo/obs VT IPD bytime,st 5_17.rds")
# saveRDS(obs.IPD_VT_18_39,"~/Desktop/pneumo/obs VT IPD bytime,st 18_39.rds")
# saveRDS(obs.IPD_VT_40_64,"~/Desktop/pneumo/obs VT IPD bytime,st 40_64.rds")
# saveRDS(obs.IPD_VT_65plus,"~/Desktop/pneumo/obs VT IPD bytime,st 65plus.rds")
# 
# saveRDS(exp.cases_all,"~/Desktop/pneumo/expected IPD bytime,st all ages.rds")
# saveRDS(exp.cases_u5,"~/Desktop/pneumo/expected IPD bytime,st u5.rds")
# saveRDS(exp.cases_5_17,"~/Desktop/pneumo/expected IPD bytime,st 5_17.rds")
# saveRDS(exp.cases_18_39,"~/Desktop/pneumo/expected IPD bytime,st 18_39.rds")
# saveRDS(exp.cases_40_64,"~/Desktop/pneumo/expected IPD bytime,st 40_64.rds")
# saveRDS(exp.cases_65plus,"~/Desktop/pneumo/expected IPD bytime,st 65plus.rds")
# 
# saveRDS(cases.gained_all_st,"~/Desktop/pneumo/cases gained IPD bytime,st all ages.rds")
# saveRDS(cases.gained_u5_st,"~/Desktop/pneumo/cases gained IPD bytime,st u5.rds")
# saveRDS(cases.gained_5_17_st,"~/Desktop/pneumo/cases gained IPD bytime,st 5_17.rds")
# saveRDS(cases.gained_18_39_st,"~/Desktop/pneumo/cases gained IPD bytime,st 18_39.rds")
# saveRDS(cases.gained_40_64_st,"~/Desktop/pneumo/cases gained IPD bytime,st 40_64.rds")
# saveRDS(cases.gained_65plus_st,"~/Desktop/pneumo/cases gained IPD bytime,st 65plus.rds")

 
plot.new()
par(mfrow=c(1,6))
par(oma = c(4, 3, 0.5, 0.5)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(1, 1, 0, 0), xpd=TRUE)
counts_allu5 <- matrix(c(tot.obs.IPD_VT_u5.mcmc, exp.cases_u5, cases.gained_u5.mcmc, unexp.cases_u5), nrow=4, byrow = T)
colnames(counts_allu5) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_allu5 ,xlab = '', ylab = '', ylim=c(0,260),
        col=c("darkblue","red", "red", "grey60"),density=c(100,100,50,90)) 
mtext('Under 5 years', side = 3, line = -2, adj = 1, cex=0.9)

counts_all5_17 <- matrix(c(tot.obs.IPD_VT_5_17.mcmc, exp.cases_5_17, cases.gained_5_17.mcmc, unexp.cases_5_17), nrow=4, byrow = T)
colnames(counts_all5_17) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_all5_17 ,xlab = '', ylab = '', ylim=c(0,260), yaxt="n",
        col=c("darkblue","red", "red", "grey60"),density=c(100,100,50,90))  
mtext('5-17 years', side = 3, line = -2, adj = 1, cex=0.9)

counts_all18_39 <- matrix(c(tot.obs.IPD_VT_18_39.mcmc, exp.cases_18_39, cases.gained_18_39.mcmc, unexp.cases_18_39), nrow=4, byrow = T)
colnames(counts_all18_39) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_all18_39 ,xlab = '', ylab = '', ylim=c(0,260), yaxt="n",
        col=c("darkblue","red", "red", "grey60"),density=c(100,100,50,90)) 
mtext('18-39 years', side = 3, line = -2, adj = 1, cex=0.9)

counts_all40_64 <- matrix(c(tot.obs.IPD_VT_40_64.mcmc, exp.cases_40_64, cases.gained_40_64.mcmc, unexp.cases_40_64), nrow=4, byrow = T)
colnames(counts_all40_64) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_all40_64 ,xlab = '', ylab = '', ylim=c(0,260), yaxt="n",
        col=c("darkblue","red", "red", "grey60"),density=c(100,100,50,90)) 
mtext('40-64 years', side = 3, line = -2, adj = 1, cex=0.9)

counts_all65plus <- matrix(c(tot.obs.IPD_VT_65plus.mcmc, exp.cases_65plus, cases.gained_65plus.mcmc, unexp.cases_65plus), nrow=4, byrow = T)
colnames(counts_all65plus) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_all65plus ,xlab = '', ylab = '', ylim=c(0,260), yaxt="n",
        col=c("darkblue","red", "red", "grey60"),density=c(100,100,50,90)) 
mtext('65+ years', side = 3, line = -2, adj = 1, cex=0.9)

mtext('Time', side = 1, outer = TRUE, line = 2.2)
mtext('Observed number of IPD cases', side = 2, outer = TRUE, line = 1.5)

plot.new()
legend("bottomright",legend = c("Observed VT cases", "Expected NVT cases \nwithout vaccines", 
                                "Estimated NVT cases \ndue to increased carriage","Estimated unexplained \nNVT cases"),
       fill=c("darkblue","red", "red", "grey60"),density=c(100,100,50,90), cex=1.2, y.intersp=2, bty="n")
dev.off()


# black and white

plot.new()
par(mfrow=c(1,6))
par(oma = c(4, 3, 0.5, 0.5)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(1, 1, 0, 0), xpd=TRUE)
counts_allu5 <- matrix(c(tot.obs.IPD_VT_u5.mcmc, exp.cases_u5, cases.gained_u5.mcmc, unexp.cases_u5), nrow=4, byrow = T)
colnames(counts_allu5) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_allu5 ,xlab = '', ylab = '', ylim=c(0,260),
        col=c("grey60",1,1,"grey70") ,angle=c(45, 45, 0, 135), density=c(100,10,30,30)) 
mtext('Under 5 years', side = 3, line = -2, adj = 1, cex=0.9)

counts_all5_17 <- matrix(c(tot.obs.IPD_VT_5_17.mcmc, exp.cases_5_17, cases.gained_5_17.mcmc, unexp.cases_5_17), nrow=4, byrow = T)
colnames(counts_all5_17) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_all5_17 ,xlab = '', ylab = '', ylim=c(0,260), yaxt="n",
        col=c("grey60",1,1,"grey70") ,angle=c(45, 45, 0, 135), density=c(100,10,30,30))  
mtext('5-17 years', side = 3, line = -2, adj = 1, cex=0.9)

counts_all18_39 <- matrix(c(tot.obs.IPD_VT_18_39.mcmc, exp.cases_18_39, cases.gained_18_39.mcmc, unexp.cases_18_39), nrow=4, byrow = T)
colnames(counts_all18_39) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_all18_39 ,xlab = '', ylab = '', ylim=c(0,260), yaxt="n",
        col=c("grey60",1,1,"grey70") ,angle=c(45, 45, 0, 135), density=c(100,10,30,30)) 
mtext('18-39 years', side = 3, line = -2, adj = 1, cex=0.9)

counts_all40_64 <- matrix(c(tot.obs.IPD_VT_40_64.mcmc, exp.cases_40_64, cases.gained_40_64.mcmc, unexp.cases_40_64), nrow=4, byrow = T)
colnames(counts_all40_64) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_all40_64 ,xlab = '', ylab = '', ylim=c(0,260), yaxt="n",
        col=c("grey60",1,1,"grey70") ,angle=c(45, 45, 0, 135), density=c(100,10,30,30)) 
mtext('40-64 years', side = 3, line = -2, adj = 1, cex=0.9)

counts_all65plus <- matrix(c(tot.obs.IPD_VT_65plus.mcmc, exp.cases_65plus, cases.gained_65plus.mcmc, unexp.cases_65plus), nrow=4, byrow = T)
colnames(counts_all65plus) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
barplot(counts_all65plus ,xlab = '', ylab = '', ylim=c(0,260), yaxt="n",
        col=c("grey60",1,1,"grey70") ,angle=c(45, 45, 0, 135), density=c(100,10,30,30)) 
mtext('65+ years', side = 3, line = -2, adj = 1, cex=0.9)

mtext('Time', side = 1, outer = TRUE, line = 2.2)
mtext('Observed number of IPD cases', side = 2, outer = TRUE, line = 1.5)

plot.new()
legend("bottomright",legend = c("Observed VT cases", "Expected NVT cases \nwithout vaccines", 
                                "Estimated NVT cases \ndue to increased carriage","Estimated unexplained \nNVT cases"),
       fill=c("grey60",1,1,"grey70") ,angle=c(45, 45, 0, 135), density=c(100,10,30,30), cex=1.2, y.intersp=2, bty="n")
dev.off()


par(mfrow=c(1,1))
counts_ALL <- matrix(c(tot.obs.IPD_VT_all.mcmc, tot.exp.cases_all.mcmc, cases.gained_all.mcmc, unexp.cases_all), nrow=4, byrow = T)
colnames(counts_ALL) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
#rownames(counts_ALL) <- c("Observed cases", "Cases gained through serotype replacement")
barplot(counts_ALL,xlab = 'Time', ylab = 'Observed number of IPD cases', ylim=c(0,750),
        col=c("grey60",1,1,"grey70") ,angle=c(45, 45, 0, 135), density=c(100,10,30,30),
        legend.text = c("Observed VT cases", "Expected NVT cases without vaccines", "Estimated additional cases due to replacement",
                        "Estimated unexplained NVT cases"), 
        args.legend = list(x=8.7, y=750, cex=.7)) 

#title("Cases over study period: All individuals")
title('All ages')
dev.off()


##########################################################################
cases.gained_all_st
obs.IPD_NVT_all


exp.cases_all    <- array(NA,dim = c(61,7,30000))
for (i in 1:30000){
  for (st in 1:61) {
    for (t in 1:7) {
      exp.cases_all[st,t,i]    <- obs.IPD_NVT_all[st,t]   /prevratiomatx[st,t,i]
    }
  }
}


cases.gained_all_st    <- array(rep(c(obs.IPD_NVT_all),30000),dim=c(61,7,30000)) - exp.cases_all
exp.cases_all[cases.gained_all_st<0] <- (array(rep(c(obs.IPD_NVT_all),30000),dim=c(61,7,30000)))[cases.gained_all_st<0]
cases.gained_all_st[cases.gained_all_st<0] <- 0


# tot.obs.IPD_VT_st.mcmc <- matrix(NA, nrow = 13, ncol = 7)
tot.exp.cases_st.mcmc <- matrix(NA, nrow = 61, ncol = 7)
tot.unexp.gained_st.mcmc <- matrix(NA, nrow = 61, ncol = 7)
tot.exp.cases_st.mcmc2 <- matrix(NA, nrow = 61, ncol = 7)
tot.cases.gained_st.mcmc <- matrix(NA, nrow = 61, ncol = 7)
for (t in 1:7){
  for (nvt in 1:61){
  tot.exp.cases_st.mcmc[nvt,t] <- mean(exp.cases_all[nvt,t,])
  tot.exp.cases_st.mcmc2[nvt,] <- rep(tot.exp.cases_st.mcmc[nvt,1],7)
  tot.unexp.gained_st.mcmc[nvt,t] <- tot.exp.cases_st.mcmc[nvt,t]-tot.exp.cases_st.mcmc2[nvt,t]
  tot.unexp.gained_st.mcmc[nvt,t][which(tot.unexp.gained_st.mcmc[nvt,t]<0)] <- 0
  tot.cases.gained_st.mcmc[nvt,t] <- mean(cases.gained_all_st[nvt,t,])
  
  #mean(exp.cases_all[nvt,t,])
  
} }

# tot.obs.IPD_VT_st.mcmc
tot.cases.gained_st.mcmc
tot.exp.cases_st.mcmc
tot.exp.cases_st.mcmc2
tot.unexp.gained_st.mcmc

zeros <- matrix(0, nrow = 3, ncol = 7)
tot.exp.cases_st.mcmc2 <- rbind(tot.exp.cases_st.mcmc2,zeros)
tot.unexp.gained_st.mcmc <- rbind(tot.unexp.gained_st.mcmc,zeros)
tot.cases.gained_st.mcmc <- rbind(tot.cases.gained_st.mcmc,zeros)

st.names <- as.character(IPD_st_NVTs[-62,1])
study.pd <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")

par(mfrow=c(4,4))
par(oma = c(4, 3, .5, 0.5)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(1, 1, 1, 0))

n=1
for (st in (((16*(n-1))+1:16))){
  counts_ALL <- matrix(c(tot.exp.cases_st.mcmc2[st,], tot.cases.gained_st.mcmc[st,],tot.unexp.gained_st.mcmc[st,]), nrow=3, byrow = T)
  colnames(counts_ALL) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
  barplot(counts_ALL, xaxt='n',yaxt='n', ylim=c(0,70),
          col=c(1,1,"grey70") ,angle=c(45, 0, 135), density=c(10,30,30)) 
  title(st.names[st])
  if (st %in% (((16*(n-1))+13:16)))
    axis(1, col = "grey40", labels = study.pd, at = 1:7)
  if (st %in% c(((16*(n-1))+1), ((16*(n-1))+5), ((16*(n-1))+9), ((16*(n-1))+13)))
    axis(2, col = "grey40", at = seq(0,50, 10))
}
mtext("Year", side = 1, outer = TRUE, cex = 0.7, line = 2.2)
mtext("Observed number of NVT IPD cases", side = 2, outer = TRUE, cex = 0.7, line = 2)

n=2
for (st in (((16*(n-1))+1:16))){
  counts_ALL <- matrix(c(tot.exp.cases_st.mcmc2[st,], tot.cases.gained_st.mcmc[st,],tot.unexp.gained_st.mcmc[st,]), nrow=3, byrow = T)
  colnames(counts_ALL) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
  barplot(counts_ALL, xaxt='n',yaxt='n', ylim=c(0,70),
          col=c(1,1,"grey70") ,angle=c(45, 0, 135), density=c(10,30,30)) 
  title(st.names[st])
  if (st %in% (((16*(n-1))+13:16)))
    axis(1, col = "grey40", labels = study.pd, at = 1:7)
  if (st %in% c(((16*(n-1))+1), ((16*(n-1))+5), ((16*(n-1))+9), ((16*(n-1))+13)))
    axis(2, col = "grey40", at = seq(0,50, 10))
}
mtext("Year", side = 1, outer = TRUE, cex = 0.7, line = 2.2)
mtext("Observed number of NVT IPD cases", side = 2, outer = TRUE, cex = 0.7, line = 2)

n=3
for (st in (((16*(n-1))+1:16))){
  counts_ALL <- matrix(c(tot.exp.cases_st.mcmc2[st,], tot.cases.gained_st.mcmc[st,],tot.unexp.gained_st.mcmc[st,]), nrow=3, byrow = T)
  colnames(counts_ALL) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
  barplot(counts_ALL, xaxt='n',yaxt='n', ylim=c(0,70),
          col=c(1,1,"grey70") ,angle=c(45, 0, 135), density=c(10,30,30)) 
  title(st.names[st])
  if (st %in% (((16*(n-1))+13:16)))
    axis(1, col = "grey40", labels = study.pd, at = 1:7)
  if (st %in% c(((16*(n-1))+1), ((16*(n-1))+5), ((16*(n-1))+9), ((16*(n-1))+13)))
    axis(2, col = "grey40", at = seq(0,50, 10))
}
mtext("Year", side = 1, outer = TRUE, cex = 0.7, line = 2.2)
mtext("Observed number of NVT IPD cases", side = 2, outer = TRUE, cex = 0.7, line = 2)

n=4
for (st in (((16*(n-1))+1:16))){
  counts_ALL <- matrix(c(tot.exp.cases_st.mcmc2[st,], tot.cases.gained_st.mcmc[st,],tot.unexp.gained_st.mcmc[st,]), nrow=3, byrow = T)
  colnames(counts_ALL) <- c("Year 0","Year 1", "Year 2", "Year 3", "Year 4", "Year 5", "Year 6")
  barplot(counts_ALL, xaxt='n',yaxt='n', ylim=c(0,70),
          col=c(1,1,"grey70") ,angle=c(45, 0, 135), density=c(10,30,30)) 
  title(st.names[st])
  if (st %in% (((16*(n-1))+13:16)))
    axis(1, col = "grey40", labels = study.pd, at = 1:7)
  if (st %in% c(((16*(n-1))+1), ((16*(n-1))+5), ((16*(n-1))+9), ((16*(n-1))+13)))
    axis(2, col = "grey40", at = seq(0,50, 10))
}
mtext("Year", side = 1, outer = TRUE, cex = 0.7, line = 2.2)
mtext("Observed number of NVT IPD cases", side = 2, outer = TRUE, cex = 0.7, line = 2)

#######################################
####################################################

# exp.cases_all_st.mcmc <- matrix(NA, nrow = 30000, ncol = 61)
# for (i in 1:30000){
#     exp.cases_all_st.mcmc[i,] <- rowSums(exp.cases_all[,,i])
#     tot.exp.cases_st.mcmc2[nvt,] <- rep(tot.exp.cases_st.mcmc[nvt,1],7)
#     # tot.unexp.gained_st.mcmc[nvt,t] <- tot.exp.cases_st.mcmc2[nvt,t]-tot.exp.cases_st.mcmc[nvt,t]
#     # tot.unexp.gained_st.mcmc[nvt,t][which(tot.unexp.gained_st.mcmc[nvt,t]<0)] <- 0
#     # tot.cases.gained_st.mcmc[nvt,t] <- mean(cases.gained_all_st[nvt,t,])
#     
#     #mean(exp.cases_all[nvt,t,])
#     
#   }
# 
# #tot.obs.IPD_VT_st.mcmc
# tot.cases.gained_st.mcmc
# tot.exp.cases_st.mcmc
# tot.exp.cases_st.mcmc2
# tot.unexp.gained_st.mcmc
#tot.exp.cases_all
#########
exp.cases_const <- array(NA,dim = c(61,7,30000))
for (i in 1:30000){
  exp.cases_const[,,i] <- matrix(rep(exp.cases_all[,1,i], each=7),nrow=61, ncol=7, byrow = T)
}

unexplained.cases <- exp.cases_all-exp.cases_const

neg_to_0 <- function(x){
  newval <- ifelse(x<0, 0, x)
  return(newval)
}


library(tictoc)

tic("sleeping")
print("falling asleep...")
unexplained.cases2 <- array(sapply(unexplained.cases, FUN = neg_to_0),dim = c(61,7,30000)) 
print("...waking up")
toc()

unexplained.cases.tot <- matrix(NA, nrow=61, ncol=30000)
unexplained.cases.tot.neg <- matrix(NA, nrow=61, ncol=30000)
unexplained.cases.all <- rep(NA, 30000)
unexplained.cases.all.neg <- rep(NA, 30000)
for (i in 1:30000){
  unexplained.cases.tot[,i] <- rowSums(unexplained.cases2[,,i])
  unexplained.cases.tot.neg[,i] <- rowSums(unexplained.cases[,,i])
  unexplained.cases.all[i] <- sum(unexplained.cases2[,,i])
  unexplained.cases.all.neg[i] <- sum(unexplained.cases[,,i])
  
}

quantile((unexplained.cases.all), probs=c(.5, .025, .975) )
quantile((unexplained.cases.all.neg), probs=c(.5, .025, .975) )

unexp.cases_quants <- matrix(NA, nrow=61,ncol=3)
unexp.cases.neg_quants <- matrix(NA, nrow=61,ncol=3)
for (v in 1:61){
  unexp.cases_quants[v,] <- quantile(unexplained.cases.tot[v,], probs = c(.5, .025,.975))
  unexp.cases.neg_quants[v,] <- quantile(unexplained.cases.tot.neg[v,], probs = c(.5, .025,.975))
}

quantile((unexplained.cases.all), probs=c(.5, .025, .975) )

rownames(unexp.cases_quants) <- st.names
unexp.cases_quants

rownames(unexp.cases.neg_quants) <- st.names
unexp.cases.neg_quants

write.csv(as.data.frame(unexp.cases_quants), "~/Desktop/pneumo/unexplained cases quants by st.csv")
write.csv(as.data.frame(unexp.cases.neg_quants), "~/Desktop/pneumo/unexplained cases quants by st w neg.csv")

unexp.cases_quants <- as.data.frame(unexp.cases_quants)
names(unexp.cases_quants) <- c("median","lCI","uCI")
unexp.cases_quants$CI <- with(unexp.cases_quants, sprintf("%.0f, %.0f",lCI,uCI))
unexp.cases_quants$Serotype <- unlist(foo6)
#unexp.cases_quants2 <- as.data.frame(left_join(as.data.frame(df1), as.data.frame(unexp.cases_quants), "Serotype"))
write.csv(unexp.cases_quants,"~/Desktop/unexplained cases.quants.AJE.csv")
