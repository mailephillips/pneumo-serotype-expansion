if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

library(rjags)
library(coda)
library(mcmcplots)
library(tidyverse)

set.seed(1234)


ds <- (read.csv(file = "carriage_final_Jewish.csv"))[,2:9]
new.steff <- read.csv("new_st_eff.csv")
st.reg <- read.csv("coef.st.regs.plot.csv")
preds <- (read.csv(file = "carriage_final_Jewish.csv"))[,c(2,10:15)]


jcode13 <-"
model{
for (t in 1:n){
N_carr[t,] ~ dmulti(prev[t,], N_swab[t])  
for (m in 1:k){
numerator[t,m] <- (exp(beta0.st[m] + beta1.st[m]*t1[t] + beta2.st[m]*t2[t]))
}

denominator[t] <- sum(numerator[t,])
for (m in 1:k){
prev[t,m] <- numerator[t,m]/denominator[t]
}}

# for (m in 1:k){
# beta2.st[m]= alpha0.st[m]+ alpha3*dens[m]}

beta0~dnorm(0, 0.0001) # mean for random intercepts
beta1~dnorm(0, 0.0001)
beta2~dnorm(0, 0.0001)# mean for random slopes
sigma0~dunif(0, 100) # SD of intercepts
sigma1~dunif(0, 100) # SD of slopes
sigma2~dunif(0, 100)

rho~dunif(-1, 1) # correlation between intercepts and slopes

Sigma[1, 1] <- sigma0^2  #var-covar matrix for the random effects
Sigma[2, 2] <- sigma2^2
Sigma[1, 2] <- rho*sigma0*sigma2
Sigma[2, 1] <- rho*sigma0*sigma2

InvSigma[1:2, 1:2] <- inverse(Sigma[,])

for (i in 1:54) {
B.hat[i, 1] <- beta0
B.hat[i, 2] <- beta2
B[i, 1:2]~dmnorm(B.hat[i, ], InvSigma[,]) # pairs of correlated random effects
beta0.st[i] <- B[i, 1] # random intercept
beta2.st[i] <- B[i, 2] # random slope
}

for(j in 1:54){
######SEROTYPE INVASIVENESS
# beta0.st[j]~dnorm( beta0 , tau.beta0)
beta1.st[j]~dnorm( beta1 , tau.beta1)
# alpha0.st[j]~dnorm(alpha0, tau.alpha0)
}

for(l in 1:55){
exp.beta0.st[l] <- exp(beta0.st[l])
exp.beta1.st[l] <- exp(beta1.st[l])
exp.beta2.st[l] <- exp(beta2.st[l])
prev.ratio.fin[l] <- prev[7,l]/prev[1,l]
log.prev.ratio.fin[l] <- log(prev.ratio.fin[l])
for(t in 1:7){
log.prevratio[l,t] <- log(prev.ratio[l,t])
prev.ratio[l,t] <- prev[t,l]/prev[1,l]}
}

#Hyperpriors for serotype-level loop
# alpha0~dnorm(0,1e-5)
# alpha3~dnorm(0,1e-5)

sd.beta1~dunif(0,100)
tau.beta1<-1/sd.beta1^2
# sd.alpha0~dunif(0,100)
# tau.alpha0<-1/sd.alpha0^2

# mu.dens~dnorm(0, 1e-5)
# sd.dens~dunif(0,100)
# tau.dens<-1/sd.dens^2

beta0.st[k] <- 0
beta1.st[k] <- 0
beta2.st[k] <- 0


}"


#same serotypes across all variables
new2 <- as.data.frame(new.steff[,2])
names(new2) <- "Pnc_Serotype"
ds1 <- left_join(new2, ds, by="Pnc_Serotype")
ds1 <- rbind(ds1,ds[37,])
ds1[is.na(ds1)] <- 0

#carriage
ds2 <- ds1[,2:8]
N_carr <- t(ds2)
N_swab <- rowSums(N_carr)
n=nrow(N_carr) #7 time periods
k=ncol(N_carr) #37 serotypes (including negative pneumo)

#define cat time 2
t0=c(1,0,0,0,0,0,0)
t1=c(0,1,1,1,1,1,1)
t2=c(0,0,0,1,1,1,1)
#t1=c(0,1,2,3,4,5,6)

# #predictors
# st.reg2 <- cbind(st.reg[,7], st.reg[,2:5])
# names(st.reg2)[1] <- "Serotype"
# st.reg.full <- rbind(left_join(new2, st.reg2, by="Serotype"), c("NEG", rep(NA,4)))
# #dens <- as.numeric(st.reg.full$coef.st.fixed.od.1.2)
# dens <- as.numeric(st.reg.full$coef.st.lag.length.1.2)
# mean.dens <- mean(dens, na.rm=T)
# sd.dens <- sd(dens, na.rm=T)
# dens <- (dens-mean.dens)/sd.dens
# preds2 <- left_join(new2, preds, by="Serotype")
# preds2 <- rbind(preds2,preds[37,])
# #cor(preds[,2:7],method = "spearman", use = "complete.obs")
# #use cfr, TCR
# cfr <- preds2$CFRHARBOE
# mean.cfr <- mean(cfr, na.rm=T)
# sd.cfr <- sd(cfr, na.rm=T)
# cfr <- (cfr-mean.cfr)/sd.cfr
# tcr <- preds2$TOTALCARBONREPEAT
# mean.tcr <- mean(tcr, na.rm=T)
# sd.tcr <- sd(tcr, na.rm=T)
# tcr <- (tcr-mean.tcr)/sd.tcr

# jdat13 <- list(n=n, t1=t1, t2=t2, N_carr=N_carr,  N_swab=N_swab, k=k, dens=dens, tcr=tcr,cfr=cfr)
jdat13 <- list(n=n, t1=t1, t2=t2, N_carr=N_carr,  N_swab=N_swab, k=k)

jmod13 <- jags.model(textConnection(jcode13), data=jdat13, n.chains=3, n.adapt=15000)
adapt(object = jmod13, n.iter = 500000)
update(jmod13,10000)


jpos13 <- coda.samples(jmod13, thin=5, c('beta0.st', 'beta1.st', 'beta2.st', 
                                         'beta0', 'beta1', 
                                         'prev', 'prev.ratio', 'log.prevratio', 'prev.ratio.fin', 'log.prev.ratio.fin',
                                         'exp.beta0.st', 'exp.beta1.st','exp.beta2.st',
                                         'sigma0','sigma1','sigma2', 'rho'), n.iter=150000)
#
# # jpos13 <- coda.samples(jmod13, thin=5, c('rho'), n.iter=150000)
sum.post13 <- summary(jpos13)
View(sum.post13$quantiles)
# -0.6830686 -0.2486812 0.04182209 0.3103973 0.6781814

# dic <- dic.samples(jmod13, n.iter=150000, thin = 5, type = "pD")
# Mean deviance:  922.1 
# penalty 71.55 
# Penalized deviance: 993.6 

prevRsamps <- as.data.frame(jpos13[[1]][,1158:1542])
dim(prevRsamps)
#[1] 30000   385
write.csv(prevRsamps, "jsamps.jewish.prevR.csv")
prevsamps <- as.data.frame(jpos13[[1]][,773:1157])
dim(prevsamps)
#[1] 30000   385
write.csv(prevsamps, "~/Desktop/pneumo/jsamps.jewish.prev.csv")


#plot(jpos13, col=c("blue", "purple", "green"))

#########################################
#MCMC Diagnostics
#########################################
test <- sort(effectiveSize(jpos13))
summary(effectiveSize(jpos13)[which(effectiveSize(jpos13)!=0)])   
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1110    5642   11726   18102   29992   91431 
which(effectiveSize(jpos13)==0) #beta0.st and beta1.st for NEG

gew <- geweke.diag(jpos13)[[1]][[1]] 
hist(geweke.diag(jpos13)[[1]][[1]])
# geweke.plot(jpos13)[[1]][[1]] 

par(mfrow=c(3,2))
autocorr.plot(as.mcmc(jpos13[1]))
gelman <- (gelman.diag(jpos13, multivariate = FALSE))
#all are 1 

#############################################################################################################
#Inference
#############################################################################################################
par(mfrow=c(1,1))

#ds <- read.csv("~/Desktop/pneumo/multinomial3.csv")
#alpha0: 2-56, beta0: 59-113, beta1: 115-169, beta2:170-224, dens: 225-279, exp(beta0): 280-334, exp(beta1): 335-389, exp(beta2): 390-444, log.prev.ratio.fin 445:499, 
#logprevR 500:884, prev 885:1269, prevR 1270:1654, prev.ratio.fin 1655:1709

sts <- (as.character(ds1[,1]))

order.b0 <- sts[order(sum.post13$quantiles[,3][2:56], decreasing = TRUE)]
order.b1 <- sts[order(sum.post13$quantiles[,3][58:112], decreasing = TRUE)]
order.b2 <- sts[order(sum.post13$quantiles[,3][113:167], decreasing = TRUE)]
#testds <- read.csv("~/Desktop/pneumo/postFINALstatsNEGwVTs.csv")
order.prevratio <- sts[order(sum.post13$quantiles[,3][1543:1597], decreasing = TRUE)]
order.logprevratio <- sts[order(sum.post13$quantiles[,3][333:387], decreasing = TRUE)]


# mean(sum.post13$statistics[,1][1158:1542])
# #[1]  1.633674
# median(sum.post13$statistics[,1][1158:1542])
# # [1] 1.384523

mean(sum.post13$quantiles[,3][1158:1542])
#[1]  1.428722
median(sum.post13$quantiles[,3][1158:1542])
# [1] 1.303473

#Auranen: 1.47859922178988

#beta0s
par(mfrow=c(1,1))
caterplot(jpos13, parms = "beta0.st", labels = order.b0, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), xlab="Log(Relative Risk Ratio)", cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", beta [0][i])))
mtext(expression(bold('Log(Relative Risk Ratio at time 0)')), side = 1, line = 3)
mtext(expression(bold('Serotype')), side =2, line = 3)
#beta1s
caterplot(jpos13, parms = "beta1.st", labels = order.b1, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", beta [1][i])))
mtext(expression(bold('Log(Relative Risk Ratio per year)')), side = 1, line = 3)
mtext(expression(bold('Serotype')), side =2, line = 3)
#beta2s
caterplot(jpos13, parms = "beta2.st", labels = order.b2, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", beta [2][i])))
mtext(expression(bold('Log(Relative Risk Ratio per year)')), side = 1, line = 3)
mtext(expression(bold('Serotype')), side =2, line = 3)
#caterpoints(runif(5, 10, 20), pch="x", col="red")


caterplot(jpos13, parms = "exp.beta0.st", labels = order.b0, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", exp(beta [0][j]))))

caterplot(jpos13, parms = "exp.beta1.st", labels = order.b1, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", exp(beta [1][j]))))

caterplot(jpos13, parms = "exp.beta2.st", labels = order.b2, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), bty="n", lwd = c(1,1), cex.labels = 1, col="grey40")
title(expression(paste("Model Estimates of ", exp(beta [2][j]))))



caterplot(jpos13, parms = "prev.ratio.fin", labels = order.prevratio, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), lwd = c(1,1), bty='n', cex.labels = 1, col="grey40")
title(expression(paste("Model serotype prevalence ratios (", italic("time") [6] / italic("time") [0], ")")), adj=0)
caterpoints(runif(55, 1, 1), type="l", lty=2, col="black", lwd=2)
caterpoints(runif(55, 1.479, 1.479), type="l", lty=2, col="red", lwd=2)
#caterpoints(runif(55, 1.364, 1.364), type="l", lty=2, col="green4", lwd=2)
#caterplot(jpos13, parms = "prev.ratio", labels = order.prevratio, style = c("plain"), quantiles = list(outer=c(0.025,0.975)), lwd = c(1,1), add=TRUE)
legend(x="bottomright",legend=c('Prevalence ratio: 1','Nurhonen constant ratio: 1.48'),lwd=2, lty=2, col=c('black','red'),bty='n', cex=.6)
#caterplot(jpos7, parms = c("beta1", "alpha0", "alpha1", 'alpha2', 'alpha3'),style = c("plain"))

caterplot(jpos13, parms = "log.prev.ratio.fin", labels = order.logprevratio, style = c("plain"), quantiles = list(inner=c(0.025,0.975)), lwd = c(1,1), bty='n', cex.labels = 1, col="grey40")
#title(expression(paste("Figure 2. Model serotype log(prevalence ratios) (", italic("time") [6] / italic("time") [0], ")")), adj=0)
title("Jewish children")
caterpoints(runif(55, 0, 0), type="l", lty=2, col="gray60", lwd=3)
caterpoints(runif(55, log(1.479), log(1.479)), type="l", lty=4, col="black", lwd=3)
#caterpoints(runif(55, log(1.153), log(1.153)), type="l", lty=2, col="green4", lwd=2)
legend(x=.65, y=4.5,legend=c('log(prev. ratio): 0','Constant log(prev. ratio): 0.39'),lwd=2, lty=c(2,4), col=c('gray60','black'),bty = "n")
mtext(expression(bold('Log(prevalence ratio) over study period')), side = 1, line = 3)
mtext(expression(bold('Serotype')), side =2, line = 3)

write.csv(sum.post13$statistics, file = "~/Desktop/pneumo/postFINALstatsNEGwVTs_jewish_cat.time.csv")
write.csv(sum.post13$quantiles, file = "~/Desktop/pneumo/postFINALquantsNEGwVTs_jewish_cat.time.csv")

# rawprev <- read.csv(file = "~/Desktop/pneumo/raw prevalences.csv", header = FALSE)
# modprev <- sum.post13$quantiles[,3][336:665]
# plot(rawprev, modprev)

# prevs <- data.frame(rawprev, modprev)
# plot(prevs, xlim=c(0, .07), ylim=c(0., .07), xlab="Observed Prevalence", ylab="Model Prevalence")
# abline(a = 0, b=1, lty=2)
# title("Model Prevalences vs. Observed (Raw) Prevalences")

# poststats <- sum.post13$quantiles[,1]
# y6prevalences <- poststats[seq(229, length(poststats), 6)]
# y6prev <- as.numeric(as.character(unlist(y6prevalences)))


prevmatx <- matrix(sum.post13$quantiles[773:1157,3], ncol = 7, byrow = TRUE)
write.csv(prevmatx, file = "~/Desktop/pneumo/prevalences y0to6 final NEGwVTs jewish cat time.csv")

# library(ggplot2)
# RRR.b2 <- sum.post13$quantiles[390:444,3][-55]
# lRRR.b2 <- sum.post13$quantiles[170:224,3][-55]
# dens.mod <- sum.post13$quantiles[225:279,3][-55]
# plot(dens.mod, RRR.b2, xlab="Relative Density",ylab="Relative Risk Ratios")
# plot(dens.mod, lRRR.b2, xlab="Relative Density",ylab="Relative Risk Ratios")
# #plot(dens[-55], RRR.b2[-55])
# 
# mod <- lm(dens.mod ~RRR.b2) 
# summary(mod)
# 
# plot(RRR.b2 ~ dens.mod, col="seagreen4", pch=16, cex=1.5, bty="n", ylim=c(0,2), xlab="Relative Density", ylab="Relative Risk Ratios")

#a3 <- 0.07775035

#text(RRR.b2 ~ dens.mod, labels=ds1, cex=0.5, font=2)

#title("Relative Density vs. Relative Risk Ratio")
# title("Jewish children")

postquants <- as.data.frame(sum.post13$quantiles[,c(1,3,5)])
names(postquants) <- c("lCI","med.est","uCI")
postquants$CrI <- with(postquants, sprintf("%.2f, %.2f",lCI,uCI))
write.csv(postquants,"~/Desktop/pneumo/postquants final Jewish w CrI.csv")


