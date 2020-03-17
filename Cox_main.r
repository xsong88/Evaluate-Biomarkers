##############################################################
# Xiao Song: 03/17/2020                                      #                                 
#                                                            #
# In the code, the standard biomarker is denoted by W        #
#                                                            #
##############################################################

dir.name <- ""
data.file <- paste(dir.name, "ECOG_simul.dat", sep="")

isForEach <- TRUE
OverAllSurv <- FALSE

source(paste(dir.name,"Cox_functions.r",sep=""))


library(doParallel)
library(foreach)
library(doRNG)
library(parallel)
library(survival)
library(cubature)


isCluster <- TRUE


n.cluster <- 12

n.boot <- 500
n.boot.interval <- 500
n.SIMEX <- 500

n.method <- 3
n.method.beta <- 3
method.name <- c("ideal", "naive", "CondScore")

n.cov <- 3
iter.max <- 100
L.c <- 4


set.seed(342, kind = "L'Ecuyer-CMRG")
   

tstart <- proc.time()

beta.all <- matrix(NA,n.cov,n.method)
se.all <- matrix(NA,n.cov,n.method)
se.all.bs <- matrix(NA,n.cov,n.method)
se.all.mad <- matrix(NA,n.cov,n.method)

T.ResMean.all <- rep(NA,n.method)
se.all.T <- rep(NA,n.method)
mad.all.T <- rep(NA,n.method)


Rate.1 <- rep(NA,n.method)


dat <- read.table(data.file, head=TRUE)

n = length(dat$id)
W.i <- cbind(dat$W.1,dat$W.2)
dat$W <- apply(W.i,1,mean,na.rm=TRUE)    # mean of W.i

rss = rep(NA,n)        # residual sum of squares
for(i in 1:n)
{
  if(dat$m.i[i]>0)
  {
    rss[i]=sum((W.i[i,]-dat$W[i])^2)
  }
}

dat$rss <- rss

dat.rep <- NULL

# sort data by descending V 
dat <- dat[sort.list(-dat$V),]

T.ResMean.A.0 <- GetRestrictedLifeTime.Nonparm.A(dat,0,rep(1,n))
T.ResMean.A.1 <- GetRestrictedLifeTime.Nonparm.A(dat,1,rep(1,n))


fit <- GetBootStrapSE.T.A(dat,0)
T.ResMean.A.se.0 <- fit$T.ResMean.se
T.ResMean.A.mad.0 <- fit$T.ResMean.mad

fit <- GetBootStrapSE.T.A(dat,1)
T.ResMean.A.se.1 <- fit$T.ResMean.se
T.ResMean.A.mad.1 <- fit$T.ResMean.mad


sigma2 <- 0
nm <- 2
beta.N.mean <- rep(0,n.cov)
beta.ini <- rep(0,n.cov)

zph.p <- matrix(NA,n.method,4)

while(nm<=n.method)
{
  fit <- GetEstimate.any(method.name[nm],beta.ini,beta.N.mean,dat,dat.rep,Z=rep(1,n),TestCox=TRUE)
  
  beta.all[,nm] <- fit$beta.est
  se.all[,nm] <- fit$beta.se
  T.ResMean.all[nm] <- fit$T.ResMean
  Rate.1[nm] <- fit$Rate.1
  zph.p[nm,] <- fit$zph.p
  
  nm <- nm+1
  
  if(nm==2)
  {
    beta.ini <- beta.all[,2]
    beta.N.mean <- beta.all[,2]
  }
  
  sigma2 <- fit$sigma2
}


fit <- GetBootStrapSE.norep.all(beta.all[,],dat,dat.rep,sigma2)

se.all.bs <- fit$beta.se
se.all.mad <- fit$beta.mad
se.all.T <- fit$T.ResMean.se
mad.all.T <- fit$T.ResMean.mad


naive.index <- 2

T.ResMean.diff <- T.ResMean.all[-naive.index]-T.ResMean.all[naive.index]
T.ResMean.diff.se <- fit$T.ResMean.diff.se
T.ResMean.diff.mad <- fit$T.ResMean.diff.mad


Rate.1.se <- fit$Rate.1.se
Rate.1.mad <- fit$Rate.1.mad

Rate.1.diff <- Rate.1[-naive.index]-Rate.1[naive.index]
Rate.diff.se <- fit$Rate.diff.se
Rate.diff.mad <- fit$Rate.diff.mad


censor.rate <- sum(dat$delta)



cat("********************************************************************************************\n")
cat("********************************************************************************************\n")


cat("\n")


{
  cat("++++ T.ResMean: all assigned to treatment 0 ++++\n")
  out.list <- cbind(T.ResMean.A.0,T.ResMean.A.mad.0)
  colnames(out.list) <- c("est","se.mad")
  print(out.list)
  cat("\n")
  
  cat("++++ T.ResMean: all assigned to treatment 1 ++++\n")
  out.list <- cbind(T.ResMean.A.1,T.ResMean.A.mad.1)
  colnames(out.list) <- c("est","se.mad")
  print(out.list)
  cat("\n")
  
  cat("++++ Estimate of beta ++++\n")
  out.list <- cbind(beta.all[,-1],se.all.mad[,-1])
  colnames(out.list) <- c("naive.est", "condscore.est", "naive.se", "condscore.se")
  rownames(out.list) <- c("beta.1","beta.2","beta.3")
  print(out.list)
  cat("\n")
  
  cat("++++ T.ResMean: optimal ++++\n")
  out.list <- cbind(T.ResMean.all,mad.all.T)[-1,]
  out.list <- rbind(out.list,c(T.ResMean.diff[-1],T.ResMean.diff.mad))
  colnames(out.list) <- c("est","se.mad")
  rownames(out.list) <- c("standard marker","new marker","diff")
  print(out.list)
  cat("\n")
  
  cat("++++ Rate to treatment 1++++\n")
  out.list <- cbind(Rate.1[-1])
  colnames(out.list) <- c("estimated rate")
  rownames(out.list) <- c("standard marker","new marker")
  print(out.list)
  cat("\n")
  
  cat("++++ difference in assignment rate to treatment 1: new-standard\n")
  out.list <- cbind(Rate.1.diff[-1],Rate.diff.mad)
  colnames(out.list) <- c("est","se.mad")
  print(out.list)
  cat("\n")
  
  
  cat("++++ estimated error variance ++++\n")
  print(sigma2)
  
  
}


cat("\n\n")


cat("++++ censoring rate ++++\n")
censor.rate <- 1-mean(dat$delta)
print(censor.rate)
cat("\n")

tend <- proc.time()

cat("++++ processing time: ++++ \n")
print(tend-tstart)


cat("\n\n\n")