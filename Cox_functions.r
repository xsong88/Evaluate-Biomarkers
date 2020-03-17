##############################################################
# Xiao Song: 03/17/2020                                      #                                 
#                                                            #
# Functions used in Cox_main.r                               #
#                                                            #
##############################################################

GetSigma2.replicate <- function(dat,Z)
{  
   sigma2 <- sum(Z*dat$rss,na.rm=TRUE)/sum(Z*(dat$m.i-1))
   return(sigma2)
}


# type=1, ideal;
# type=2, naive;
# type=3, conditional score
# Z: weight
GetEstimate.any <- function(type.name,beta.ini,beta.N.mean,dat,dat.rep, Z=1, TestCox=FALSE)
{
    beta.est <- rep(NA,n.cov)
    beta.se <- rep(NA,n.cov)
    T.ResMean <- NA
    T.ResMean.nonlinear <- NA
    T.ResMean.p <- NA
    T.ResMean.p.2 <- NA
    Rate.1 <- NA
    D.miss <- rep(NA,6)
    sigma2 <- NA
    zph.p <- rep(NA,4)
  
    death.time <- sort(dat$V[dat$delta==1])
    
    dat$id <- seq_along(dat$V)
    dat$Z <- Z/dat$m.i
    data.temp <- dat[rep(row.names(dat),dat$m.i),]
    data.temp$W.s <- data.temp[cbind(row.names(data.temp),paste("W.",as.character(as.integer(sapply(strsplit(paste(row.names(data.temp),".0"), "\\."), "[", 2))+1),sep=""))]
    dat$Z <- Z
    

   
    Lambda0.est <- NA
    T.ResMean <- NA
    
    
   
    if(type.name %in% c("ideal", "naive"))
    {
        var.name <- ifelse(type.name=="ideal","X","W")
        formula.string <- paste("Surv(V, delta) ~ ", var.name,"*A",sep="")
        fit <- coxph(as.formula(formula.string),data=dat,weights=Z)
        zph <- cox.zph(fit)
        zph.p <- zph$table[,3]
        
        beta.est <- fit$coef
        beta.se <- sqrt(diag(fit$var))
        
        T.ResMean <- GetRestrictedLifeTime.Nonparm.optimal(beta.est,dat,type.name,Z=Z)
        Rate.1 <- GetAssignmentRate.mean(beta.est,dat,Z=Z)
      
    }    
    else if(type.name=="CondScore")
    {
        sigma2 <- GetSigma2.replicate(dat,Z)
        fit <- GetEstimate.cox(beta.ini,dat,dat.rep,sigma2,GetScoreInformation.X2,GetScore.X2.square,GetLambda0,Z=Z)
        
        beta.est <- fit$beta
        beta.se <- fit$beta.se
        Lambda0.est <- fit$Lambda0
        
        death.time <- death.time[death.time<=L.c]
        Lambda0.est <- Lambda0.est[1:(length(death.time)+1)]

        if(!is.na(beta.est[1]))
        {
            TD <-  GetRestrictedLifeTime.Nonparm.SIMEX(beta.est,beta.N.mean,dat,data.temp,sigma2,Z=Z,TestCox=TestCox)
            T.ResMean <- TD$T.mean.square 
            T.ResMean.nonlinear <- TD$T.mean.nonlinear 
            Rate.1 <- TD$Rate.1
            zph.p <- TD$zph.p
        }
    }


    
    if(TestCox)
    {
      return.list <- list(beta.est=beta.est,beta.se=beta.se,
                          T.ResMean=T.ResMean,
                          sigma2=sigma2,
                          Rate.1=Rate.1,
                          zph.p=zph.p)
    }
    else
    {
      return.list <- list(beta.est=beta.est,beta.se=beta.se,
                          T.ResMean=T.ResMean,
                          Rate.1=Rate.1)
    }

    return(return.list)
}


GetEstimate.cox <- function(beta.ini,data,data.rep,error.var,GetScoreInformation.X2,GetScore.X2.square,GetLambda0,Z,bIncInteraction=TRUE)
{
    max.iter <- 100
    beta0 <- beta.ini
    bSuccess <- FALSE
    iter <- 0
    while(!bSuccess && iter < max.iter)
    {
        scoreinfor <- GetScoreInformation.X2(beta0,data,error.var/data$m.i,Z=Z,bIncInteraction=bIncInteraction)
        scoreinfor.inv <- try(solve(scoreinfor$information),TRUE)
        if(is.na(scoreinfor.inv[1]) || is.nan(scoreinfor.inv[1]) || is.character(scoreinfor.inv[1]) || is.na(scoreinfor$score[1]))
        {
            beta <- rep(NA,n.cov)
            break
        }
        else
        {
            beta <- beta0-as.vector(scoreinfor.inv%*%scoreinfor$score)
            bSuccess <- all(abs(beta-beta0)<=(1e-5)*abs(beta))
            beta0 <- beta
        }
        iter <- iter+1
    }
    if(!bSuccess || any(beta.ini*beta<0) || max(abs(beta-beta.ini))>10*max(abs(beta.ini)))
    {
    
        fit <- optim(beta.ini,GetScore.X2.square,method="Nelder-Mead",dat=data,sigma2=error.var/data$m.i,Z=Z)
        

        if(fit$convergence==0)
        {
            beta <- fit$par
            seLambda0 <- GetLambda0(beta,data,data.rep,error.var/data$m.i,Z=Z)
            Lambda0 <- seLambda0$Lambda0
            beta.se <- rep(NA,n.cov)
        }else
        {
            beta <- rep(NA,n.cov)
            beta.se <- rep(NA,n.cov)
            Lambda0 <- NA
        }
    }
    else
    {
        seLambda0 <- GetLambda0(beta,data,data.rep,error.var/data$m.i,Z=Z,bIncInteraction=bIncInteraction)
        beta.se <- seLambda0$beta.se
        Lambda0 <- seLambda0$Lambda0
    }

    return(list(beta=beta,beta.se=beta.se,Lambda0=Lambda0))
}



# conditional score
GetScoreInformation.X2 <- function(beta,dat,sigma2,Z,bIncInteraction=TRUE)
{
  n <- nrow(dat)
  if(bIncInteraction)
  { 
    W.ext <- as.matrix(cbind(dat$W,dat$A,dat$W*dat$A))
    pred.der.X <- cbind(1,0,dat$A)
    n.cov.t <- n.cov
  }
  else
  {
    W.ext <- as.matrix(cbind(dat$W,dat$A))
    pred.der.X <- cbind(rep(1,n),0)
    n.cov.t <- n.cov-1
  }
  
  corr.term <- rep(0,n)
  Omega.beta <- matrix(rep(0,n*n.cov.t),nrow=n)         # The first column is zero, reserved for the intercept
  Omega.mat <- matrix(rep(0,n*n.cov.t*n.cov.t),nrow=n)    # The first row and first column is zero
  for(i in 1:n)
  {
    Omega <- sigma2[i]*pred.der.X[i,]%*%t(pred.der.X[i,])
    Omega.beta[i,] <- Omega%*%beta
    corr.term[i] <-t(beta)%*%Omega.beta[i,]
    Omega.mat[i,] <- c(Omega)
  }
  
  temp0 <- c(exp(W.ext%*%beta-corr.term/2)*Z)
  temp1 <- W.ext*temp0
  temp0.c <- c(exp(W.ext%*%beta-corr.term/2+corr.term*dat$delta)*Z)
  W.ext.c <- W.ext+Omega.beta*dat$delta
  temp1.c <- W.ext.c*temp0.c
  E0 <- cumsum(temp0)
  E1 <- apply(temp1,2,cumsum)
  E0 <- E0+(temp0.c-temp0)*dat$delta
  E1 <- E1+(temp1.c-temp1)*dat$delta
  temp1.der <- matrix(0,n,n.cov.t*n.cov.t)
  temp1.der.c <- matrix(0,n,n.cov.t*n.cov.t)
  E1.2 <- matrix(0,n,n.cov.t*n.cov.t)
  temp0.der.beta <- (W.ext-Omega.beta)*temp0
  temp0.der.beta.c <- (W.ext+Omega.beta*(2*dat$delta-1))*temp0.c
  E0.der <- apply(temp0.der.beta,2,cumsum)+(temp0.der.beta.c-temp0.der.beta)*dat$delta
  for(i in 1:n)
  {
    temp1.der[i,] <- as.vector(temp1[i,]%*%t(W.ext[i,]-Omega.beta[i,]))
    temp1.der.c[i,] <- as.vector(temp1.c[i,]%*%t(W.ext[i,]+Omega.beta[i,]*(2*dat$delta[i]-1)))
    E1.2[i,] <- as.vector(E1[i,]%*%t(E0.der[i,]))
  }
  temp1.der.c <- temp1.der.c+Omega.mat*dat$delta*temp0.c
  E2 <- apply(temp1.der.c,2,cumsum)
  E2 <- E2+(temp1.der.c-temp1.der)*dat$delta
  score <- colMeans((W.ext.c-E1/E0)*dat$delta*Z)
  E0 <- as.vector(E0)
  information <- matrix(colMeans((Omega.mat-(E2/E0-E1.2/(E0*E0)))*dat$delta*Z),ncol=n.cov.t)
  return(list(score=score,information=information))
}




GetScore.X2.square <- function(beta,dat,sigma2,Z,bIncInteraction=TRUE)
{
  n <- nrow(dat)
  if(bIncInteraction)
  { 
    W.ext <- as.matrix(cbind(dat$W,dat$A,dat$W*dat$A))
    pred.der.X <- cbind(1,0,dat$A)
    n.cov.t <- n.cov
  }
  else
  {
    W.ext <- as.matrix(cbind(dat$W,dat$A))
    pred.der.X <- cbind(rep(1,n),0)
    n.cov.t <- n.cov-1
  }
  
  corr.term <- rep(0,n)
  Omega.beta <- matrix(rep(0,n*n.cov.t),nrow=n)         # The first column is zero, reserved for the intercept
  Omega.mat <- matrix(rep(0,n*n.cov.t*n.cov.t),nrow=n)    # The first row and first column is zero
  for(i in 1:n)
  {
    Omega <- sigma2[i]*pred.der.X[i,]%*%t(pred.der.X[i,])
    Omega.beta[i,] <- Omega%*%beta
    corr.term[i] <-t(beta)%*%Omega.beta[i,]
    Omega.mat[i,] <- c(Omega)
  }
  
  temp0 <- c(exp(W.ext%*%beta-corr.term/2)*Z)
  temp1 <- W.ext*temp0
  temp0.c <- c(exp(W.ext%*%beta-corr.term/2+corr.term*dat$delta)*Z)
  W.ext.c <- W.ext+Omega.beta*dat$delta
  temp1.c <- W.ext.c*temp0.c
  E0 <- cumsum(temp0)
  E1 <- apply(temp1,2,cumsum)
  E0 <- E0+(temp0.c-temp0)*dat$delta
  E1 <- E1+(temp1.c-temp1)*dat$delta
  score <- colMeans((W.ext.c-E1/E0)*dat$delta*Z)
  obj <- sum(score*score)
  return(obj)
}




GetLambda0 <- function(beta,dat,dat.rep,sigma2,Z,bIncInteraction=TRUE)
{
    n <- nrow(dat)
    n.rep <- nrow(dat.rep)

    if(bIncInteraction)
    { 
      W.ext <- as.matrix(cbind(dat$W,dat$A,dat$W*dat$A))
      pred.der.X <- cbind(1,0,dat$A)
      n.cov.t <- n.cov
    }
    else
    {
      W.ext <- as.matrix(cbind(dat$W,dat$A))
      pred.der.X <- cbind(rep(1,n),0)
      n.cov.t <- n.cov-1
    }
    corr.term <- rep(0,n)
    Omega.beta <- matrix(rep(0,n*n.cov.t),nrow=n)         # The first column is zero, reserved for the intercept
    Omega.mat <- matrix(rep(0,n*n.cov.t*n.cov.t),nrow=n)    # The first row and first column is zero
    for(i in 1:n)
    {
        Omega <- sigma2[i]*pred.der.X[i,]%*%t(pred.der.X[i,])
        Omega.beta[i,] <- Omega%*%beta
        corr.term[i] <-t(beta)%*%Omega.beta[i,]
        Omega.mat[i,] <- c(Omega)
    }


    temp0 <- c(exp(W.ext%*%beta-corr.term/2)*Z)
    temp0.c <- c(exp(W.ext%*%beta-corr.term/2+corr.term*dat$delta)*Z)
    E0.corr <- cumsum(temp0)
    E0 <- E0.corr+(temp0.c-temp0)*dat$delta
    E0 <- as.vector(E0)
    
    
    beta.se <- rep(NA,n.cov)
    Lambda0 <- NA
   
    Lambda0 <- c(0,cumsum(rev((Z*dat$delta/E0)[dat$delta==1])))
    
    
    return(list(beta.se=beta.se,Lambda0=Lambda0))
}





GetCumHaz.Nonparm <- function(dat,D,A,Z)
{
    temp0 <- as.integer((dat$A==A)&(D==A))*Z
    E0 <- rev(cumsum(temp0))
    temp1 <- rev(temp0*dat$delta)
    Lambda <- cumsum(temp1/E0)[rev(dat$delta==1&dat$A==A&D==A)]

    return(Lambda)
}


GetAssignmentRate.mean <- function(beta,dat,Z)
{
  
  D <- mean(as.integer(Z*(beta[2]+beta[3]*dat$W) <0 ))
  
  return(D)
}


GetRestrictedLifeTime.Nonparm.optimal <- function(beta,dat,type.name,Z,bPrint=FALSE)
{
    n <- nrow(dat)
    T.mean <- NA

    if(type.name=="ideal")  # ideal
    {
        D <- as.integer((beta[2]+beta[3]*dat$X) <0 )
    }
    else
    {
        D <- as.integer((beta[2]+beta[3]*dat$W) <0 )

    }
    
    
    v.1 <- dat$V[dat$A==1 & D==1]
    v.0 <- dat$V[dat$A==0 & D==0]
    
    b.continue <- FALSE
    if(length(v.1)>0 & length(v.0)>0)
    {
        v.1 <- v.1[1]
        v.0 <- v.0[1]
        
        if((v.1>L.c) & (v.0>L.c))
        {
          b.continue <- TRUE
        }
        
    }
    
    if(b.continue)
    {
      p.1 <- sum(Z*D)/sum(Z)
      death.time <- rev(dat$V[dat$delta==1 & dat$A==1 & D==1])
      death.time <- death.time[death.time<=L.c]  
      if(length(death.time)>0)
      {
        Lambda <- c(0,GetCumHaz.Nonparm(dat,D,A=1,Z=Z)[1:length(death.time)])
        exp.Lambda <- exp(-Lambda)
        death.time.diff <- diff(c(0,death.time,L.c))
    
        exp.Lambda.ave <- (exp.Lambda+c(exp.Lambda[-1],exp.Lambda[length(exp.Lambda)]))/2
        T.mean <- p.1*sum(exp.Lambda.ave*death.time.diff)
      }
      else 
      {
        T.mean <- p.1*L.c
      }

      death.time <- rev(dat$V[dat$delta==1 & dat$A==0 & D==0])
      death.time <- death.time[death.time<=L.c]
      if(length(death.time)>0)
      {
        Lambda <- c(0,GetCumHaz.Nonparm(dat,D,A=0,Z=Z)[1:length(death.time)])
        exp.Lambda <- exp(-Lambda)
        death.time.diff <- diff(c(0,death.time,L.c))
    
        exp.Lambda.ave <- (exp.Lambda+c(exp.Lambda[-1],exp.Lambda[length(exp.Lambda)]))/2
        T.mean <- T.mean+(1-p.1)*sum(exp.Lambda.ave*death.time.diff)
      }
      else
      {
          T.mean <- T.mean+(1-p.1)*L.c
     
      }
    }
    

    return(T.mean)
}


GetRestrictedLifeTime.Nonparm.A <- function(dat,A,Z)
{
    death.time <- rev(dat$V[dat$delta==1 & dat$A==A])
    death.time <- death.time[death.time<=L.c]
    Lambda <- c(0,GetCumHaz.Nonparm(dat,A,A,Z=Z)[1:length(death.time)])
    exp.Lambda <- exp(-Lambda)
    death.time.diff <- diff(c(0,death.time,L.c))

    exp.Lambda.ave <- (exp.Lambda+c(exp.Lambda[-1],exp.Lambda[length(exp.Lambda)]))/2
    T.mean <- sum(exp.Lambda.ave*death.time.diff)

    return(T.mean)
}

logit <- function(p)
{
  x <- 1/(1+exp(-p))
  return(x)
}

GetRestrictedLifeTime.Nonparm.SIMEX <- function(beta,beta.N.mean,dat,dat.temp,sigma2,Z,TestCox=FALSE,extrotype=1)
{
    T.mean.square <- NA
    T.mean.nonlinear <- NA
    Rate.1 <- NA
    D.miss <- rep(NA,6)
    
    zph.p <- rep(NA,4)
    
    n <- nrow(dat)
    dat.old <- dat
    dat.temp.old <- dat.temp
    ksi <- seq(0,2,0.2)
    T.mean.kb <- matrix(0,n.SIMEX,length(ksi)-1)  
    Rate.1.kb <- matrix(0,n.SIMEX,length(ksi)-1)  
    
    if(TestCox)
    {
      zph.p.kb <- array(0,c(4,n.SIMEX,length(ksi)-1))
    }
    
    
    T.mean.naive <- GetRestrictedLifeTime.Nonparm.optimal(beta,dat,type="naive",Z=Z)
   
    Rate.1.naive <- GetAssignmentRate.mean(beta,dat,Z=Z)
    
    if(TestCox)
    {
        fit <- coxph(Surv(V, delta) ~ W*A, data=dat)
        zph <- cox.zph(fit)
        zph.naive <- zph$table[,3]
    }
     
    for(k in 2:length(ksi))
    {       
        for(b in 1:n.SIMEX)
        {
            dat$W <- dat$W+rnorm(n,mean=0,sd=sqrt(ksi[k]*sigma2/dat$m.i))
        
            T.mean.kb[b,k-1] <- GetRestrictedLifeTime.Nonparm.optimal(beta,dat,type="naive",Z=Z)
        
            
            Rate.1.kb[b,k-1] <- GetAssignmentRate.mean(beta,dat,Z=Z)
           
            if(TestCox)
            {
                fit <- coxph(Surv(V, delta) ~ W*A, data=dat)
                zph <- cox.zph(fit)
                zph.p.kb[,b,k-1] <- zph$table[,3]
            }
            dat <- dat.old   
            dat.temp <- dat.temp.old 
        }      
    }
    
    T.mean.square <- Extrapolate(T.mean.naive,T.mean.kb,ksi,"T.mean.square",type=1)
    T.mean.nonlinear <- Extrapolate(T.mean.naive,T.mean.kb,ksi,"T.miss.nonlinear",type=2)
    Rate.1 <- Extrapolate(logit(Rate.1.naive),logit(Rate.1.kb),ksi,"Rate.1",type=1)
    Rate.1 <- log(Rate.1/(1-Rate.1))
    
    if(TestCox)
    {
      for(j in 1:4)
      {
        zph.p[j] <- Extrapolate(zph.naive[j],zph.p.kb[j,,],ksi,"zph.p",type=1)
      }
    }
    
    
    return(list(T.mean.square=T.mean.square,Rate.1=Rate.1,D.miss=D.miss,zph.p=zph.p))
}


#type=1, least square, type=2, nonlinear
Extrapolate <- function(T.mean.naive,T.mean.kb,ksi,var.name,type=1)
{
  T.mean <- NA
  
  if(!is.na(T.mean.naive))
  {
    T.mean.b <- c(T.mean.naive,colMeans(T.mean.kb,na.rm=TRUE))
    
    ksi.2 <- ksi^2
    ksi.3 <- ksi^3
    
    
    fit <- try(lm(T.mean.b~ksi+ksi.2))
    if(!is.character(fit))
    {
      y.fit <- fit$fitted.values
      
      T.mean <- predict(fit,list(ksi=-1,ksi.2=1)) 
      
      if(type>1)
      {
        # nonlinear least square
        ksi.s <- c(0,1,2)
        T.mean.fitted <- predict(fit,list(ksi=ksi.s,ksi.2=ksi.s^2))
        d <- diff(T.mean.fitted)
        gamma3 <- (d[2]*ksi.s[3]*(ksi.s[2]-ksi.s[1])-ksi.s[1]*d[1]*(ksi.s[3]-ksi.s[2]))/(d[1]*(ksi.s[3]-ksi.s[2])-d[2]*(ksi.s[2]-ksi.s[1]))
        gamma2 <- d[2]*(gamma3+ksi.s[2])*(gamma3+ksi.s[3])/(ksi.s[3]-ksi.s[2])
        gamma1 <- T.mean.fitted[1]-gamma2/(gamma3+ksi.s[1])
        
        fit.nls <-try(nls(T.mean.b~gamma1+gamma2/(gamma3+ksi),start=list(gamma1=gamma1,gamma2=gamma2,gamma3=gamma3)),silent=TRUE)
        
        if(!is.character(fit.nls))
        {
          T.mean.temp <- predict(fit.nls,list(ksi=-1))
          if(sign(T.mean.temp-T.mean.naive)==sign(T.mean.naive-T.mean.b[2]))
          {
             T.mean <- T.mean.temp
          }
        }
      }
      
    }
    else
    {
      cat("error in SIMEX for ",var.name,"\n")
      
    }
  }
  
  return(T.mean)
}



GetBootStrapSE.norep.all <- function(beta,dat,dat.rep,error.var)
{
    n <- nrow(dat)
    n.rep <- nrow(dat.rep)

    #D.miss.bs <- matrix(NA,2,n.boot)

    getDoParWorkers()
    mycluster = makeCluster(n.cluster)
    clusterSetRNGStream(mycluster, 123)
    registerDoParallel(mycluster)
    getDoParWorkers()

    
    fit.all <- NULL

    bs.result <- foreach(i=1:n.boot.interval, .combine=rbind, .multicombine=TRUE,
                .export=c('GetEstimate.any','GetRestrictedLifeTime.Nonparm.optimal','GetCumHaz.Nonparm','L.c',
                'GetEstimate.cox',
                'GetScoreInformation.X2',
                'n.cov',
                'GetScore.X2.square',
                'GetLambda0',
                'GetRestrictedLifeTime.Nonparm.SIMEX',
                'n.SIMEX','n.method','method.name', 'beta.ini', 
                'coxph','Surv','basehaz','GetSigma2.replicate',
                'sigma2','cluster','cox.zph',
                 'GetAssignmentRate.mean',
               'Extrapolate','isForEach','logit')) %dorng% #%dopar%
    {
        
        Z <- (rbeta(n,sqrt(2)-1,1)-(1-sqrt(2)/2))*sqrt(2*(3+2*sqrt(2)))+1
        
        if(isForEach)
        {
          fit.all <- NULL
        }
        
        beta.ini <- rep(NA,n.cov)
        beta.N.mean <- rep(NA,n.cov)
      
        nt <- 2
        
        while(nt<=n.method)
        {    
        
            
            fit <- list(beta.est=rep(NA,n.cov),beta.se=rep(NA,n.cov),  T.ResMean=NA, sigma2=sigma2,  Rate.1=NA)
            if(nt>1)
            {
              fit <- GetEstimate.any(method.name[nt],beta.ini,beta.N.mean,dat,dat.rep,Z=Z)
            }
            
            nt <- nt+1
            
            if(nt==3)
            {
              beta.ini <- fit$beta.est
              beta.N.mean <- fit$beta.est
            }
            

            fit.all <- c(fit.all,unlist(fit))

            if(length(unlist(fit))!=16)
            {
              cat("i",length(unlist(fit)),"\n")
            }
        }
        
        
        return(fit.all)
        
    }  

    stopCluster(mycluster)

    beta.se <- matrix(NA,n.cov,n.method)
    beta.mad <- matrix(NA,n.cov,n.method)
    T.ResMean.se <- rep(NA,n.method)
    T.ResMean.mad <- rep(NA,n.method)
    beta.ci <- array(NA,c(2,n.cov,n.method))
    T.ci <- matrix(NA,2,n.method)
    Rate.1.se <- rep(NA,n.method)
    Rate.1.mad <- rep(NA,n.method)
    Rate.1.ci <- matrix(NA,2,n.method)
   
    
    i <- 2
    j <- 1
    
    T.ResMean.all <- NULL
    rate.1.all <- NULL
    while(i<=n.method)
    {
      if(isForEach)      
      {
        beta.bs <- bs.result[,8*(j-1)+(1:3)]
        T.ResMean.bs <- bs.result[,8*(j-1)+7] 
        rate.1.bs <- bs.result[,8*(j-1)+8]
      }  
      else
      {
        beta.bs <- bs.result[seq(i,n.boot.interval*n.method.beta,n.method.beta),1:3]
        T.ResMean.bs <- bs.result[seq(i,n.boot.interval*n.method.beta,n.method.beta),7]
        rate.1.bs <- bs.result[seq(i,n.boot.interval*n.method.beta,n.method.beta),8]
      }
      
      beta.bs <- matrix(unlist(beta.bs),nrow=n.boot.interval)
      T.ResMean.bs <- unlist(T.ResMean.bs)
      rate.1.bs <- matrix(unlist(rate.1.bs),nrow=n.boot.interval)
      
      # Compute SE
      beta.se[,i] <- apply(beta.bs,2,sd,na.rm=TRUE)
      beta.mad[,i] <- apply(beta.bs,2,mad,na.rm=TRUE)
      T.ResMean.se[i] <- sd(T.ResMean.bs,na.rm=TRUE)
      T.ResMean.mad[i] <- mad(T.ResMean.bs,na.rm=TRUE)
      
      Rate.1.se[i] <- sd(rate.1.bs,na.rm=TRUE)
      Rate.1.mad[i] <- mad(rate.1.bs,na.rm=TRUE)
      
      # Compute confidence interval
      beta.ci[,,i] <- apply(beta.bs,2,quantile,c(0.025,0.975),na.rm=TRUE)
      T.ci[,i] <- quantile(T.ResMean.bs,c(0.025,0.975),na.rm=TRUE)
      Rate.1.ci[,i] <- quantile(rate.1.bs,c(0.025,0.975),na.rm=TRUE)
      
      T.ResMean.all <- cbind(T.ResMean.all,T.ResMean.bs)
      rate.1.all <- cbind(rate.1.all,rate.1.bs)
      
      i <- i+1
      j <- j+1
    }
 
    naive.index <- 1
    T.ResMean.diff.mean.bs <- T.ResMean.all[,-naive.index]-T.ResMean.all[,naive.index]
    rate.1.diff.mean.bs <- rate.1.all[,-naive.index]-rate.1.all[,naive.index]
    
    T.ResMean.diff.se <- sd(T.ResMean.diff.mean.bs,na.rm=TRUE)
    T.ResMean.diff.mad <- mad(T.ResMean.diff.mean.bs,na.rm=TRUE)
    T.ResMean.diff.ci <- quantile(T.ResMean.diff.mean.bs,c(0.025,0.975),na.rm=TRUE)
    
    Rate.diff.se <- sd(rate.1.diff.mean.bs,na.rm=TRUE)
    Rate.diff.mad <- mad(rate.1.diff.mean.bs,na.rm=TRUE)
    Rate.diff.ci <- quantile(rate.1.diff.mean.bs,c(0.025,0.975),na.rm=TRUE)
    
    return(list(beta.se=beta.se,beta.mad=beta.mad,
                T.ResMean.se=T.ResMean.se, T.ResMean.mad=T.ResMean.mad,
                beta.ci=beta.ci,T.ci=T.ci,
                T.ResMean.diff.se=T.ResMean.diff.se,T.ResMean.diff.mad=T.ResMean.diff.mad,T.ResMean.diff.ci=T.ResMean.diff.ci,
                Rate.1.se=Rate.1.se, Rate.1.mad=Rate.1.mad,Rate.1.ci=Rate.1.ci,
                Rate.diff.se=Rate.diff.se, Rate.diff.mad=Rate.diff.mad, Rate.diff.ci=Rate.diff.ci
                ))
}



GetBootStrapSE.T.A <- function(dat,A)
{
    n <- nrow(dat)

    T.ResMean.bs <- rep(NA,n.boot)
    
    getDoParWorkers()
    mycluster = makeCluster(n.cluster)
    clusterSetRNGStream(mycluster, 123)
    registerDoParallel(mycluster)
    getDoParWorkers()
    
   
    bs.result <-  foreach(i=1:n.boot, .combine=rbind, .multicombine=TRUE,
                 .export=c('GetRestrictedLifeTime.Nonparm.A','GetCumHaz.Nonparm','L.c')) %dorng% #%dopar%
    {
        Z <- (rbeta(n,sqrt(2)-1,1)-(1-sqrt(2)/2))*sqrt(2*(3+2*sqrt(2)))+1
        index <- sample(1:n,n,replace=TRUE)
        dat.bs <- dat[index,]
        dat.bs <- dat.bs[sort.list(-dat.bs$V),]

        T.ResMean.i <- GetRestrictedLifeTime.Nonparm.A(dat,A,Z=Z)
        return(c(T.ResMean.i))
    }
    
    T.ResMean.bs <- bs.result[,1]
    stopCluster(mycluster)
    
    T.ResMean.bs <- unlist(T.ResMean.bs)
   
    # Compute SE
    T.ResMean.se <- sd(T.ResMean.bs,na.rm=TRUE)
    T.ResMean.mad <- mad(T.ResMean.bs,na.rm=TRUE)
    
    T.ci <- quantile(T.ResMean.bs,c(0.025,0.975),na.rm=TRUE)
   
  
    return(list(T.ResMean.se=T.ResMean.se,T.ResMean.mad=T.ResMean.mad,T.ci=T.ci))
 }
