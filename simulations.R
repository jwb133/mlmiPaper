#This R program performs simulations for the paper "Maximum likelihood multiple imputation: Imputation can work without posterior draws"
#by von Hippel and Bartlett
#It is intended to be called on a high performance cluster. If instead you want to run it on a single machine, you
#need to comment out the lines which are indicated in the file. A separate R program combines the resulting outputted
#datasets into a single file. An RMarkdown file then loads the results and tabulates them.

library(MASS)
library(mlmi)
library(bootImpute)
library(Matrix)

#set working directory
#setwd("C:...")

#Set-up for running on high performance cluster and setting random number seed. Delete this tabbed section
#if you want to run on a single machine
  slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  # coerce the value to an integer
  batch <- as.numeric(slurm_arrayid)
  
  #find seed for this run
  library(parallel)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(69012365) # something
  s <- .Random.seed
  for (i in 1:batch) {
    s <- nextRNGStream(s)
  }
  .GlobalEnv$.Random.seed <- s
  
  print(batch)

  
#Define functions that will perform analyses
yonx <- function(inputData,vcovtype) {
  fitmod <- lm(y~x, data=inputData)
  n <- dim(inputData)[1]
  if (vcovtype==1) {
    #just the two regression coefficients
    list(est=fitmod$coef, var=vcov(fitmod))
  } else if (vcovtype==2) {
    #now include the model residual variance
    myvcov <- as.matrix(bdiag(vcov(fitmod), 2*sigma(fitmod)^4/(fitmod$df.residual)))
    list(est=c(fitmod$coef, sigma(fitmod)^2), var=myvcov)
  } else {
    #now all 5 parameters of the bivariate normal model
    myvcov <- as.matrix(bdiag(vcov(fitmod), 2*sigma(fitmod)^4/(fitmod$df.residual),var(inputData$x)/(n-1), 2*var(inputData$x)^2/(n-1)))
    list(est=c(fitmod$coef, sigma(fitmod)^2, mean(inputData$x), var(inputData$x)), var=myvcov)
  }
}

yonxscore <- function(inputData, parm, vcovtype) {
  beta0 <- parm[1]
  beta1 <- parm[2]
  res <- inputData$y - beta0 - beta1*inputData$x
  if (vcovtype==1) {
    cbind(res, (res*inputData$x))
  } else if (vcovtype==2) {
    sigmasq <- parm[3]
    cbind(res/sigmasq, (res*inputData$x)/sigmasq, res^2/(2*sigmasq^2)-1/(2*sigmasq))
  } else {
    sigmasq <- parm[3]
    mux <- parm[4]
    sigmaxsq <- parm[5]
    xres <- inputData$x - mux
    cbind(res/sigmasq, (res*inputData$x)/sigmasq, res^2/(2*sigmasq^2)-1/(2*sigmasq),
          xres/sigmaxsq,xres^2/(2*sigmaxsq^2)-1/(2*sigmaxsq))
  }
}

xonyscore <- function(inputData, parm, vcovtype) {
  beta0 <- parm[1]
  beta1 <- parm[2]
  res <- inputData$x - beta0 - beta1*inputData$y
  if (vcovtype==1) {
    cbind(res, (res*inputData$y))
  } else if (vcovtype==2) {
    sigmasq <- parm[3]
    cbind(res/sigmasq, (res*inputData$y)/sigmasq, res^2/(2*sigmasq^2)-1/(2*sigmasq))
  } else {
    sigmasq <- parm[3]
    muy <- parm[4]
    sigmaysq <- parm[5]
    yres <- inputData$y - muy
    cbind(res/sigmasq, (res*inputData$y)/sigmasq, res^2/(2*sigmasq^2)-1/(2*sigmasq),
          yres/sigmaysq,yres^2/(2*sigmaysq^2)-1/(2*sigmaysq))
  }
}

yonxboot <- function(inputData) {
  fitmod <- lm(y~x, data=inputData)
  fitmod$coef[2]
}

xony <- function(inputData,vcovtype) {
  fitmod <- lm(x~y, data=inputData)
  n <- dim(inputData)[1]
  if (vcovtype==1) {
    #just the two regression coefficients
    list(est=fitmod$coef, var=vcov(fitmod))
  } else if (vcovtype==2) {
    #now include the model residual variance
    myvcov <- as.matrix(bdiag(vcov(fitmod), 2*sigma(fitmod)^4/(fitmod$df.residual)))
    list(est=c(fitmod$coef, sigma(fitmod)^2), var=myvcov)
  } else {
    #now all 5 parameters of the bivariate normal model
    myvcov <- as.matrix(bdiag(vcov(fitmod), 2*sigma(fitmod)^4/(fitmod$df.residual),var(inputData$y)/(n-1), 2*var(inputData$y)^2/(n-1)))
    list(est=c(fitmod$coef, sigma(fitmod)^2, mean(inputData$y), var(inputData$y)), var=myvcov)
  }
}

xonyboot <- function(inputData) {
  fitmod <- lm(x~y, data=inputData)
  fitmod$coef[2]
}

#set number of simulations. This is 10 here because when running on an HPC each independent instance runs
#10 simulations. If running on a single machine, increase this number to something sensible, like 1000
nSim <- 10
N <- 500

numMethods <- 10
missPropLevels <- c(25,50)
missMecLevels <- c("MCAR", "MAR")
MLevels <- c(10, 50, 200)
ests <- array(0, dim=c(2,length(missPropLevels),length(missMecLevels),length(MLevels),nSim,numMethods))
se <- array(0, dim=c(2,length(missPropLevels),length(missMecLevels),length(MLevels),nSim,numMethods))
ci <- array(0, dim=c(2,length(missPropLevels),length(missMecLevels),length(MLevels),nSim,numMethods,2))
df <- array(0, dim=c(2,length(missPropLevels),length(missMecLevels),length(MLevels),nSim,numMethods))

i <- 0
for (yonxOrxony in c(0,1)) {
  i <- i + 1
  j <- 0
  for (missProp in missPropLevels) {
    j <- j + 1
    k <- 0
    print(missProp)
    for (missMec in missMecLevels) {
      k <- k + 1
      l <- 0
      print(missMec)
      for (M in MLevels) {
        l <- l + 1
        p <- missProp/100
        print(M)
  
        res <- array(0, dim=c(nSim,4))
  
        for (sim in 1:nSim) {
          print(sim)
          #simulate data
          simData <- mvrnorm(n=N, c(0,0), Sigma=matrix(c(1,0.5,0.5,1), nrow=2))
          colnames(simData) <- c("x", "y")
          simData <- data.frame(simData)
  
          #make data missing
          if (missMec=="MCAR") {
            simData$y[(runif(N)<p)] <- NA
          } else {
            simData$y[(runif(N)<(2*p*pnorm(simData$x)))] <- NA
          }
  
          methodNum <- 1
          
          #PDMI
          imps <- mlmi::normUniImp(simData, y~x, M=M, pd=TRUE)
          
            #WB
            if (yonxOrxony==1) {
              fit <- mlmi::withinBetween(imps, yonx, vcovtype=1, dfComplete=rep(N-2,2))
            } else {
              fit <- mlmi::withinBetween(imps, xony, vcovtype=1, dfComplete=rep(N-2,2))
            }
            ests[i,j,k,l,sim,methodNum] <- fit$est[2]
            se[i,j,k,l,sim,methodNum] <- fit$var[2,2]^0.5
            df[i,j,k,l,sim,methodNum] <- fit$df[2]
            ci[i,j,k,l,sim,methodNum,] <- c(fit$est[2]-qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5, fit$est[2]+qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5)
            methodNum <- methodNum + 1
            
            #SB including residual variance
            if (yonxOrxony==1) {
              fit <- mlmi::scoreBased(imps, analysisFun=yonx, scoreFun=yonxscore, vcovtype=2, dfComplete=rep(N-2,3))
            } else {
              fit <- mlmi::scoreBased(imps, analysisFun=xony, scoreFun=xonyscore, vcovtype=2, dfComplete=rep(N-2,3))
            }
            ests[i,j,k,l,sim,methodNum] <- fit$est[2]
            se[i,j,k,l,sim,methodNum] <- fit$var[2,2]^0.5
            df[i,j,k,l,sim,methodNum] <- fit$df[2]
            ci[i,j,k,l,sim,methodNum,] <- c(fit$est[2]-qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5, fit$est[2]+qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5)
            methodNum <- methodNum + 1
            
            #SB full bivariate
            if (yonxOrxony==1) {
              fit <- mlmi::scoreBased(imps, analysisFun=yonx, scoreFun=yonxscore, vcovtype=3, dfComplete=rep(N-2,5))
            } else {
              fit <- mlmi::scoreBased(imps, analysisFun=xony, scoreFun=xonyscore, vcovtype=3, dfComplete=rep(N-2,5))
            }
            ests[i,j,k,l,sim,methodNum] <- fit$est[2]
            se[i,j,k,l,sim,methodNum] <- fit$var[2,2]^0.5
            df[i,j,k,l,sim,methodNum] <- fit$df[2]
            ci[i,j,k,l,sim,methodNum,] <- c(fit$est[2]-qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5, fit$est[2]+qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5)
            methodNum <- methodNum + 1
  
          #MLMI
          imps <- mlmi::normUniImp(simData, y~x, M=M, pd=FALSE)
          
            #WB
            if (yonxOrxony==1) {
              fit <- mlmi::withinBetween(imps, yonx, vcovtype=1, dfComplete=rep(N-2,2))
            } else {
              fit <- mlmi::withinBetween(imps, xony, vcovtype=1, dfComplete=rep(N-2,2))
            }
            ests[i,j,k,l,sim,methodNum] <- fit$est[2]
            se[i,j,k,l,sim,methodNum] <- fit$var[2,2]^0.5
            df[i,j,k,l,sim,methodNum] <- fit$df[2]
            ci[i,j,k,l,sim,methodNum,] <- c(fit$est[2]-qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5, fit$est[2]+qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5)
            methodNum <- methodNum + 1
          
            #WB including residual variance
            if (yonxOrxony==1) {
              fit <- mlmi::withinBetween(imps, yonx, vcovtype=2, dfComplete=rep(N-2,3))
            } else {
              fit <- mlmi::withinBetween(imps, xony, vcovtype=2, dfComplete=rep(N-2,3))
            }
            ests[i,j,k,l,sim,methodNum] <- fit$est[2]
            se[i,j,k,l,sim,methodNum] <- fit$var[2,2]^0.5
            df[i,j,k,l,sim,methodNum] <- fit$df[2]
            ci[i,j,k,l,sim,methodNum,] <- c(fit$est[2]-qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5, fit$est[2]+qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5)
            methodNum <- methodNum + 1
    
            #MLMI WB full bivariate
            if (yonxOrxony==1) {
              fit <- mlmi::withinBetween(imps, yonx, vcovtype=3, dfComplete=rep(N-2,5))
            } else {
              fit <- mlmi::withinBetween(imps, xony, vcovtype=3, dfComplete=rep(N-2,5))
            }
            ests[i,j,k,l,sim,methodNum] <- fit$est[2]
            se[i,j,k,l,sim,methodNum] <- fit$var[2,2]^0.5
            df[i,j,k,l,sim,methodNum] <- fit$df[2]
            ci[i,j,k,l,sim,methodNum,] <- c(fit$est[2]-qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5, fit$est[2]+qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5)
            methodNum <- methodNum + 1
    
            #MLMI SB including residual variance
            if (yonxOrxony==1) {
              fit <- mlmi::scoreBased(imps, analysisFun=yonx, scoreFun=yonxscore, vcovtype=2, dfComplete=rep(N-2,3))
            } else {
              fit <- mlmi::scoreBased(imps, analysisFun=xony, scoreFun=xonyscore, vcovtype=2, dfComplete=rep(N-2,3))
            }
            ests[i,j,k,l,sim,methodNum] <- fit$est[2]
            se[i,j,k,l,sim,methodNum] <- fit$var[2,2]^0.5
            df[i,j,k,l,sim,methodNum] <- fit$df[2]
            ci[i,j,k,l,sim,methodNum,] <- c(fit$est[2]-qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5, fit$est[2]+qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5)
            methodNum <- methodNum + 1
    
            #MLMI SB full bivariate
            if (yonxOrxony==1) {
              fit <- mlmi::scoreBased(imps, analysisFun=yonx, scoreFun=yonxscore, vcovtype=3, dfComplete=rep(N-2,5))
            } else {
              fit <- mlmi::scoreBased(imps, analysisFun=xony, scoreFun=xonyscore, vcovtype=3, dfComplete=rep(N-2,5))
            }
            ests[i,j,k,l,sim,methodNum] <- fit$est[2]
            se[i,j,k,l,sim,methodNum] <- fit$var[2,2]^0.5
            df[i,j,k,l,sim,methodNum] <- fit$df[2]
            ci[i,j,k,l,sim,methodNum,] <- c(fit$est[2]-qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5, fit$est[2]+qt(0.975,df=fit$df[2])*fit$var[2,2]^0.5)
            methodNum <- methodNum + 1
  
          #PDMI boot
          imps <- bootImpute::bootImpute(simData, normUniImp, nBoot=M/2, nImp=2, pd=TRUE, impFormula=y~x, M=1)
          if (yonxOrxony==1) {
            fit <- bootImpute::bootImputeAnalyse(imps, yonxboot, quiet=TRUE)
          } else {
            fit <- bootImpute::bootImputeAnalyse(imps, xonyboot, quiet=TRUE)
          }
          ests[i,j,k,l,sim,methodNum]  <- fit$ests
          se[i,j,k,l,sim,methodNum] <- fit$var^0.5
          df[i,j,k,l,sim,methodNum] <- fit$df
          ci[i,j,k,l,sim,methodNum,] <- c(fit$ci[1], fit$ci[2])
          methodNum <- methodNum + 1
  
          #MLMI boot
          imps <- bootImpute::bootImpute(simData, normUniImp, nBoot=M/2, nImp=2, pd=FALSE, impFormula=y~x, M=1)
          if  (yonxOrxony==1) {
            fit <- bootImpute::bootImputeAnalyse(imps, yonxboot, quiet=TRUE)
          } else {
            fit <- bootImpute::bootImputeAnalyse(imps, xonyboot, quiet=TRUE)
          }
          ests[i,j,k,l,sim,methodNum]  <- fit$ests
          se[i,j,k,l,sim,methodNum] <- fit$var^0.5
          df[i,j,k,l,sim,methodNum] <- fit$df
          ci[i,j,k,l,sim,methodNum,] <- c(fit$ci[1], fit$ci[2])
          methodNum <- methodNum + 1
  
        }
      }
    }
  }
}


save(ests,se,df,ci, file=(paste("./results/simRes_", batch, ".RData", sep="")))
