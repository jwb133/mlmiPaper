#This R program combines the datasets outputted from running the simulations

library(abind)
#setwd("...")

load("simRes_1.RData")
estsC <- ests
seC <- se
ciC <- ci
dfC <- df

nSets <- 1000
for (i in 2:nSets) {
  print(i)
  load(paste("simRes_",i,".RData",sep=""))
  estsC <- abind(estsC, ests, along=5)
  seC <- abind(seC, se, along=5)
  ciC <- abind(ciC, ci, along=5)
  dfC <- abind(dfC, df, along=5)
}

ests <- estsC[1,,,,,]
se <- seC[1,,,,,]
df <- dfC[1,,,,,]
ci <- ciC[1,,,,,,]
save(ests,se,df,ci, file=paste("simResults_yonx_",0, ".RData", sep=""))

ests <- estsC[2,,,,,]
se <- seC[2,,,,,]
df <- dfC[2,,,,,]
ci <- ciC[2,,,,,,]
save(ests,se,df,ci, file=paste("simResults_yonx_",1, ".RData", sep=""))
