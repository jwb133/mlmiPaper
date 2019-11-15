#This R program performs illustrative analyses of wave 2 of the Millenium Cohort Study data
#The datasets that are used can be obtained from the UK's Data Archive: https://www.data-archive.ac.uk/
#The specific collection of datasets that is used is http://doi.org/10.5255/UKDA-SN-5350-4
#The analyses use the CRAN packages mlmi, bootImpute, sandwich, xtable

library(xtable)
#set working directory
#setwd("C:/...")

#construct dataset for analysis from the various wave 2 Millenium Cohort Study datasets

#w2 child assessment dataset
w2_child_assessment <- read.delim("./w2_delim/UKDA-5350-tab/tab/mcs2_child_assessment_data.tab")

#age in categories for Bracken test
#we will treat this as continuous since the categories of equal width of age, and are quite narrow (2.5 months)
table(w2_child_assessment$bbrknage, useNA=c("ifany"))
w2_child_assessment$bbrknage[w2_child_assessment$bbrknage<0] <- NA

#Bracken score-measure of readiness for school
hist(w2_child_assessment$bdsrcs00)
table(w2_child_assessment$bdsrcs00, useNA="ifany")
w2_child_assessment$bdsrcs00[w2_child_assessment$bdsrcs00<0] <- NA
colnames(w2_child_assessment)[which(colnames(w2_child_assessment)=="bdsrcs00")] <- "bracken"

#create subset of variables we will use
w2_child_assessment2 <- subset(w2_child_assessment, select=c("mcsid", "bracken"))
summary(w2_child_assessment2)

#w2 parent interview dataset
w2_parent_interview <- read.delim("./w2_delim/UKDA-5350-tab/tab/mcs2_parent_interview.tab")

#w2 total family income (couple)
table(w2_parent_interview$bminco00, useNA=c("ifany"))
#w2 total family income (lone parent)
table(w2_parent_interview$bmincm00, useNA=c("ifany"))
#combine these two variables
w2_parent_interview$income <- w2_parent_interview$bminco00
w2_parent_interview$income[is.na(w2_parent_interview$income)==TRUE] <- w2_parent_interview$bmincm00[is.na(w2_parent_interview$income)==TRUE]
w2_parent_interview$income[w2_parent_interview$income>18] <- NA

#w2 MAIN Marital status
table(w2_parent_interview$bmfcin00, useNA=c("ifany"))
w2_parent_interview$bmfcin00[w2_parent_interview$bmfcin00<0] <- NA
#combine into married and everything else
w2_parent_interview$bmfcin00[w2_parent_interview$bmfcin00>3] <- 1
w2_parent_interview$bmfcin00[w2_parent_interview$bmfcin00==3] <- 2
colnames(w2_parent_interview)[which(colnames(w2_parent_interview)=="bmfcin00")] <- "maritalStatus"
w2_parent_interview$maritalStatus <- factor(w2_parent_interview$maritalStatus, 
                                            labels=c("Not married", "Married"))

#tenure of current home
table(w2_parent_interview$bmroow00, useNA=c("ifany"))
#combine categories to give: owned, rented, and other
w2_parent_interview$bmroow00[w2_parent_interview$bmroow00<0] <- NA
w2_parent_interview$bmroow00 <- cut(w2_parent_interview$bmroow00,
                                    breaks=c(-Inf, 0, 3.5, 6.5,Inf), 
                                    labels=FALSE)
table(w2_parent_interview$bmroow00, useNA=c("ifany"))
colnames(w2_parent_interview)[which(colnames(w2_parent_interview)=="bmroow00")] <- "tenure"
w2_parent_interview$tenure <- w2_parent_interview$tenure - 1
w2_parent_interview$tenure <- factor(w2_parent_interview$tenure,
                                     labels=c("Owned", "Rented", "Other"))

#any problems with hearing, wave 2
table(w2_parent_interview$bmherpa0, useNA=c("ifany"))
#not applicable to missing
w2_parent_interview$bmherpa0[w2_parent_interview$bmherpa0<0] <- NA
#don't know to missing
w2_parent_interview$bmherpa0[w2_parent_interview$bmherpa0>2] <- NA
colnames(w2_parent_interview)[which(colnames(w2_parent_interview)=="bmherpa0")] <- "hearing"
w2_parent_interview$hearing <- factor(w2_parent_interview$hearing,
                                      labels=c("Yes", "No"))

w2_parent_interview2 <- subset(w2_parent_interview, select=c("mcsid", "income", "maritalStatus", "tenure", "hearing"))
summary(w2_parent_interview2)


#w2 derived dataset
w2_derived <- read.delim("./w2_delim/UKDA-5350-tab/tab/mcs2_derived_variables.tab")

#S2 MAIN DV Respondent's Ethnic Group - 8 category classification
table(w2_derived$BMD08E00, useNA=c("ifany"))
w2_derived$BMD08E00[w2_derived$BMD08E00<0] <- NA
colnames(w2_derived)[which(colnames(w2_derived)=="BMD08E00")] <- "ethnicity"
#dichotomise to binary white/non-white to simplify modelling
w2_derived$ethnicity[w2_derived$ethnicity>1] <- 2
w2_derived$ethnicity <- factor(w2_derived$ethnicity, labels=c("White", "Non-white"))

#main respondent is in work or not
table(w2_derived$BMDWRK00, useNA=c("ifany"))
w2_derived$BMDWRK00[w2_derived$BMDWRK00<0] <- NA
colnames(w2_derived)[which(colnames(w2_derived)=="BMDWRK00")] <- "inWork"
w2_derived$inWork <- factor(w2_derived$inWork, labels=c("In work", "Not in work"))

#number of siblings
table(w2_derived$BDOTHS00, useNA=c("ifany"))
w2_derived$BDOTHS00[w2_derived$BDOTHS00<0] <- NA
#categorise into a 4-level variable (0, 1, 2, 3 or more)
w2_derived$noSiblings <- cut(w2_derived$BDOTHS00, breaks=c(-Inf, 0.5, 1.5, 2.5,Inf), 
                             labels=FALSE)
w2_derived$noSiblings <- factor(w2_derived$noSiblings, labels=c("0", "1", "2", "3 or more"))

table(w2_derived$BMDAGI00, useNA=c("ifany"))
w2_derived$respAge <- w2_derived$BMDAGI00
w2_derived$respAge[w2_derived$respAge<0] <- NA

#rename mcsid to lower case
colnames(w2_derived)[which(colnames(w2_derived)=="MCSID")] <- "mcsid"

w2_derived2 <- subset(w2_derived, select=c("mcsid", "ethnicity", "inWork", "noSiblings", "respAge"))
summary(w2_derived2)

#merge datasets
mergedData <- merge(w2_child_assessment2, w2_parent_interview2, by="mcsid")
mergedData <- merge(mergedData, w2_derived2, by="mcsid")

#remove ID
mergedData <- mergedData[,-c(1)]

dim(mergedData)

#move categorical variables to be first
summary(mergedData)
mergedData <- mergedData[,c(3,4,5,6,7,8,1,2,9)]

#convert factors to numerics for the purposes of later imputation using MLMI
mergedDataNumerics <- mergedData
for (i in 1:6) {
  mergedDataNumerics[,i] <- as.numeric(mergedData[,i])
}

#summarize data with descriptive statistics
colnames(mergedData) <- c("Marital status", "Housing tenure", "History of hearing loss", 
                          "Ethnicity", "Respondent in work", "No. siblings", "Bracken score",
                          "Family income code", "Respondent age")

#find number of rows needed
nRows <- 0
for (i in 1:6) {
  nRows <- nRows + length(levels(mergedData[,i]))+1
}
desTable <- array("", dim=c(nRows,4))
rowNum <- 0
n <- dim(mergedData)[1]
for (i in 1:6) {
  rowNum <- rowNum + 1
  #get counts for this variable
  varcounts <-  table(mergedData[,i],useNA="ifany")
  desTable[rowNum,1] <- colnames(mergedData)[i]
  desTable[rowNum,2] <- levels(mergedData[,i])[1]
  desTable[rowNum,3] <- as.numeric(varcounts[1])
  desTable[rowNum,4] <- signif(100*(as.numeric(varcounts[1])/n),3)
  for (j in 1:(length(levels(mergedData[,i]))-1)) {
    rowNum <- rowNum+1
    desTable[rowNum,2] <- levels(mergedData[,i])[j+1]
    desTable[rowNum,3] <- as.numeric(varcounts[j+1])
    desTable[rowNum,4] <- signif(100*(as.numeric(varcounts[j+1])/n),3)
  }
  #add row for missing values
  rowNum <- rowNum+1
  desTable[rowNum,2] <- "Missing"
  desTable[rowNum,3] <- as.numeric(varcounts[length(varcounts)])
  desTable[rowNum,4] <- signif(100*(as.numeric(varcounts[length(varcounts)])/n),3)
}
colnames(desTable) <- c("Variable", "Level", "Number", "%")
print(xtable(desTable, caption="Distribution of categorical variables in Millenium Cohort Study data"), include.rownames=FALSE)

#now continuous variables
desTable <- array("", dim=c(3,3))
colnames(desTable) <- c("Variable", "Mean (SD)", "No. missing (%)")
for (i in 1:3) {
  desTable[i,1] <- colnames(mergedData[,7:9])[i]
  desTable[i,2] <- paste(round(mean(mergedData[,7:9][,i], na.rm=T),1), " (",
                         round(sd(mergedData[,7:9][,i], na.rm=T),1), ")", sep="")
  desTable[i,3] <- paste(sum(is.na(mergedData[,7:9][,i])), " (", 
                         signif(100*(sum(is.na(mergedData[,7:9][,i]))/n),3), ")", sep="")
}
print(xtable(desTable, caption="Distribution of continuous variables in Millenium Cohort Study data"), include.rownames=FALSE)

#define linear predictor formula of outcome model
outcomeModFormula <- bracken~respAge+income+factor(tenure)+factor(hearing)+factor(ethnicity)+factor(noSiblings)

cca <- lm(outcomeModFormula, data=mergedDataNumerics)

#add sandwich SE version
library(sandwich)
vcovHC(cca)

# MLMI and PDMI
library(mlmi)

#define a function to run the analyses and print results table
runAnalysis <- function(M, B) {
  
  #MLMI
  #we will time all the analyses
  mlmiTimes <- array(0, dim=c(5,2,5))
  mlmiTimes[1,1,] <- proc.time()
  mlmi_imp <- mixImp(mergedDataNumerics, nCat=6, M=M, pd=FALSE, marginsType=1, designType=1, rseed=NULL)
  mlmiTimes[1,2,] <- proc.time()
  
  #WB
  myAnalysis <- function(inputData) {
    mod <- lm(outcomeModFormula, data=inputData)
    myvcov <- as.matrix(Matrix::bdiag(vcov(mod), 2*sigma(mod)^4/(mod$df.residual)))
    list(est=c(mod$coef, sigma(mod)^2), var=myvcov)
  }
  mlmiTimes[2,1,] <- proc.time()
  mlmiWB <- withinBetween(mlmi_imp, myAnalysis)
  mlmiTimes[2,2,] <- proc.time()
  
  #SB
  regscore <- function(inputData, parm) {
    beta <- parm[1:(length(parm)-1)]
    sigmasq <- parm[length(parm)]
    #construct design matrix
    designMat <- model.matrix(formula(paste("~",strsplit(as.character(outcomeModFormula),"~")[[3]])), data=inputData)
    res <- inputData$bracken - designMat %*% beta
    cbind(sweep(designMat, MARGIN=1, res, `*`)/sigmasq, res^2/(2*sigmasq^2)-1/(2*sigmasq))
  }
  
  mlmiTimes[3,1,] <- proc.time()
  mlmiSB <- scoreBased(mlmi_imp, myAnalysis, regscore)
  mlmiTimes[3,2,] <- proc.time()
  
  #bootstrap
  myAnalysisBoot <- function(inputData) {
    mod <- lm(outcomeModFormula, data=inputData)
    c(mod$coef, sigma(mod)^2)
  }
  
  library(bootImpute)
  mlmiTimes[4,1,] <- proc.time()
  #mix::rngseed(412896)
  #set.seed(65423)
  mlmi_boot_imp <- bootImpute(mergedDataNumerics, mixImp, nBoot=B, nImp=2, nCat=6, M=1, pd=FALSE, marginsType=1, designType=1, rseed=NULL)
  mlmiTimes[4,2,] <- proc.time()
  mlmiTimes[5,1,] <- proc.time()
  mlmi_boot_analyse <- bootImputeAnalyse(mlmi_boot_imp, myAnalysisBoot)
  mlmiTimes[5,2,] <- proc.time()
  
  #PDMI
  pdmiTimes <- array(0, dim=c(5,2,5))
  pdmiTimes[1,1,] <- proc.time()
  pdmi_imp <- mixImp(mergedDataNumerics, nCat=6, M=M, pd=TRUE, marginsType=1, designType=1, rseed=NULL)
  pdmiTimes[1,2,] <- proc.time()
  pdmiTimes[2,1,] <- proc.time()
  pdmiWB <- withinBetween(pdmi_imp, myAnalysis)
  pdmiTimes[2,2,] <- proc.time()
  pdmiTimes[3,1,] <- proc.time()
  pdmiSB <- scoreBased(pdmi_imp, myAnalysis, regscore)
  pdmiTimes[3,2,] <- proc.time()
  
  pdmiTimes[4,1,] <- proc.time()
  #mix::rngseed(412896)
  #set.seed(65423)
  pdmi_boot_imp <- bootImpute(mergedDataNumerics, mixImp, nBoot=B, nImp=2, nCat=6, M=1, pd=TRUE, marginsType=1, designType=1, rseed=NULL)
  pdmiTimes[4,2,] <- proc.time()
  pdmiTimes[5,1,] <- proc.time()
  pdmi_boot_analyse <- bootImputeAnalyse(pdmi_boot_imp, myAnalysisBoot)
  pdmiTimes[5,2,] <- proc.time()
  
  #print results tables in LaTeX
  
  #compare timings
  timingTable <- cbind(mlmiTimes[,2,3]-mlmiTimes[,1,3], pdmiTimes[,2,3]-pdmiTimes[,1,3])
  rownames(timingTable) <- c("Imputation", "Within-between SE calculation", "Score-based SE calculation",
                             "Bootstrap imputation", "Bootstrap SE calculation")
  colnames(timingTable) <- c("MLMI", "PDMI")
  xtable(timingTable, digits=1)
  
  estTable <- cbind(mlmiWB$est,pdmiWB$est,mlmi_boot_analyse$ests,pdmi_boot_analyse$ests,myAnalysisBoot(mergedDataNumerics))
  #remove final row which corresponds to residual variance
  estTable <- estTable[1:10,]
  colnames(estTable) <- c("MLMI WB", "PDMI WB", "MLMI Boot", "PDMI Boot", "Complete case")
  rownames(estTable) <- c("Intercept", "Respondent age (years)", "Family income", "Rented housing", "Other housing",
                          "Child hearing loss", "Non-white", "1 sibling", "2 siblings", "3 or more siblings")
  
  seTable <- cbind(diag(mlmiSB$var)^0.5, diag(pdmiSB$var)^0.5, 
                   diag(mlmiWB$var)^0.5, diag(pdmiWB$var)^0.5, mlmi_boot_analyse$var^0.5, pdmi_boot_analyse$var^0.5)
  #remove last row corresponding to residual variance
  seTable <- seTable[1:10,]
  colnames(seTable) <- c("MLMI SB", "PDMI SB", "MLMI WB", "PDMI WB", "MLMI Boot", "PDMI Boot")
  rownames(seTable) <- rownames(estTable)
  
  list(timingTable=timingTable, estTable=estTable, seTable=seTable)
}

#set seeds only once at start
mix::rngseed(412896)
set.seed(65423)

run1 <- runAnalysis(M=100,B=50)
xtable(run1$timingTable, digits=1)
xtable(run1$estTable)
xtable(run1$seTable)

#re-run with larger M and B
run2 <- runAnalysis(M=1000,B=500)
xtable(run2$timingTable, digits=1)
xtable(run2$estTable)
xtable(run2$seTable)