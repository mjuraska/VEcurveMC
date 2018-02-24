library(osDesign)
library(np)
library(chngpt)

# 'tpsPredict' returns predicted values from a model fitted by tps
# columns of newMatrix in the same order as the coefficient vector from tps
tpsPredict <- function(fit, newMatrix){
  linPred <- newMatrix %*% fit$coef
  return(drop(1/(1+exp(-linPred))))
}

# 'hNum' returns function values at s0 of the integrand in the numerator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1, x.male, x.age, x.country are scalars
hNum <- function(s0, s1, x.age, x.country, tpsFit, changePoint, npcdensFit1, npcdensFit2){
  phat.s0 <- tpsPredict(tpsFit, cbind(1, ifelse(s0>changePoint,s0-changePoint,0), as.numeric(x.age==">11"), 
                                      as.numeric(x.country=="COL"), as.numeric(x.country=="HND"), as.numeric(x.country=="IND"), as.numeric(x.country=="MEX"),
                                      as.numeric(x.country=="MYS"), as.numeric(x.country=="PHL"), as.numeric(x.country=="PRI"),as.numeric(x.country=="THA"), 
                                      as.numeric(x.country=="VNM")))
  fhat.s0 <- predict(npcdensFit1, newdata=data.frame(bAUC=s0, AGE=x.age, COUNTRY=x.country, IMPSTLOG.AUCMB=s1))
  ghat.s0 <- predict(npcdensFit2, newdata=data.frame(AGE=x.age, COUNTRY=x.country, IMPSTLOG.AUCMB=s0))
  return(phat.s0*fhat.s0*ghat.s0)
}

# 'phNum' returns function values at s0 of the integrand in the numerator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1, x.age, x.country are scalars
# 'phNum' uses parametric (Gaussian) density estimation, whereas 'hNum' uses nonparametric kernel density estimation
phNum <- function(s0, s1, x.male, x.age, x.country, tpsFit, lmFit1, lmFit2){
  phat.s0 <- tpsPredict(tpsFit, cbind(1, s0, as.numeric(x.male==1), as.numeric(x.age==">11"), 
                                      as.numeric(x.country=="COL"), as.numeric(x.country=="HND"), as.numeric(x.country=="MEX"),
                                      as.numeric(x.country=="PRI")))
  fMean.s0 <- predict(lmFit1, newdata=data.frame(bAUC=s0, MALE=x.male, AGE=x.age, COUNTRY=x.country))
  fSD.s0 <- summary(lmFit1)$sigma
  fhat.s0 <- dnorm(s1, mean=fMean.s0, sd=fSD.s0)
  gMean.s0 <- predict(lmFit2, newdata=data.frame(MALE=x.male, AGE=x.age, COUNTRY=x.country))
  gSD.s0 <- summary(lmFit2)$sigma
  ghat.s0 <- dnorm(s0, mean=gMean.s0, sd=gSD.s0)
  return(phat.s0*fhat.s0*ghat.s0)
}

# 'hDen' returns function values at s0 of the integrand in the denominator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1, x.age, x.country are scalars
hDen <- function(s0, s1, x.age, x.country, npcdensFit1, npcdensFit2){
  fhat.s0 <- predict(npcdensFit1, newdata=data.frame(bAUC=s0, AGE=x.age, COUNTRY=x.country, IMPSTLOG.AUCMB=s1))
  ghat.s0 <- predict(npcdensFit2, newdata=data.frame(AGE=x.age, COUNTRY=x.country, IMPSTLOG.AUCMB=s0))
  return(fhat.s0*ghat.s0)
}

# 'phDen' returns function values at s0 of the integrand in the denominator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1, x.male, x.age, x.country are scalars
# 'phDen' uses parametric (Gaussian) density estimation, whereas 'hDen' uses nonparametric kernel density estimation
phDen <- function(s0, s1, x.male, x.age, x.country, lmFit1, lmFit2){
  fMean.s0 <- predict(lmFit1, newdata=data.frame(bAUC=s0, MALE=x.male, AGE=x.age, COUNTRY=x.country))
  fSD.s0 <- summary(lmFit1)$sigma
  fhat.s0 <- dnorm(s1, mean=fMean.s0, sd=fSD.s0)
  gMean.s0 <- predict(lmFit2, newdata=data.frame(MALE=x.male, AGE=x.age, COUNTRY=x.country))
  gSD.s0 <- summary(lmFit2)$sigma
  ghat.s0 <- dnorm(s0, mean=gMean.s0, sd=gSD.s0)
  return(fhat.s0*ghat.s0)
}

# 'propX' returns the sample proportion of subjects in the (x.age, x.country) category in 'data'
propX <- function(x.age, x.country, data){
  return(NROW(subset(data, AGE==x.age & COUNTRY==x.country))/NROW(data))
}

# 'riskP' returns the value of risk_{(0)}(s1)
# s1 is a scalar
riskP <- function(s1, data, tpsFit, npcdensFit1, npcdensFit2, changePoint){
  den <- num <- 0
  UL <- max(data$IMPSTLOG.AUCMB, na.rm=TRUE) + 0.2 # if integration over (0,Inf) fails, use (0,UL)
  
  for (age in levels(as.factor(data$AGE))){
    for (country in levels(as.factor(data$COUNTRY))){
      pX <- propX(x.age=age, x.country=country, data=data)
      hNumInt <- try(integrate(hNum, 0, Inf, s1=s1, x.age=age, x.country=country, tpsFit=tpsFit, changePoint=changePoint, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000)$value, silent=TRUE)
      if (inherits(hNumInt, 'try-error')){
        num <- num + pX*integrate(hNum, 0, UL, s1=s1, x.age=age, x.country=country, tpsFit=tpsFit, changePoint=changePoint, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000)$value
      } else {
        num <- num + pX*hNumInt
      }
      den <- den + pX*integrate(hDen, 0, UL, s1=s1, x.age=age, x.country=country, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000, rel.tol=30*.Machine$double.eps^0.25)$value
    }
  }
  
  return(num/den)
}

# 'riskP' returns the value of risk_{(0)}(s1)
# s1 is a scalar
# 'priskP' uses parametric (Gaussian) density estimation, whereas 'riskP' uses nonparametric kernel density estimation
priskP <- function(s1, data, tpsFit, lmFit1, lmFit2){
  den <- num <- 0
  UL <- max(data$IMPSTLOG.AUCMB, na.rm=TRUE) + 0.2 # if integration over (0,Inf) fails, use (0,UL)
  for (male in 0:1){
    for (age in levels(as.factor(data$AGE))){
      for (country in levels(as.factor(data$COUNTRY))){
        pX <- propX(x.male=male, x.age=age, x.country=country, data=data)
        hNumInt <- try(integrate(phNum, 0, Inf, s1=s1, x.male=male, x.age=age, x.country=country, tpsFit=tpsFit, lmFit1=lmFit1, lmFit2=lmFit2, subdivisions=2000)$value, silent=TRUE)
        if (inherits(hNumInt, 'try-error')){
          num <- num + pX*integrate(phNum, 0, UL, s1=s1, x.male=male, x.age=age, x.country=country, tpsFit=tpsFit, lmFit1=lmFit1, lmFit2=lmFit2, subdivisions=2000)$value
        } else {
          num <- num + pX*hNumInt
        }
        den <- den + pX*integrate(phDen, 0, UL, s1=s1, x.male=male, x.age=age, x.country=country, lmFit1=lmFit1, lmFit2=lmFit2, subdivisions=2000, rel.tol=30*.Machine$double.eps^0.25)$value
      }
    }
  }
  return(num/den)
}

# 'riskV' returns the value of risk_{(1)}(s1)
# s1 is a scalar
riskV <- function(s1, data, dataI, markerName, changePoint){
  if (markerName %in% c("AUC","Min")){
    nVControls <- NROW(subset(data, VACC==1 & ofstatus_m13==0))    
  } else {
    nVControls <- NROW(subset(data, VACC==1 & ofstatus_m0==0))    
  }
  nVCases <- NROW(subset(data, VACC==1 & ofstatus_m13==1))
  group <- rep(1, NROW(subset(dataI, VACC==1 & !is.na(ofstatus_m13))))
  fit <- tps(ofstatus_m13 ~ IMPSTLOG.AUCMB.trunc, data=subset(dataI, VACC==1 & !is.na(ofstatus_m13)), nn0=nVControls, nn1=nVCases, group=group, method="PL", cohort=TRUE)
  return(tpsPredict(fit, cbind(1, ifelse(s1>changePoint,s1-changePoint,0))))
}

# 'risk' returns the estimates of risk in each treatment arm for a given s1
# s1 is a scalar
risk <- function(s1, data, dataI, tpsFit, npcdensFit1, npcdensFit2, markerName, changePoint){
  risk1 <- riskV(s1, data, dataI, markerName, changePoint)
  risk0 <- riskP(s1, data, tpsFit, npcdensFit1, npcdensFit2, changePoint)
  return(list(plaRisk=risk0, txRisk=risk1))
}

# 'pVE' returns the value of VE(s1)
# s1 is a scalar
pVE <- function(s1, data, dataI, tpsFit, lmFit1, lmFit2, markerName){ 1 - riskV(s1, data, dataI, markerName)/priskP(s1, data, tpsFit, lmFit1, lmFit2) }

# 'bRiskCurve' returns a list; one list component is a matrix (or vector in case iter=1) with rows 
# being the estimated VE(s1) curves based on bootstrap samples
# 'data' is assumed to be the ITT set at-risk at month 13 with no prior infection
# 'markerName' is one of "AUC", "S1", "S2", "S3", "S4"
# 'iter' is the number of bootstrap iterations
# 'saveFile' is the name of the .RData file where the output list will be stored
# 'saveDir' specifies the output directory
# 'seed' is an integer for set.seed()
bRiskCurve <- function(data, markerName, iter, saveFile=NULL, saveDir=NULL, seed=NULL){
  # so that the below generic code can be used for each marker
  if (markerName!="AUC"){
    if (markerName=="Min"){  # assumes that there exist variables 'bMin' and 'IMPSTLOG.Min'
      data$bAUC <- data$bMin
      data$IMPSTLOG.AUCMB <- data$IMPSTLOG.Min
    } else {
      data$ofstatus_m0 <- data[,paste0(tolower(markerName),"fstatus_m0")]
      data$ofstatus_m13 <- data[,paste0(tolower(markerName),"fstatus_m13")]
      data$oftime_m13 <- data[,paste0(tolower(markerName),"ftime_m13")]
      data$bAUC <- data[,paste0("b",markerName)]
      data$IMPSTLOG.AUCMB <- data[,paste0("IMPSTLOG.Sero",substr(markerName, start=2, stop=2))] 
    }       
  }
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  
  # compute the change point for the association of the marker with the dengue disease risk in both the placebo and vaccine groups and consider their minimum
  # as the change point in subsequent logistic regression models
  # use the same change point for each bootstrap sample
  cpointP <- chngptm(formula.1=ofstatus_m13 ~ factor(AGE) + factor(COUNTRY), formula.2=~IMPSTLOG.AUCMB, data=subset(dataI, VACC==0 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
  cpointV <- chngptm(formula.1=ofstatus_m13 ~ 1, formula.2=~IMPSTLOG.AUCMB, data=subset(dataI, VACC==1 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
  cpoint <- min(cpointP, cpointV)
  
  # the truncated version of IMPSTLOG.AUCMB based on the change point 'cpoint' is used in all subsequent logistic regression models of the dengue disease risk
  data$IMPSTLOG.AUCMB.trunc <- with(data, ifelse(IMPSTLOG.AUCMB>cpoint, IMPSTLOG.AUCMB-cpoint, 0))
  
  # re-extract the immunogenicity set to include IMPSTLOG.AUCMB.trunc
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  
  # extract subsets of controls ('dataControls') and cases ('dataCases') to be used for resampling
  # in addition, within each treatment group in the immunogenicity set, delete cases to recover
  # the case:control ratio in the ITT set at-risk at month 13 with no prior infection
  dataControls <- subset(data, ofstatus_m0==0)
  nPControlsI <- NROW(dataPControlsI <- subset(dataI, VACC==0 & ofstatus_m0==0))
  nVControlsI <- NROW(dataVControlsI <- subset(dataI, VACC==1 & ofstatus_m0==0))
  
  dataCases <- subset(data, ofstatus_m13==1)  
  nPCasesI <- NROW(dataPCasesI <- subset(dataI, VACC==0 & ofstatus_m13==1))
  nVCasesI <- NROW(dataVCasesI <- subset(dataI, VACC==1 & ofstatus_m13==1))
  
  nPControls <- NROW(subset(dataControls, VACC==0))
  nVControls <- NROW(subset(dataControls, VACC==1))
  nPCases <- NROW(subset(dataCases, VACC==0))
  nVCases <- NROW(subset(dataCases, VACC==1))  
  
  # within each treatment group, calculate the number of cases in the immunogenicity set needed to achieve
  # the correct case:control ratio
  nPCasesInew <- nPCases * nPControlsI / nPControls
  nVCasesInew <- nVCases * nVControlsI / nVControls
  
  # within each treatment group, sample as many cases in the immunogenicity set as needed to achieve 
  # the correct case:control ratio
  dataPIcorrectRatio <- rbind(dataPControlsI, dataPCasesI[sample(1:nPCasesI, nPCasesInew),])
  dataVIcorrectRatio <- rbind(dataVControlsI, dataVCasesI[sample(1:nVCasesI, nVCasesInew),])
  rm(dataPControlsI); rm(dataPCasesI); rm(dataVControlsI); rm(dataVCasesI)
  
  # the overall numbers of controls and cases for resampling
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
  # estimate the optimal bandwidths for kernel density estimation and use these in each bootstrap run
  # RUNS SLOWLY!
  fbw <- npcdensbw(IMPSTLOG.AUCMB ~ bAUC + factor(AGE) + factor(COUNTRY), data=dataVIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  gbw <- npcdensbw(IMPSTLOG.AUCMB ~ factor(AGE) + factor(COUNTRY), data=dataPIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")  
  rm(dataPIcorrectRatio); rm(dataVIcorrectRatio)
  
  if(!is.null(seed)){ set.seed(seed) }
  bSampleControls <- matrix(sample(1:nControls, nControls*iter, replace=TRUE), nrow=nControls, ncol=iter)
  bSampleCases <- matrix(sample(1:nCases, nCases*iter, replace=TRUE), nrow=nCases, ncol=iter)
  
  markerVals <- seq(min(dataI$IMPSTLOG.AUCMB), max(dataI$IMPSTLOG.AUCMB), by=0.05)
  
  # 'bVE' is a list each of whose components is also a list with components:
  # 'VEcurve' - a vector with VE(s) estimates (a single curve)
  # 'bnI'     -  the size of the bootstrapped immunogenicity set
  bRiskCurveList <- lapply(1:iter, function(i){
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,i],], dataCases[bSampleCases[,i],])
    # extract the bootstrapped immunogenicity set
    bdataI <- subset(bdata, !is.na(IMPSTLOG.AUCMB))
    
    # compute the change point for the association of the marker with the dengue disease risk in both the placebo and vaccine groups and consider their minimum
    # as the change point in subsequent logistic regression models
    # use the same change point for each bootstrap sample
    bcpointP <- chngptm(formula.1=ofstatus_m13 ~ factor(AGE) + factor(COUNTRY), formula.2=~IMPSTLOG.AUCMB, data=subset(bdataI, VACC==0 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
    bcpointV <- chngptm(formula.1=ofstatus_m13 ~ 1, formula.2=~IMPSTLOG.AUCMB, data=subset(bdataI, VACC==1 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
    bcpoint <- min(bcpointP, bcpointV)
    
    # the truncated version of IMPSTLOG.AUCMB based on the change point 'cpoint' is used in all subsequent logistic regression models of the dengue disease risk
    bdata$IMPSTLOG.AUCMB.trunc <- with(bdata, ifelse(IMPSTLOG.AUCMB>bcpoint, IMPSTLOG.AUCMB-bcpoint, 0))
    
    # re-extract the immunogenicity set to include IMPSTLOG.AUCMB.trunc
    bdataI <- subset(bdata, !is.na(IMPSTLOG.AUCMB))
    
    bdataControls <- subset(bdata, ofstatus_m0==0)
    nPControlsI <- NROW(bdataPControlsI <- subset(bdataI, VACC==0 & ofstatus_m0==0))
    nVControlsI <- NROW(bdataVControlsI <- subset(bdataI, VACC==1 & ofstatus_m0==0))
    
    bdataCases <- subset(bdata, ofstatus_m13==1)    
    nPCasesI <- NROW(bdataPCasesI <- subset(bdataI, VACC==0 & ofstatus_m13==1))
    nVCasesI <- NROW(bdataVCasesI <- subset(bdataI, VACC==1 & ofstatus_m13==1))
    
    nPControls <- NROW(subset(bdataControls, VACC==0))
    nVControls <- NROW(subset(bdataControls, VACC==1))
    nPCases <- NROW(subset(bdataCases, VACC==0))
    nVCases <- NROW(subset(bdataCases, VACC==1))    
    
    # within each treatment group, calculate the number of cases in the bootstrapped immunogenicity set 
    # needed to achieve the correct case:control ratio
    nPCasesInew <- nPCases * nPControlsI / nPControls
    nVCasesInew <- nVCases * nVControlsI / nVControls
    
    # within each treatment group, sample as many cases in the bootstrapped immunogenicity set as needed 
    # to achieve the correct case:control ratio
    bdataPIcorrectRatio <- rbind(bdataPControlsI, bdataPCasesI[sample(1:nPCasesI, nPCasesInew),])
    bdataVIcorrectRatio <- rbind(bdataVControlsI, bdataVCasesI[sample(1:nVCasesI, nVCasesInew),])
    rm(bdataPControlsI); rm(bdataPCasesI); rm(bdataVControlsI); rm(bdataVCasesI)
    
    group <- rep(1, NROW(subset(bdataI, VACC==0)))
    
    # weighted logistic regression model using the placebo group in the bootstrapped immunogenicity set
    fit1 <- tps(ofstatus_m13 ~ IMPSTLOG.AUCMB.trunc + factor(AGE) + factor(COUNTRY), data=subset(bdataI, VACC==0), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
    
    # kernel density estimator for f(s1|S_base=sbase, X=x) using the vaccine group in the bootstrapped immunogenicity set
    bfbw <- npcdensbw(IMPSTLOG.AUCMB ~ bAUC + factor(AGE) + factor(COUNTRY), data=bdataVIcorrectRatio, bws=fbw, bandwidth.compute=FALSE)
    fhat <- npcdens(bfbw)
    
    # kernel density estimator for g(s0|X=x) using the placebo group in the bootstrapped immunogenicity set
    bgbw <- npcdensbw(IMPSTLOG.AUCMB ~ factor(AGE) + factor(COUNTRY), data=bdataPIcorrectRatio, bws=gbw, bandwidth.compute=FALSE)    
    ghat <- npcdens(bgbw)
    
    # a single bootstrap VE(s) curve, risk1(s) curve, and risk0(s) curve
    curves <- lapply(markerVals, function(s){ risk(s, data, dataI, fit1, fhat, ghat, markerName, cpoint) })
    plaRiskCurve <- sapply(curves, "[[", "plaRisk")
    txRiskCurve <- sapply(curves, "[[", "txRisk")
    
    return(list(plaRiskCurve=plaRiskCurve, txRiskCurve=txRiskCurve))
  })
  
  # cbind all bootstrap risk curves
  plaRiskCurveBootEst <- drop(do.call(cbind, lapply(bRiskCurveList,"[[","plaRiskCurve")))
  txRiskCurveBootEst <- drop(do.call(cbind, lapply(bRiskCurveList,"[[","txRiskCurve")))
  # bList <- list(markerVals=markerVals, plaRiskCurveBootEst=plaRiskCurveBootEst, txRiskCurveBootEst=txRiskCurveBootEst, bIdxControls=bSampleControls, bIdxCases=bSampleCases, 
  #               cpointP=cpointP, cpointV=cpointV, fOptBandwidths=fbw, gOptBandwidths=gbw, seed=seed)
  bList <- list(markerVals=markerVals, plaRiskCurveBootEst=plaRiskCurveBootEst, txRiskCurveBootEst=txRiskCurveBootEst)
  
  if (!is.null(saveFile)){
    save(bList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n")
  }
  
  return(invisible(bList))
}

# 'riskCurve' returns the estimated P(Y(0)=1|S(1)=s1) and P(Y(1)=1|S(1)=s1) on a grid of s1 values
# 'data' is assumed to be the ITT set at-risk at month 13 with no prior infection (the target population)
riskCurve <- function(data, markerName, saveFile=NULL, saveDir=NULL){
  if (markerName!="AUC"){
    if (markerName=="Min"){  # assumes that there exist variables 'bMin' and 'IMPSTLOG.Min'
      data$bAUC <- data$bMin
      data$IMPSTLOG.AUCMB <- data$IMPSTLOG.Min
    } else {
      data$ofstatus_m0 <- data[,paste0(tolower(markerName),"fstatus_m0")]
      data$ofstatus_m13 <- data[,paste0(tolower(markerName),"fstatus_m13")]
      data$oftime_m13 <- data[,paste0(tolower(markerName),"ftime_m13")]
      data$bAUC <- data[,paste0("b",markerName)]
      data$IMPSTLOG.AUCMB <- data[,paste0("IMPSTLOG.Sero",substr(markerName, start=2, stop=2))]    
    }    
  }
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  
  # extract subsets of controls ('dataControls') and cases ('dataCases') to be used for resampling
  # in addition, within each treatment group in the immunogenicity set, delete cases to recover
  # the case:control ratio in the ITT set at-risk at month 13 with no prior infection
  dataControls <- subset(data, ofstatus_m0==0)
  nPControlsI <- NROW(dataPControlsI <- subset(dataI, VACC==0 & ofstatus_m0==0))
  nVControlsI <- NROW(dataVControlsI <- subset(dataI, VACC==1 & ofstatus_m0==0))  
  
  dataCases <- subset(data, ofstatus_m13==1)  
  nPCasesI <- NROW(dataPCasesI <- subset(dataI, VACC==0 & ofstatus_m13==1))
  nVCasesI <- NROW(dataVCasesI <- subset(dataI, VACC==1 & ofstatus_m13==1))
  
  nPControls <- NROW(subset(dataControls, VACC==0))
  nVControls <- NROW(subset(dataControls, VACC==1))
  nPCases <- NROW(subset(dataCases, VACC==0))
  nVCases <- NROW(subset(dataCases, VACC==1))
  
  # within each treatment group, calculate the number of cases in the immunogenicity set needed to achieve
  # the correct case:control ratio
  nPCasesInew <- nPCases * nPControlsI / nPControls
  nVCasesInew <- nVCases * nVControlsI / nVControls
  
  # within each treatment group, sample as many cases in the immunogenicity set as needed to achieve 
  # the correct case:control ratio
  dataPIcorrectRatio <- rbind(dataPControlsI, dataPCasesI[sample(1:nPCasesI, nPCasesInew),])
  dataVIcorrectRatio <- rbind(dataVControlsI, dataVCasesI[sample(1:nVCasesI, nVCasesInew),])
  rm(dataPControlsI); rm(dataPCasesI); rm(dataVControlsI); rm(dataVCasesI)
  
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
  # estimate the optimal bandwidths
  # 'AGE' is categorical: <=11 and >11
  fbw <- npcdensbw(IMPSTLOG.AUCMB ~ bAUC + factor(AGE) + factor(COUNTRY), data=dataVIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  gbw <- npcdensbw(IMPSTLOG.AUCMB ~ factor(AGE) + factor(COUNTRY), data=dataPIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  
  group <- rep(1, NROW(subset(dataI, VACC==0 & !is.na(ofstatus_m13))))
  
  # compute the change point for the association of the marker with the dengue disease risk in both the placebo and vaccine group and consider their minimum
  # as the change point in subsequent logistic regression models
  cpointP <- chngptm(formula.1=ofstatus_m13 ~ factor(AGE) + factor(COUNTRY), formula.2=~IMPSTLOG.AUCMB, data=subset(dataI, VACC==0 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
  cpointV <- chngptm(formula.1=ofstatus_m13 ~ 1, formula.2=~IMPSTLOG.AUCMB, data=subset(dataI, VACC==1 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
  cpoint <- min(cpointP, cpointV)
  
  # the truncated version of IMPSTLOG.AUCMB based on the change point 'cpoint' is used in all subsequent logistic regression models of the dengue disease risk
  data$IMPSTLOG.AUCMB.trunc <- with(data, ifelse(IMPSTLOG.AUCMB>cpoint, IMPSTLOG.AUCMB-cpoint, 0))
  
  # re-extract the immunogenicity set to include IMPSTLOG.AUCMB.trunc
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  
  # weighted logistic regression model using the placebo group in the immunogenicity set
  fit1 <- tps(ofstatus_m13 ~ IMPSTLOG.AUCMB.trunc + factor(AGE) + factor(COUNTRY), data=subset(dataI, VACC==0 & !is.na(ofstatus_m13)), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
  
  # kernel density estimator for f(s1|S_base=sbase, X=x) using the vaccine group in the immunogenicity set
  fhat <- npcdens(fbw)
  
  # kernel density estimator for g(s0|X=x) using the placebo group in the immunogenicity set
  ghat <- npcdens(gbw)
  
  # a single VE(s) curve
  markerVals <- seq(min(dataI$IMPSTLOG.AUCMB), max(dataI$IMPSTLOG.AUCMB), by=0.05)
  
  curves <- lapply(markerVals, function(s){ risk(s, data, dataI, fit1, fhat, ghat, markerName, cpoint) })
  plaRiskCurve <- sapply(curves, "[[", "plaRisk")
  txRiskCurve <- sapply(curves, "[[", "txRisk")
  
  # the output list
  oList <- list(markerVals=markerVals, plaRiskCurve=plaRiskCurve, txRiskCurve=txRiskCurve, cpointP=cpointP, cpointV=cpointV, fOptBandwidths=fbw, gOptBandwidths=gbw)
  
  if (!is.null(saveFile)){
    save(oList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n")
  }
  
  return(invisible(oList))
}

VEcurveWithLimitedBaselineMarkersInVaccinees <- function(data, markerVals=NULL, markerName, nVaccineesWithBaselineMarker, nMCruns, contrast="multiplicative", 
                                                         saveFile=NULL, saveDir=NULL, seed=NULL){
  if (!is.null(seed)){ set.seed(seed) }
  
  if (markerName!="AUC"){
    if (markerName=="Min"){  # assumes that there exist variables 'bMin' and 'IMPSTLOG.Min'
      data$bAUC <- data$bMin
      data$IMPSTLOG.AUCMB <- data$IMPSTLOG.Min
    } else {
      data$ofstatus_m0 <- data[,paste0(tolower(markerName),"fstatus_m0")]
      data$ofstatus_m13 <- data[,paste0(tolower(markerName),"fstatus_m13")]
      data$oftime_m13 <- data[,paste0(tolower(markerName),"ftime_m13")]
      data$bAUC <- data[,paste0("b",markerName)]
      data$IMPSTLOG.AUCMB <- data[,paste0("IMPSTLOG.Sero",substr(markerName, start=2, stop=2))]    
    }    
  }
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  
  # extract subsets of controls ('dataControls') and cases ('dataCases') to be used for resampling
  # in addition, within each treatment group in the immunogenicity set, delete cases to recover
  # the case:control ratio in the ITT set at-risk at month 13 with no prior infection
  dataControls <- subset(data, ofstatus_m0==0)
  nPControlsI <- NROW(dataPControlsI <- subset(dataI, VACC==0 & ofstatus_m0==0))
  nVControlsI <- NROW(dataVControlsI <- subset(dataI, VACC==1 & ofstatus_m0==0))  
  
  dataCases <- subset(data, ofstatus_m13==1)  
  nPCasesI <- NROW(dataPCasesI <- subset(dataI, VACC==0 & ofstatus_m13==1))
  nVCasesI <- NROW(dataVCasesI <- subset(dataI, VACC==1 & ofstatus_m13==1))
  
  nPControls <- NROW(subset(dataControls, VACC==0))
  nVControls <- NROW(subset(dataControls, VACC==1))
  nPCases <- NROW(subset(dataCases, VACC==0))
  nVCases <- NROW(subset(dataCases, VACC==1))
  
  # within each treatment group, calculate the number of cases in the immunogenicity set needed to achieve
  # the correct case:control ratio
  nPCasesInew <- nPCases * nPControlsI / nPControls
  nVCasesInew <- nVCases * nVControlsI / nVControls
  
  # within each treatment group, sample as many cases in the immunogenicity set as needed to achieve 
  # the correct case:control ratio
  dataPIcorrectRatio <- rbind(dataPControlsI, dataPCasesI[sample(1:nPCasesI, nPCasesInew),])
  dataVIcorrectRatio <- rbind(dataVControlsI, dataVCasesI[sample(1:nVCasesI, nVCasesInew),])
  rm(dataPControlsI); rm(dataPCasesI); rm(dataVControlsI); rm(dataVCasesI)
  
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
  # estimate the optimal bandwidths
  gbw <- npcdensbw(IMPSTLOG.AUCMB ~ MALE + factor(AGE) + factor(COUNTRY), data=dataPIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  
  # kernel density estimator for g(s0|X=x) using the placebo group in the immunogenicity set
  ghat <- npcdens(gbw)
  
  group <- rep(1, NROW(subset(dataI, VACC==0 & !is.na(ofstatus_m13))))
  
  # compute the change point for the association of the marker with the dengue disease risk in both the placebo and vaccine group and consider their minimum
  # as the change point in subsequent logistic regression models
  cpointP <- chngptm(formula.1=ofstatus_m13 ~ MALE + factor(AGE) + factor(COUNTRY), formula.2=~IMPSTLOG.AUCMB, data=subset(dataI, VACC==0 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
  cpointV <- chngptm(formula.1=ofstatus_m13 ~ 1, formula.2=~IMPSTLOG.AUCMB, data=subset(dataI, VACC==1 & !is.na(ofstatus_m13)), family="binomial", type="hinge", prob.weights=wts)$coefficients["chngpt"]
  cpoint <- min(cpointP, cpointV)
  
  # the truncated version of IMPSTLOG.AUCMB based on the change point 'cpoint' is used in all subsequent logistic regression models of the dengue disease risk
  data$IMPSTLOG.AUCMB.trunc <- with(data, ifelse(IMPSTLOG.AUCMB>cpoint, IMPSTLOG.AUCMB-cpoint, 0))
  
  # re-extract the immunogenicity set
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  
  # the complement of the IS
  dataNonI <- subset(data, is.na(IMPSTLOG.AUCMB))
  
  # weighted logistic regression model using the placebo group in the immunogenicity set
  fit1 <- tps(ofstatus_m13 ~ IMPSTLOG.AUCMB.trunc + MALE + factor(AGE) + factor(COUNTRY), data=subset(dataI, VACC==0 & !is.na(ofstatus_m13)), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
  
  # a single VE(s) curve
  if (is.null(markerVals)){ markerVals <- seq(min(dataI$IMPSTLOG.AUCMB), max(dataI$IMPSTLOG.AUCMB), by=0.1) }
  
  # calculate the new number of cases in the vaccine group with a measured baseline marker maintaining the case:control ratio of measured baseline markers
  # in the originanl data
  nVCasesWithBaselineMarker <- NROW(dataVCasesWithBaselineMarker <- subset(dataI, VACC==1 & ofstatus_m13==1 & !is.na(bAUC)))
  nVControlsWithBaselineMarker <- NROW(dataVControlsWithBaselineMarker <- subset(dataI, VACC==1 & ofstatus_m0==0 & !is.na(bAUC)))
  if (nVaccineesWithBaselineMarker > nVCasesWithBaselineMarker + nVControlsWithBaselineMarker){ stop("The 'nVaccineesWithBaselineMarker' argument exceeds the number of measured baseline markers among vaccine recipients in the original data.") }
  nVCasesWithBaselineMarker.new <- ceiling(nVaccineesWithBaselineMarker * nVCasesWithBaselineMarker/(nVCasesWithBaselineMarker + nVControlsWithBaselineMarker))
  nVControlsWithBaselineMarker.new <- nVaccineesWithBaselineMarker - nVCasesWithBaselineMarker.new
  
  # initialize the output matrices
  Risk0MC <- Risk1MC <- VEcurvesMC <- matrix(0, nrow=nMCruns, ncol=length(markerVals))
  
  for (b in 1:nMCruns){
    # sample subsets of cases and controls in the vaccine group with measured baseline markers
    idxSampledCases <- sample(1:nVCasesWithBaselineMarker, size=nVCasesWithBaselineMarker.new)
    idxSampledControls <- sample(1:nVControlsWithBaselineMarker, size=nVControlsWithBaselineMarker.new)
    dataVCasesWithBaselineMarker.new <- dataVCasesWithBaselineMarker[idxSampledCases,]
    dataVControlsWithBaselineMarker.new <- dataVControlsWithBaselineMarker[idxSampledControls,]
    
    # take the complement of the sampled subsets among vaccinees with measured baseline markers and add noise to their baseline markers
    dataVaccineesWithMismeasuredBaselineMarker.new <- rbind(dataVCasesWithBaselineMarker[-idxSampledCases,], dataVControlsWithBaselineMarker[-idxSampledControls,])
    # calculate the standard deviation of the error term that corresponds to R^2=0.7
    x <- dataVaccineesWithMismeasuredBaselineMarker.new$bAUC
    sigma <- 1.2 * sqrt((1-0.7) * var(x))
    dataVaccineesWithMismeasuredBaselineMarker.new$bAUC <- dataVaccineesWithMismeasuredBaselineMarker.new$bAUC + rnorm(NROW(dataVaccineesWithMismeasuredBaselineMarker.new), 0, sd=sigma)
    # replace values below the LLOQ with log10(5)
    dataVaccineesWithMismeasuredBaselineMarker.new$bAUC <- ifelse(dataVaccineesWithMismeasuredBaselineMarker.new$bAUC < 1, log10(5), dataVaccineesWithMismeasuredBaselineMarker.new$bAUC)
    
    # assemble the new immunogenicity data
    dataVaccineesWithoutBaselineMarker <- subset(dataI, VACC==1 & is.na(bAUC))
    dataPlacebos <- subset(dataI, VACC==0)
    dataI.new <- rbind(dataVCasesWithBaselineMarker.new, dataVControlsWithBaselineMarker.new, dataVaccineesWithMismeasuredBaselineMarker.new, dataVaccineesWithoutBaselineMarker, dataPlacebos)
    data.new <- rbind(dataI.new, dataNonI)
    
    # within the vaccine group in the immunogenicity set, delete cases to recover
    # the case:control ratio in the ITT set at-risk at month 13 with no prior infection
    dataVControlsI.new <- subset(dataI.new, VACC==1 & ofstatus_m0==0)
    dataVCasesI.new <- subset(dataI.new, VACC==1 & ofstatus_m13==1)
    
    # within the vaccine group, sample as many cases in the immunogenicity set as needed to achieve the correct case:control ratio
    dataVIcorrectRatio.new <- rbind(dataVControlsI.new, dataVCasesI.new[sample(1:nVCasesI, nVCasesInew),])
    
    # estimate the optimal bandwidths
    fbw <- npcdensbw(IMPSTLOG.AUCMB ~ bAUC + MALE + factor(AGE) + factor(COUNTRY), data=dataVIcorrectRatio.new, cxkertype="epanechnikov", cykertype="epanechnikov")
    
    # kernel density estimator for f(s1|S_base=sbase, X=x) using the vaccine group in the immunogenicity set
    fhat <- npcdens(fbw)
    
    curves <- lapply(markerVals, function(s){ VE(s, data.new, dataI.new, fit1, fhat, ghat, markerName, cpoint, contrast) })
    VEcurve <- sapply(curves, "[[", "VE")
    risk1 <- sapply(curves, "[[", "risk1")
    risk0 <- sapply(curves, "[[", "risk0")
    
    VEcurvesMC[b,] <- VEcurve
    Risk1MC[b,] <- risk1
    Risk0MC[b,] <- risk0
  }
  
  oList <- list(markerVals=markerVals, risk1=Risk1MC, risk0=Risk0MC, VE=VEcurvesMC, cpointP=cpointP, cpointV=cpointV, gOptBandwidths=gbw)
  
  if (!is.null(saveFile)){
    save(oList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n")
  } else {
    return(oList)
  }  
}


# 'pVEcurve' is a "parametric" version of 'VEcurve' with parametric estimates of conditional densities
# 'data' is assumed to be the ITT set
pVEcurve <- function(data, markerName, saveFile=NULL, saveDir=NULL){
  if (markerName!="AUC"){
    data$ofstatus_m0 <- data[,paste0(tolower(markerName),"fstatus_m0")]
    data$ofstatus_m13 <- data[,paste0(tolower(markerName),"fstatus_m13")]
    data$oftime_m13 <- data[,paste0(tolower(markerName),"ftime_m13")]
    data$bAUC <- data[,paste0("b",markerName)]
    data$IMPSTLOG.AUCMB <- data[,paste0("IMPSTLOG.Sero",substr(markerName, start=2, stop=2))]    
  }
  # extract the ITT set at-risk at month 13, which is the target population for inference about VE(s1)
  data <- subset(data, oftime_m13>0)
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  
  # calculate the sampling weights
  if (markerName=="AUC"){
    nPControls <- NROW(subset(data, VACC==0 & ofstatus_m13==0))
    wtPControls <- NROW(subset(data, VACC==0 & ofstatus_m13==0))/NROW(subset(dataI, VACC==0 & ofstatus_m13==0))
    wtVControls <- NROW(subset(data, VACC==1 & ofstatus_m13==0))/NROW(subset(dataI, VACC==1 & ofstatus_m13==0))
  } else {
    nPControls <- NROW(subset(data, VACC==0 & ofstatus_m0==0))
    wtPControls <- NROW(subset(data, VACC==0 & ofstatus_m0==0))/NROW(subset(dataI, VACC==0 & ofstatus_m0==0))
    wtVControls <- NROW(subset(data, VACC==1 & ofstatus_m0==0))/NROW(subset(dataI, VACC==1 & ofstatus_m0==0))
  }  
  nPCases <- NROW(subset(data, VACC==0 & ofstatus_m13==1))
  wtPCases <- NROW(subset(data, VACC==0 & ofstatus_m13==1))/NROW(subset(dataI, VACC==0 & ofstatus_m13==1))
  wtVCases <- NROW(subset(data, VACC==1 & ofstatus_m13==1))/NROW(subset(dataI, VACC==1 & ofstatus_m13==1))
  group <- rep(1, NROW(subset(dataI, VACC==0 & !is.na(ofstatus_m13))))
  
  # weighted logistic regression model using the placebo group in the immunogenicity set
  fit1 <- tps(ofstatus_m13 ~ IMPSTLOG.AUCMB + MALE + factor(AGE) + factor(COUNTRY), data=subset(dataI, VACC==0 & !is.na(ofstatus_m13)), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
  
  # normal density estimator for f(s1|S_base=sbase, X=x) using the vaccine group in the immunogenicity set with baseline markers
  # sampling weights are incorporated
  dataB <- subset(dataI, VACC==1 & !is.na(bAUC))
  fLM <- lm(IMPSTLOG.AUCMB ~ ns(bAUC) + MALE + factor(AGE) + factor(COUNTRY), data=dataB, weights=ifelse(dataB$ofstatus_m13==1, wtVCases, wtVControls))
  
  # kernel density estimator for g(s0|X=x) using the placebo group in the immunogenicity set
  dataB <- subset(dataI, VACC==0)
  gLM <- lm(IMPSTLOG.AUCMB ~ MALE + factor(AGE) + factor(COUNTRY), data=dataB, weights=ifelse(dataB$ofstatus_m13==1, wtPCases, wtPControls))
  
  # a single VE(s) curve
  markerVals <- seq(ifelse(markerName=="S3", 1, 0.8), 4.3, by=0.1)
  VEcurve <- sapply(markerVals, function(s){ pVE(s, data, dataI, fit1, fLM, gLM, markerName) })
  
  oList <- list(VE=VEcurve, markerVals=markerVals)
  
  if (!is.null(saveFile)){
    save(oList, file=file.path(saveDir, saveFile))
    cat("Output saved in:\n", file.path(saveDir, saveFile), "\n\n")
  } else {
    return(oList)
  }  
}

# 'plotVEcurve' plots point and interval estimates of VE(s), where s is a continuous biomarker, together
# with the histogram of s measured in the vaccine arm, with incorporated sampling weights to recover the
# case:control ratio in the target population
# data                  - complete data
# dataI                 - immunogenicity data
# markerName            - one of "AUC", "Min", "S1", "S2", "S3", "S4"
# VEcurveFile           - a character string; an .RData file name storing the oList object returned by VEcurve (point estimate)
# bVEcurveFile          - a character string; an .RData file name storing the bList object returned by bVEcurve (interval estimate)
# hingePoint            - shall the minimum of the two estimated hinge points from the placebo and vaccine models be plotted?
# percentileBootstrap   - shall percentile bootstrap pointwise CIs be plotted?
# simultCI              - shall the simultaneous confidence band be plotted?
# xMax                  - the max value of the x-axis
# showLegend            - shall the legend be plotted?
# leaveOutPropBootRuns  - [ignore]
# saveDir               - a directory path for VEcurveFile and bVEcurveFile
plotVEcurve <- function(data, dataI, markerName, histogramOnly=FALSE, VEcurveFile, bVEcurveFile=NULL, hingePoint=TRUE, percentileBootstrap=FALSE, 
                        simultCI=TRUE, xMax=NULL, showLegend=TRUE, leaveOutPropBootRuns=-1, saveDir, verbose=FALSE){
  markerName2 <- switch(markerName, Min="Min", AUC="AUCMB", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4", verbose=FALSE)
  markerName3 <- switch(markerName, Min="the minimum $\\log_{10}$ serotype-specific titer", AUC="the AUC-MB", S1="$\\log_{10}$ serotype 1 titers", S2="$\\log_{10}$ serotype 2 titers", S3="$\\log_{10}$ serotype 3 titers", S4="$\\log_{10}$ serotype 4 titers")
  markerVar <- paste0("IMPSTLOG.", markerName2)
  fstatusVar <- paste0(switch(markerName, Min="o", AUC="o", S1="s1", S2="s2", S3="s3", S4="s4"),"fstatus_m13") 
  markerX <- switch(markerName, Min="Minimum Log10 Titer", AUC="Average Titer", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  
  
  #xLim <- range(dataI[,markerVar], na.rm=TRUE)
  xLim <- c(0.699,5.172)
  par(mar=c(4,5,4,4.2), las=1, cex.axis=0.9, cex.lab=1)
  
  if (histogramOnly){
    plot(0,0, xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=xLim, ylim=c(-1,1.3), xaxt="n", yaxt="n", lwd=2.5, type="n")
    axis(side=1, at=seq(0.5,5.5,by=0.5))
    text(mean(xLim), 0.9, "Insufficient Data", cex=1)
  } else {
    load(file.path(saveDir, VEcurveFile))  
    x <- oList$markerVals
    y <- oList$VE
    if (!is.null(xMax)){
      y <- y[x<=xMax]
      x <- x[x<=xMax]
    }
    
    if (verbose){
      write(c(markerName,x[1],y[1]), file=file.path(saveDir,"verbose.txt"), append=TRUE)
    }
    
    plot(x, y, type="l", xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=xLim, ylim=c(-1,1.3), xaxt="n", yaxt="n", col="red3", lwd=2.5)
    #axis(side=1, at=seq(0.5,5.5,by=0.5))
    axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
    axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
    mtext("Vaccine Efficacy", side=2, las=0, line=3.6, cex=1.1)
    
    if (!is.null(bVEcurveFile)){
      # a list named 'results' (the parallelized version)
      load(file.path(saveDir, bVEcurveFile))
      bVE <- NULL
      for (i in 1:length(results)){
        if (!is.null(names(results[[i]]))){
          bVE <- rbind(bVE, results[[i]]$bVE)
        }
      }   
      #bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
      
      # this code assembles 'bVE' when the bootstrap .RData list is generated in a non-parallelized fashion
      #bVE <- bList$bVE
      #s <- bList$markerVals
      
      # extract the grid on which bootstrap VE is evaluated
      for (i in 1:length(results)){
        if (!is.null(names(results[[i]]))){
          s <- results[[i]]$markerVals
          break
        }
      }
      
      s <- s[s<=max(x)]
      if (length(s)<NCOL(bVE)){ bVE <- bVE[,-((length(s)+1):NCOL(bVE))] }
      
      oVE <- oList$VE
      #oVE <- ifelse(oVE>0.99,0.99,oVE)
      logRR <- log(1-oVE)
      
      #bVE <- ifelse(bVE>0.99,0.99,bVE)
      bLogRR <- log(1-bVE)
      
      # remove all nonfinite values in bLogRR
      #bLogRR <- ifelse(is.infinite(bLogRR),NA,bLogRR)
      
      # remove bootstrap runs with the largest sum of squared differences in VE
      #     if (leaveOutPropBootRuns>0){
      #       if (markerName=="S3"){ leaveOutRange <- 1:NCOL(bVE) }
      #       if (markerName=="S4"){ leaveOutRange <- 30:33 }
      #       V <- apply(t(t(bVE[,leaveOutRange])-oVE[leaveOutRange])^2, 1, sum, na.rm=TRUE)
      #       idx <- which(V < quantile(V, probs=1-leaveOutPropBootRuns, na.rm=TRUE))
      #       bVE <- bVE[idx,]
      #       bLogRR <- bLogRR[idx,]
      #     }
      
      # bootstrap SE of log RR estimates
      bSE <- apply(bLogRR, 2, sd, na.rm=TRUE)
      
      if (NCOL(bLogRR)!=length(logRR)){ stop("The point and bootstrap VE estimates have mismatching sets of marker values.") }
      
      if (percentileBootstrap){
        if (is.null(xMax)){
          lines(s, apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE), lty="dashed", lwd=2, col="black")
          lines(s, apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE), lty="dashed", lwd=2, col="black")    
        } else {
          lines(s[s<=xMax], apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)[s<=xMax], lty="dashed", lwd=2, col="black")
          lines(s[s<=xMax], apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)[s<=xMax], lty="dashed", lwd=2, col="black")
        }
      } else {
        LB.VE <- 1 - exp(logRR + qnorm(0.975) * bSE)
        UB.VE <- 1 - exp(logRR - qnorm(0.975) * bSE)
        if (is.null(xMax)){
          lines(s, LB.VE, lty="dashed", lwd=2)
          lines(s, UB.VE, lty="dashed", lwd=2)  
          
          if (verbose){
            write(c(s[1],LB.VE[1],UB.VE[1]), file=file.path(saveDir,"verbose.txt"), append=TRUE)
          }
        } else {
          lines(s[s<=xMax], LB.VE[s<=xMax], lty="dashed", lwd=2)
          lines(s[s<=xMax], UB.VE[s<=xMax], lty="dashed", lwd=2)
          
          if (verbose){
            write(c(s[1],LB.VE[1],UB.VE[1]), file=file.path(saveDir,"verbose.txt"), append=TRUE)
          }
        }
      }
      
      if (simultCI){
        supAbsZ <- NULL
        for (i in 1:NROW(bLogRR)){ 
          Z <- abs((bLogRR[i,]-logRR)/bSE)
          supAbsZ <- c(supAbsZ, max(Z, na.rm=!all(is.na(Z)))) 
        }
        qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
        
        LB.VE <- 1 - exp(logRR + qSupAbsZ * bSE)
        UB.VE <- 1 - exp(logRR - qSupAbsZ * bSE)
        if (is.null(xMax)){
          lines(s, LB.VE, lty="dotted", lwd=3)
          lines(s, UB.VE, lty="dotted", lwd=3)  
        } else {
          lines(s[s<=xMax], LB.VE[s<=xMax], lty="dotted", lwd=3)
          lines(s[s<=xMax], UB.VE[s<=xMax], lty="dotted", lwd=3)  
        }
      }
    }
  }
  
  colHist <- c(col2rgb("olivedrab3"))
  colHist <- rgb(colHist[1], colHist[2], colHist[3], alpha=255*0.4, maxColorValue=255)
  xLeg1 <- ifelse(markerName=="S3", 3.2, ifelse(markerName=="S4", 2.7, 0.7))
  xLeg2 <- 3.5
  yLeg1 <- ifelse(markerName %in% c("S3","S4"), 0.3, 1.3)
  yLeg2 <- ifelse(markerName %in% c("S3","S4"), 0.1, 1.1)
  if (showLegend){
    legend(xLeg1, yLeg1, lty=c("dashed","dotted"), lwd=c(2,3), col="black", legend=paste0("95% ",c("Pointwise","Simultaneous")," CI"), bty="n", cex=0.7, y.intersp=0.8)
    #legend(xLeg, yLeg1, lty="dashed", lwd=2, col="black", legend="95% Pointwise CI", bty="n", cex=0.7, y.intersp=0.8)
    legend(xLeg1, yLeg2, fill=colHist, border=colHist, legend="Vaccine Group", bty="n", cex=0.7)    
  }
  if (hingePoint){ text(xLim[2], 1.17, paste0("Hinge Point = ",round(10^min(oList$cpointP, oList$cpointV, na.rm=TRUE),0)), cex=0.7, pos=2) }
  
  # sampling weights in hist
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  yLim <- c(0, max(hist(S.V, plot=FALSE, breaks=14)$density)*2.5)
  par(new=TRUE)
  hist(S.V, col=colHist, axes=FALSE, labels=FALSE, main="", xlab="", ylab="", border=0, freq=FALSE, xlim=xLim, ylim=yLim, breaks=ifelse(markerName=="S3",16,14))
  if (markerName %in% c("Min","AUC","S4")){ ticks <- seq(0,1,by=0.25) }
  if (markerName %in% c("S2","S3")){ ticks <- seq(0,0.75,by=0.25) }
  if (markerName=="S1"){ ticks <- seq(0,0.5,by=0.25) }
  axis(side=4, at=ticks)
  mtext("Density                          ", side=4, las=0, line=3, cex=1.1)  
  
  epoint <- switch(markerName, Min="primary", AUC="primary", S1="matched serotype 1", S2="matched serotype 2", S3="matched serotype 3", S4="matched serotype 4")  
  #cat("\\caption[Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint with the 95\\% pointwise confidence band]{Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint with the 95\\% pointwise confidence band based on $10^3$ bootstrap iterations and the two phase-adjusted histogram of ",markerName3," measured in the vaccine group at month 13 in the CYD15 trial}\\label{Fig: VE(s1) ",markerName2,"}", sep="")
}

plotVEcurveWithLimitedBaselineMarkersInVaccinees <- function(data, dataI, markerName, histogramOnly=FALSE, VEcurveFile, VEcurveFileMC, hingePoint=TRUE, xMax=NULL, 
                                                             showLegend=TRUE, VEcurveDir, VEcurveMCdir, main){
  markerName2 <- switch(markerName, Min="Min", AUC="AUCMB", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, Min="the minimum $\\log_{10}$ serotype-specific titer", AUC="the AUC-MB", S1="$\\log_{10}$ serotype 1 titers", S2="$\\log_{10}$ serotype 2 titers", S3="$\\log_{10}$ serotype 3 titers", S4="$\\log_{10}$ serotype 4 titers")
  markerVar <- paste0("IMPSTLOG.", markerName2)
  fstatusVar <- paste0(switch(markerName, Min="o", AUC="o", S1="s1", S2="s2", S3="s3", S4="s4"),"fstatus_m13") 
  markerX <- switch(markerName, Min="Minimum Log10 Titer", AUC="Average Titer", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  
  par(mar=c(4.5,5,2.5,1), las=1, cex.axis=0.9, cex.lab=1)
  
  if (histogramOnly){
    plot(0,0, xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=xLim, ylim=c(-1,1.3), xaxt="n", yaxt="n", lwd=2.5, type="n")
    axis(side=1, at=seq(0.5,5.5,by=0.5))
    text(mean(xLim), 0.9, "Insufficient Data", cex=1)
  } else {
    load(file.path(VEcurveMCdir, VEcurveFileMC))
    
    plot(oList$markerVals, oList$VE[1,], xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=range(oList$markerVals[oList$markerVals<=4.4]), ylim=c(0,1), xaxt="n", yaxt="n", type="n",
         main=main)
    axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
    axis(side=2, at=seq(-0.5,1,by=0.25), labels=seq(-50,100,by=25))
    mtext("Vaccine Efficacy (%)", side=2, las=0, line=3, cex=1.1)
    
    for (b in 1:NROW(oList$VE)){ lines(oList$markerVals[oList$markerVals<=4.4], oList$VE[b,][oList$markerVals<=4.4], lwd=1.5, col="red") }
    
    load(file.path(VEcurveDir, VEcurveFile))
    
    lines(oList$markerVals[oList$markerVals<=4.4], oList$VE[oList$markerVals<=4.4], lwd=2.5)
  }
  
  colHist <- c(col2rgb("olivedrab3"))
  colHist <- rgb(colHist[1], colHist[2], colHist[3], alpha=255*0.4, maxColorValue=255)
  xLeg1 <- ifelse(markerName=="S3", 3.2, ifelse(markerName=="S4", 2.7, 0.7))
  xLeg2 <- 3.5
  yLeg1 <- ifelse(markerName %in% c("S3","S4"), 0.3, 1.3)
  yLeg2 <- ifelse(markerName %in% c("S3","S4"), 0.1, 1.1)
  if (showLegend){
    #legend(xLeg1, yLeg1, lty=c("dashed","dotted"), lwd=c(2,3), col="black", legend=paste0("95% ",c("Pointwise","Simultaneous")," CI"), bty="n", cex=0.7, y.intersp=0.8)
    #legend(xLeg, yLeg1, lty="dashed", lwd=2, col="black", legend="95% Pointwise CI", bty="n", cex=0.7, y.intersp=0.8)
    legend(xLeg1, 0.9, fill=colHist, border=colHist, legend="Vaccine group", bty="n", cex=0.7) 
    legend(xLeg1, 1, lwd=1.5, col="red", legend="Monte-Carlo estimates", bty="n", cex=0.7)
  }
  #if (hingePoint){ text(xLim[2], 1.17, paste0("Hinge Point = ",round(10^min(oList$cpointP, oList$cpointV, na.rm=TRUE),0)), cex=0.7, pos=2) }
  
  # sampling weights in hist
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  yLim <- c(0, max(hist(S.V, plot=FALSE, breaks=14)$density)*2.5)
  par(new=TRUE)
  hist(S.V, col=colHist, axes=FALSE, labels=FALSE, main="", xlab="", ylab="", border=0, freq=FALSE, xlim=range(oList$markerVals[oList$markerVals<=4.4]), ylim=yLim, breaks=ifelse(markerName=="S3",16,14))
  if (markerName %in% c("Min","AUC","S4")){ ticks <- seq(0,1,by=0.25) }
  if (markerName %in% c("S2","S3")){ ticks <- seq(0,0.75,by=0.25) }
  if (markerName=="S1"){ ticks <- seq(0,0.5,by=0.25) }
  #   axis(side=4, at=ticks)
  #   mtext("Density                          ", side=4, las=0, line=3, cex=1.1)  
}


plotVEcurve_vLondonTalk <- function(data, dataI, markerName, histogramOnly=FALSE, VEcurveFile, bVEcurveFile=NULL, hingePoint=TRUE, percentileBootstrap=FALSE, simultCI=TRUE, xMax=NULL, showLegend=TRUE, leaveOutPropBootRuns=-1, saveDir){
  markerName2 <- switch(markerName, Min="Min", AUC="AUCMB", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, Min="the minimum $\\log_{10}$ serotype-specific titer", AUC="the AUC-MB", S1="$\\log_{10}$ serotype 1 titers", S2="$\\log_{10}$ serotype 2 titers", S3="$\\log_{10}$ serotype 3 titers", S4="$\\log_{10}$ serotype 4 titers")
  markerVar <- paste0("IMPSTLOG.", markerName2)
  fstatusVar <- paste0(switch(markerName, Min="o", AUC="o", S1="s1", S2="s2", S3="s3", S4="s4"),"fstatus_m13") 
  markerX <- switch(markerName, Min="Minimum Log10 Titer", AUC="Average Titer", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  
  
  #xLim <- range(dataI[,markerVar], na.rm=TRUE)
  xLim <- c(0.699,5.172)
  par(mar=c(4.3,3.5,1,1), cex.axis=0.9, cex.lab=1)
  
  if (histogramOnly){
    plot(0,0, xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=xLim, ylim=c(-1,1.3), xaxt="n", yaxt="n", lwd=2.5, type="n")
    axis(side=1, at=seq(0.5,5.5,by=0.5))
    text(mean(xLim), 0.9, "Insufficient Data", cex=1)
  } else {
    load(file.path(saveDir, VEcurveFile))  
    x <- oList$markerVals
    y <- oList$VE
    if (!is.null(xMax)){
      y <- y[x<=xMax]
      x <- x[x<=xMax]
    }
    plot(x, y, type="l", xlab=expression(paste("Month 13 ",PRNT[50])), ylab="", xlim=xLim, ylim=0:1, xaxt="n", yaxt="n", lwd=2.5)
    #axis(side=1, at=seq(0.5,5.5,by=0.5))
    axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
    axis(side=2, at=seq(0,1,by=0.25), labels=paste0(seq(0,100,by=25),"%"))
    mtext("Estimated VE(s)", side=2, las=0, line=2, cex=1.1)
    
    if (!is.null(bVEcurveFile)){
      # a list named 'results' (the parallelized version)
      load(file.path(saveDir, bVEcurveFile))
      bVE <- NULL
      for (i in 1:length(results)){
        if (!is.null(names(results[[i]]))){
          bVE <- rbind(bVE, results[[i]]$bVE)
        }
      }   
      #bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
      
      # this code assembles 'bVE' when the bootstrap .RData list is generated in a non-parallelized fashion
      #bVE <- bList$bVE
      #s <- bList$markerVals
      
      # extract the grid on which bootstrap VE is evaluated
      for (i in 1:length(results)){
        if (!is.null(names(results[[i]]))){
          s <- results[[i]]$markerVals
          break
        }
      }
      
      s <- s[s<=max(x)]
      if (length(s)<NCOL(bVE)){ bVE <- bVE[,-((length(s)+1):NCOL(bVE))] }
      
      oVE <- oList$VE
      #oVE <- ifelse(oVE>0.99,0.99,oVE)
      logRR <- log(1-oVE)
      
      #bVE <- ifelse(bVE>0.99,0.99,bVE)
      bLogRR <- log(1-bVE)
      
      # remove all nonfinite values in bLogRR
      #bLogRR <- ifelse(is.infinite(bLogRR),NA,bLogRR)
      
      # remove bootstrap runs with the largest sum of squared differences in VE
      #     if (leaveOutPropBootRuns>0){
      #       if (markerName=="S3"){ leaveOutRange <- 1:NCOL(bVE) }
      #       if (markerName=="S4"){ leaveOutRange <- 30:33 }
      #       V <- apply(t(t(bVE[,leaveOutRange])-oVE[leaveOutRange])^2, 1, sum, na.rm=TRUE)
      #       idx <- which(V < quantile(V, probs=1-leaveOutPropBootRuns, na.rm=TRUE))
      #       bVE <- bVE[idx,]
      #       bLogRR <- bLogRR[idx,]
      #     }
      
      # bootstrap SE of log RR estimates
      bSE <- apply(bLogRR, 2, sd, na.rm=TRUE)
      
      if (NCOL(bLogRR)!=length(logRR)){ stop("The point and bootstrap VE estimates have mismatching sets of marker values.") }
      
      if (percentileBootstrap){
        if (is.null(xMax)){
          lines(s, apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE), lty="dashed", lwd=2, col="black")
          lines(s, apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE), lty="dashed", lwd=2, col="black")    
        } else {
          lines(s[s<=xMax], apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)[s<=xMax], lty="dashed", lwd=2, col="black")
          lines(s[s<=xMax], apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)[s<=xMax], lty="dashed", lwd=2, col="black")
        }
      } else {
        LB.VE <- 1 - exp(logRR + qnorm(0.975) * bSE)
        UB.VE <- 1 - exp(logRR - qnorm(0.975) * bSE)
        if (is.null(xMax)){
          lines(s, LB.VE, lty="dashed", lwd=2)
          lines(s, UB.VE, lty="dashed", lwd=2)  
        } else {
          lines(s[s<=xMax], LB.VE[s<=xMax], lty="dashed", lwd=2)
          lines(s[s<=xMax], UB.VE[s<=xMax], lty="dashed", lwd=2)
        }
      }
      
      if (simultCI){
        supAbsZ <- NULL
        for (i in 1:NROW(bLogRR)){ 
          Z <- abs((bLogRR[i,]-logRR)/bSE)
          supAbsZ <- c(supAbsZ, max(Z, na.rm=!all(is.na(Z)))) 
        }
        qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
        
        LB.VE <- 1 - exp(logRR + qSupAbsZ * bSE)
        UB.VE <- 1 - exp(logRR - qSupAbsZ * bSE)
        if (is.null(xMax)){
          lines(s, LB.VE, lty="dotted", lwd=3)
          lines(s, UB.VE, lty="dotted", lwd=3)  
        } else {
          lines(s[s<=xMax], LB.VE[s<=xMax], lty="dotted", lwd=3)
          lines(s[s<=xMax], UB.VE[s<=xMax], lty="dotted", lwd=3)  
        }
      }
    }
  }
}


plotBvsM13inP <- function(data, markerName){
  require(robustbase)
  markerName2 <- switch(markerName, AUC="AUCMB", Min="Min", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, AUC="AUC-MB", Min="minimum titer", S1="$\\log_{10}$ serotype 1 titer", S2="$\\log_{10}$ serotype 2 titer", S3="$\\log_{10}$ serotype 3 titer", S4="$\\log_{10}$ serotype 4 titer")
  markerX <- switch(markerName, AUC="AUC-MB", Min="Minimum Titer", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  if (markerName!="AUC"){
    if (markerName=="Min"){  # assumes that there exist variables 'bMin' and 'IMPSTLOG.Min'
      data$bAUC <- data$bMin
      data$IMPSTLOG.AUCMB <- data$IMPSTLOG.Min
    } else {
      data$ofstatus_m0 <- data[,paste0(tolower(markerName),"fstatus_m0")]
      data$ofstatus_m13 <- data[,paste0(tolower(markerName),"fstatus_m13")]
      data$oftime_m13 <- data[,paste0(tolower(markerName),"ftime_m13")]
      data$bAUC <- data[,paste0("b",markerName)]
      data$IMPSTLOG.AUCMB <- data[,paste0("IMPSTLOG.Sero",substr(markerName, start=2, stop=2))] 
    }       
  }
  
  dataI <- subset(data, !is.na(IMPSTLOG.AUCMB))
  dataP <- subset(dataI, VACC==0 & !is.na(bAUC)) # subset of the placebo group in the imm set with baseline markers
  par(mar=c(4,5,3,4.2), las=1, cex.axis=0.9, cex.lab=1)
  with(dataP, plot(bAUC, IMPSTLOG.AUCMB, xlab=paste0("Baseline ",markerX," in Placebo Recipients"),
                   ylab=paste0("Month 13 ",markerX,ifelse(markerName=="AUC","","\n")," in Placebo Recipients"), 
                   cex=0.7))
  abline(0,1, lty="dotted", lwd=2)
  abline(lmrob(IMPSTLOG.AUCMB ~ bAUC, data=dataP, maxit.scale=500), lwd=2.5, col="darkgoldenrod2")
  with(dataP, lines(lowess(bAUC, IMPSTLOG.AUCMB), col="red", lwd=2, lty="dashed"))
  legend("topleft", lwd=2, lty=c("solid","dashed","dotted"), col=c("darkgoldenrod2", "red", "black"),
         legend=c("Robust linear reg (Yohai, 1987)", "LOWESS", "y=x"), bty="n", cex=0.7)
  legend("bottomright", legend=paste0("Spearman's r=",with(dataP, round(cor(bAUC, IMPSTLOG.AUCMB, method="spearman"), 3))), 
         cex=0.8, bty="n")
  #cat("\\caption{Association between $S(0)$ and $S_b$ measurements of ",markerName3," in arm $Z=0$ in the CYD14 trial}\\label{Fig: Sbase vs S0 ",markerName2,"}", sep="")
}

# the 'plotVEcurve' version for the age-stratified analysis
# 'stratum'   - one of "le7" (less than or equal to 7 years) or "g7" (greater than 7 years) because 
#               'bVEcurveFile' includes both
plotVEcurveA <- function(data, dataI, markerName, VEcurveFile, bVEcurveFile, stratum,
                         pVEcurveFile=NULL, VEcurveMin=0.8, legPos="left", saveDir){
  markerName2 <- switch(markerName, AUC="AUCMB", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, AUC="the AUC-MB", S1="$\\log_{10}$ serotype 1 titers", S2="$\\log_{10}$ serotype 2 titers", S3="$\\log_{10}$ serotype 3 titers", S4="$\\log_{10}$ serotype 4 titers")
  markerVar <- paste0("IMPSTLOG.", markerName2)
  fstatusVar <- paste0(switch(markerName, AUC="o", S1="s1", S2="s2", S3="s3", S4="s4"),"fstatus_m13") 
  markerX <- switch(markerName, AUC="AUC-MB", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  
  load(file.path(saveDir, VEcurveFile))
  
  if (markerName=="S1"){ xLim <- c(min(dataI[,markerVar], na.rm=TRUE) - 0.1, max(4.5, max(dataI[,markerVar], na.rm=TRUE))) }
  if (markerName=="S2"){ xLim <- c(min(dataI[,markerVar], na.rm=TRUE), min(4.7, max(dataI[,markerVar], na.rm=TRUE))) }
  if (markerName %in% c("S3","S4","AUC")){ xLim <- c(min(dataI[,markerVar], na.rm=TRUE), max(4.5, max(dataI[,markerVar], na.rm=TRUE))) } 
  par(mar=c(4,5,3,4.2), las=1, cex.axis=0.9, cex.lab=1)
  plot(oList$markerVals[oList$markerVals>=VEcurveMin], oList$VE[oList$markerVals>=VEcurveMin], type="l", 
       xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=xLim, ylim=c(-1,1), xaxt="n", yaxt="n", 
       col="red3", lwd=2)
  axis(side=1, at=seq(1,4.5,by=0.5))
  axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
  mtext(expression(paste("Estimated VE(",s[1],")"), sep=""), side=2, las=0, line=3.6, cex=1.1)
  
  load(file.path(saveDir, bVEcurveFile))
  if (stratum=="le7"){
    results <- resultsle7
  } else {
    results <- resultsg7
  }
  if (markerName=="S3"){
    bVE <- NULL
    for (i in 1:length(results)){
      if (!is.null(names(results[[i]]))){
        bVE <- rbind(bVE, results[[i]]$bVE)
      }
    }
  } else {
    bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
  }
  s <- results[[1]]$markerVals
  lines(s[s>=VEcurveMin], apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)[s>=VEcurveMin], lty="dashed", lwd=2, col="black")
  lines(s[s>=VEcurveMin], apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)[s>=VEcurveMin], lty="dashed", lwd=2, col="black")
  
  xLeg <- ifelse(legPos=="right", 3.1, 0.7)
  yLeg1 <- ifelse(legPos=="right", 0.4, 1)
  yLeg2 <- ifelse(legPos=="right", 0.3, 0.9)
  legend(xLeg, yLeg1, lty="dashed", lwd=2, col="black", legend="95% pointwise CI", bty="n", cex=0.7)
  legend(xLeg, yLeg2, fill=rgb(0,1,0,alpha=0.23), legend=paste0("Vaccine ",ifelse(stratum=="le7","2-7","8-14")," years"), bty="n", cex=0.7)  
  
  if (!is.null(pVEcurveFile)){
    load(file.path(saveDir, pVEcurveFile))
    lines(oList$markerVals[oList$markerVals>=VEcurveMin], oList$VE[oList$markerVals>=VEcurveMin], lwd=2,
          lty="dotdash", col="orangered")
  }
  
  # sampling weights in hist
  if (stratum=="le7"){
    data <- subset(data, AGEYRS<=7)
    dataI <- subset(dataI, AGEYRS<=7)
  } else {
    data <- subset(data, AGEYRS>7)
    dataI <- subset(dataI, AGEYRS>7)
  }
  
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  ylim <- c(0, max(hist(S.V, plot=FALSE, breaks=14)$density)+ifelse(markerName=="S1",0.05,0))
  par(new=TRUE)
  hist(S.V, col=rgb(0,1,0,alpha=0.23), axes=FALSE, labels=FALSE, 
       main="", xlab="", ylab="", border=0, freq=FALSE, xlim=xLim, ylim=ylim, breaks=ifelse(markerName=="S1",18,14))
  axis(side=4, at=axTicks(side=2))
  #par(new=TRUE)
  #hist(subset(dataI, VACC==0 & IMPSTLOG.AUCMB<=4)$IMPSTLOG.AUCMB, col=rgb(0,0,1,alpha=0.22), axes=FALSE, labels=FALSE, main="",xlab="", ylab="", border=0, freq=TRUE, ylim=ylim, breaks=10)
  mtext("Density", side=4, las=0, line=3, cex=1.1)
  
  epoint <- switch(markerName, AUC="overall", S1="serotype 1", S2="serotype 2", S3="serotype 3", S4="serotype 4")  
  if (stratum=="le7"){
    cat("\\caption{Estimated $VE(s_1)$ against ",epoint, " dengue \\emph{among 2--7-year-olds} with the 95\\% pointwise confidence band based on $10^3$ bootstrap iterations and the two phase-adjusted histogram of ",markerName3," measured in the vaccine group at month 13 in the CYD14 trial}\\label{Fig: VE(s1) ",markerName2," ",stratum,"}", sep="")
  } else {
    cat("\\caption{Estimated $VE(s_1)$ against ",epoint, " dengue \\emph{among 8--14-year-olds} with the 95\\% pointwise confidence band based on $10^3$ bootstrap iterations and the two phase-adjusted histogram of ",markerName3," measured in the vaccine group at month 13 in the CYD14 trial}\\label{Fig: VE(s1) ",markerName2," ",stratum,"}", sep="")
  }
  
}

plotVEcurveOverlayA <- function(dataFull, dataFullI, markerName, VEcurveFileAle11, VEcurveFileAg11, bVEcurveFileAle11, bVEcurveFileAg11, legOut=FALSE, saveDir){
  markerName2 <- switch(markerName, AUC="AUCMB", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, AUC="the AUC-MB", S1="$\\log_{10}$ serotype 1 titers", S2="$\\log_{10}$ serotype 2 titers", S3="$\\log_{10}$ serotype 3 titers", S4="$\\log_{10}$ serotype 4 titers")
  markerVar <- paste0("IMPSTLOG.", markerName2)
  fstatusVar <- paste0(switch(markerName, AUC="o", S1="s1", S2="s2", S3="s3", S4="s4"),"fstatus_m13") 
  markerX <- switch(markerName, AUC="AUC-MB", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  
  # initiate the plot
  load(file.path(saveDir, VEcurveFileAle11))
  xLim <- range(dataFullI[,markerVar], na.rm=TRUE)
  
  par(mar=c(4+ifelse(legOut,3.5,0),5,3,4.2), las=1, cex.axis=0.9, cex.lab=1)
  plot(oList$markerVals, oList$VE, type="n", xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=xLim, ylim=c(-1,1), xaxt="n", yaxt="n")
  axis(side=1, at=seq(0.5,5.5,by=0.5))
  axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
  mtext(expression(paste("Estimated VE(",s[1],")"), sep=""), side=2, las=0, line=3.6, cex=1.1)
  
  colRGBle11 <- c(col2rgb("red1"))
  colRGBle11 <- rgb(colRGBle11[1], colRGBle11[2], colRGBle11[3], alpha=255*0.25, maxColorValue=255)
  colRGBg11 <- c(col2rgb("gray40"))
  colRGBg11 <- rgb(colRGBg11[1], colRGBg11[2], colRGBg11[3], alpha=255*0.5, maxColorValue=255)
  
  # interval estimates of VE(s_1)
  if (!is.null(bVEcurveFileAle11)){
    load(file.path(saveDir, bVEcurveFileAle11))  
    bVE <- NULL
    for (i in 1:length(results)){
      if (!is.null(names(results[[i]]))){
        bVE <- rbind(bVE, results[[i]]$bVE)
      }
    }    
    #bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
    
    s <- results[[1]]$markerVals
    UB <- apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)
    LB <- apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)
    idxNA <- which(is.na(UB))
    s <- s[-idxNA]
    UB <- UB[-idxNA]
    LB <- LB[-idxNA]
    xCoord <- c(s, rev(s))
    yCoord <- c(UB, rev(LB))    
    polygon(xCoord, yCoord, col=colRGBle11, border=NA)
  }  
  
  if (!is.null(bVEcurveFileAg11)){  
    load(file.path(saveDir, bVEcurveFileAg11))
    bVE <- NULL
    for (i in 1:length(results)){
      if (!is.null(names(results[[i]]))){
        bVE <- rbind(bVE, results[[i]]$bVE)
      }
    }    
    #bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
    
    s <- results[[1]]$markerVals
    UB <- apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)
    LB <- apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)
    idxNA <- which(is.na(UB))
    if (length(idxNA)>0){
      s <- s[-idxNA]
      UB <- UB[-idxNA]
      LB <- LB[-idxNA]
    }    
    xCoord <- c(s, rev(s))
    yCoord <- c(UB, rev(LB))
    polygon(xCoord, yCoord, col=colRGBg11, border=NA)    
  }
  
  # point estimates of VE(s_1)
  load(file.path(saveDir, VEcurveFileAle11))
  lines(oList$markerVals, oList$VE, col="red3", lwd=2.5, lty="longdash")
  load(file.path(saveDir, VEcurveFileAg11))
  lines(oList$markerVals, oList$VE, col="black", lwd=2.5)
  
  # sampling weights in hist
  data <- subset(dataFull, AGEYRS<=11)
  dataI <- subset(dataFullI, AGEYRS<=11)  
  
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  yLim <- c(0, max(hist(S.V, plot=FALSE, breaks=14)$density)*2.5)
  colRGB1 <- c(col2rgb("orange"))
  colRGB1 <- rgb(colRGB1[1], colRGB1[2], colRGB1[3], alpha=255*0.23, maxColorValue=255)
  par(new=TRUE)
  hist(S.V, col=colRGB1, axes=FALSE, labels=FALSE, main="", xlab="", ylab="", border=0, freq=FALSE, xlim=xLim, ylim=yLim, breaks=ifelse(markerName=="S3",16,14))
  if (markerName %in% c("AUC","S4")){ ticks <- seq(0,1,by=0.25) }
  if (markerName %in% c("S2","S3")){ ticks <- seq(0,0.75,by=0.25) }
  if (markerName=="S1"){ ticks <- seq(0,0.5,by=0.25) }
  axis(side=4, at=ticks)
  mtext("Density                          ", side=4, las=0, line=3, cex=1.1)
  
  data <- subset(dataFull, AGEYRS>11)
  dataI <- subset(dataFullI, AGEYRS>11)
  
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  colRGB2 <- c(col2rgb("olivedrab3"))
  colRGB2 <- rgb(colRGB2[1], colRGB2[2], colRGB2[3], alpha=255*0.4, maxColorValue=255)
  par(new=TRUE)
  hist(S.V, col=colRGB2, axes=FALSE, labels=FALSE, main="", xlab="", ylab="", border=0, xlim=xLim, ylim=yLim, freq=FALSE, breaks=ifelse(markerName=="S3",16,14))
  
  if (legOut){
    legend("bottomleft", lty=c("longdash","solid"), lwd=2.5, col=c("red3","black"), legend=c("9-11 years","12-16 years"), y.intersp=0.9, inset=c(0.07,-0.55), xpd=NA, bty="n", cex=0.7)
    legend("bottom", fill=c(colRGBle11,colRGBg11), border=c(colRGBle11,colRGBg11), legend=c("\n95% pointwise CIs  ",""), bty="n", y.intersp=0.3, inset=c(0,-0.533), xpd=NA, cex=0.7)
    legend("bottomright", fill=c(colRGB1, colRGB2), border=c(colRGB1, colRGB2), legend=c("Vaccine 9-11 years","Vaccine 12-16 years"), bty="n", y.intersp=0.9, inset=c(0.02,-0.55), xpd=NA, cex=0.7)
  } else {
    xLeg <- ifelse(markerName=="S3",3.5,3)
    y1Leg <- switch(markerName, AUC=1.65, S1=0.75, S2=1.1, S3=0.95, S4=1.62)
    y2Leg <- switch(markerName, AUC=1.45, S1=0.67, S2=1, S3=0.85, S4=1.4)
    y3Leg <- switch(markerName, AUC=1.1, S1=0.52, S2=0.8, S3=0.65, S4=1.05)
    legend(xLeg, y1Leg, lty=c("longdash","solid"), lwd=2.5, col=c("red3","black"), legend=c("9-11 years","12-16 years"), bty="n", cex=0.7, y.intersp=0.8)
    legend(xLeg, y2Leg, fill=c(colRGBle11,colRGBg11), border=c(colRGBle11,colRGBg11), legend=c("\n95% pointwise CIs",""), bty="n", cex=0.7, y.intersp=0.3)
    legend(xLeg, y3Leg, fill=c(colRGB1, colRGB2), border=c(colRGB1, colRGB2), legend=c("Vaccine 9-11 years","Vaccine 12-16 years"), bty="n", cex=0.7, y.intersp=0.8)
  }  
  
  epoint <- switch(markerName, AUC="overall", S1="serotype 1", S2="serotype 2", S3="serotype 3", S4="serotype 4")
  if (legOut){ cat("\\vspace*{-3mm}") }
  cat("\\caption[Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint comparing subgroups of 9--11 versus 12--16-year-olds, with the 95\\% pointwise confidence bands]{Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint comparing subgroups of 9--11 versus 12--16-year-olds, with the 95\\% pointwise confidence bands based on $10^3$ bootstrap iterations and the two phase-adjusted histograms of ",markerName3," measured in the vaccine group at month 13 in the CYD15 trial}\\label{Fig: VE(s1) ",markerName2,"_ageStrata}", sep="")
}

plotVEcurveOverlayA2 <- function(dataFull, dataFullI, markerName, VEcurveFileAle11, VEcurveFileAg11, bVEcurveFileAle11, bVEcurveFileAg11, legOut=FALSE, saveDir){
  markerName2 <- switch(markerName, Min="Min", AUC="AUCMB", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, Min="the minimum $\\log_{10}$ serotype-specific titer", AUC="the AUC-MB", S1="$\\log_{10}$ serotype 1 titers", S2="$\\log_{10}$ serotype 2 titers", S3="$\\log_{10}$ serotype 3 titers", S4="$\\log_{10}$ serotype 4 titers")
  markerVar <- paste0("IMPSTLOG.", markerName2)
  fstatusVar <- paste0(switch(markerName, Min="o", AUC="o", S1="s1", S2="s2", S3="s3", S4="s4"),"fstatus_m13") 
  markerX <- switch(markerName, Min="Min Log10 Titer", AUC="AUC-MB", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  
  # initiate the plot
  load(file.path(saveDir, VEcurveFileAle11))
  if (markerName=="S1"){ sLimAle11 <- oList$markerVals[c(4,37)] }
  load(file.path(saveDir, VEcurveFileAg11))
  if (markerName=="S1"){ sLimAg11 <- oList$markerVals[c(1,38)] }
  
  xLim <- range(dataFullI[,markerVar], na.rm=TRUE)
  
  par(mar=c(4+ifelse(legOut,3.5,0),5,3,4.2), las=1, cex.axis=0.9, cex.lab=1)
  plot(oList$markerVals, oList$VE, type="n", xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=xLim, ylim=c(-1,1), xaxt="n", yaxt="n")
  axis(side=1, at=seq(0.5,5.5,by=0.5))
  axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
  mtext(expression(paste("Estimated VE(",s[1],")"), sep=""), side=2, las=0, line=3.6, cex=1.1)
  
  colRGBle11 <- c(col2rgb("red1"))
  colRGBle11 <- rgb(colRGBle11[1], colRGBle11[2], colRGBle11[3], alpha=255*0.25, maxColorValue=255)
  colRGBg11 <- c(col2rgb("gray40"))
  colRGBg11 <- rgb(colRGBg11[1], colRGBg11[2], colRGBg11[3], alpha=255*0.5, maxColorValue=255)
  
  if (!is.null(bVEcurveFileAle11)){
    if (markerName=="Min"){
      load(file.path(saveDir, bVEcurveFileAle11))
      bVE <- NULL
      for (i in 1:length(results)){
        if (!is.null(names(results[[i]]))){
          bVE <- rbind(bVE, results[[i]]$bVE)
        }
      }   
      
      s <- results[[1]]$markerVals
      LB <- apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)
      UB <- apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)
      idxNA <- which(is.na(UB))
      if (length(idxNA)>0){
        s <- s[-idxNA]
        UB <- UB[-idxNA]
        LB <- LB[-idxNA]
      }
    } else {
      # interval estimates of VE(s_1)      
      load(file.path(saveDir, bVEcurveFileAle11))
      
      s <- bList$s
      UB <- bList$UB
      LB <- bList$LB
      idxNA <- which(is.na(UB))
      if (length(idxNA)>0){
        s <- s[-idxNA]
        UB <- UB[-idxNA]
        LB <- LB[-idxNA]
      }
      if (markerName=="S1"){
        keep <- which(s>=sLimAle11[1] & s<=sLimAle11[2])
        s <- s[keep]
        UB <- UB[keep]
        LB <- LB[keep]
      }
    }
    if (markerName=="S4"){
      xCoord <- c(s[1:2], rev(s[1:2]))
      yCoord <- c(UB[1:2], rev(LB[1:2]))    
      polygon(xCoord, yCoord, col=colRGBle11, border=NA)
      xCoord <- c(s[3:34], rev(s[3:34]))
      yCoord <- c(UB[3:34], rev(LB[3:34]))    
      polygon(xCoord, yCoord, col=colRGBle11, border=NA)
    } else {
      xCoord <- c(s, rev(s))
      yCoord <- c(UB, rev(LB))    
      polygon(xCoord, yCoord, col=colRGBle11, border=NA)
    }    
  }  
  
  if (!is.null(bVEcurveFileAg11)){
    if (markerName=="Min"){
      load(file.path(saveDir, bVEcurveFileAg11))
      bVE <- NULL
      for (i in 1:length(results)){
        if (!is.null(names(results[[i]]))){
          bVE <- rbind(bVE, results[[i]]$bVE)
        }
      }   
      
      s <- results[[1]]$markerVals
      LB <- apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)
      UB <- apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)
      idxNA <- which(is.na(UB))
      if (length(idxNA)>0){
        s <- s[-idxNA]
        UB <- UB[-idxNA]
        LB <- LB[-idxNA]
      }
    } else {
      load(file.path(saveDir, bVEcurveFileAg11))
      
      s <- bList$s
      UB <- bList$UB
      LB <- bList$LB
      idxNA <- which(is.na(UB))
      if (length(idxNA)>0){
        s <- s[-idxNA]
        UB <- UB[-idxNA]
        LB <- LB[-idxNA]
      }
    }    
    if (markerName=="S1"){
      keep <- which(s>=sLimAg11[1] & s<=sLimAg11[2])
      s <- s[keep]
      UB <- UB[keep]
      LB <- LB[keep]
    }    
    xCoord <- c(s, rev(s))
    yCoord <- c(UB, rev(LB))
    polygon(xCoord, yCoord, col=colRGBg11, border=NA)    
  }
  
  # point estimates of VE(s_1)
  load(file.path(saveDir, VEcurveFileAle11))
  lines(oList$markerVals[!is.na(oList$VE)], na.omit(oList$VE), col="red3", lwd=2.5, lty="longdash")
  load(file.path(saveDir, VEcurveFileAg11))
  lines(oList$markerVals[!is.na(oList$VE)], na.omit(oList$VE), col="black", lwd=2.5)
  
  # sampling weights in hist
  data <- subset(dataFull, AGEYRS<=11)
  dataI <- subset(dataFullI, AGEYRS<=11)  
  
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  yLim <- c(0, max(hist(S.V, plot=FALSE, breaks=14)$density)*2.5)
  colRGB1 <- c(col2rgb("orange"))
  colRGB1 <- rgb(colRGB1[1], colRGB1[2], colRGB1[3], alpha=255*0.23, maxColorValue=255)
  par(new=TRUE)
  hist(S.V, col=colRGB1, axes=FALSE, labels=FALSE, main="", xlab="", ylab="", border=0, freq=FALSE, xlim=xLim, ylim=yLim, breaks=ifelse(markerName=="S3",16,14))
  if (markerName %in% c("Min","AUC","S4")){ ticks <- seq(0,1,by=0.25) }
  if (markerName %in% c("S2","S3")){ ticks <- seq(0,0.75,by=0.25) }
  if (markerName=="S1"){ ticks <- seq(0,0.5,by=0.25) }
  axis(side=4, at=ticks)
  mtext("Density                          ", side=4, las=0, line=3, cex=1.1)
  
  data <- subset(dataFull, AGEYRS>11)
  dataI <- subset(dataFullI, AGEYRS>11)
  
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  colRGB2 <- c(col2rgb("olivedrab3"))
  colRGB2 <- rgb(colRGB2[1], colRGB2[2], colRGB2[3], alpha=255*0.4, maxColorValue=255)
  par(new=TRUE)
  hist(S.V, col=colRGB2, axes=FALSE, labels=FALSE, main="", xlab="", ylab="", border=0, xlim=xLim, ylim=yLim, freq=FALSE, breaks=ifelse(markerName=="S3",16,14))
  
  if (legOut){
    legend("bottomleft", lty=c("longdash","solid"), lwd=2.5, col=c("red3","black"), legend=c("9-11 years","12-16 years"), y.intersp=0.9, inset=c(0.07,-0.55), xpd=NA, bty="n", cex=0.7)
    legend("bottom", fill=c(colRGBle11,colRGBg11), border=c(colRGBle11,colRGBg11), legend=c("\n95% pointwise CIs  ",""), bty="n", y.intersp=0.3, inset=c(0,-0.533), xpd=NA, cex=0.7)
    legend("bottomright", fill=c(colRGB1, colRGB2), border=c(colRGB1, colRGB2), legend=c("Vaccine 9-11 years","Vaccine 12-16 years"), bty="n", y.intersp=0.9, inset=c(0.02,-0.55), xpd=NA, cex=0.7)
  } else {
    xLeg <- switch(markerName, Min=2.75, AUC=3, S1=3.3, S2=3, S3=3.5, S4=3)
    y1Leg <- switch(markerName, Min=1.65, AUC=1.65, S1=0.75, S2=1.1, S3=0.85, S4=1.62)
    y2Leg <- switch(markerName, Min=1.45, AUC=1.45, S1=0.67, S2=1, S3=0.75, S4=1.4)
    y3Leg <- switch(markerName, Min=1.1, AUC=1.1, S1=0.52, S2=0.8, S3=0.55, S4=1.05)
    legend(xLeg, y1Leg, lty=c("longdash","solid"), lwd=2.5, col=c("red3","black"), legend=c("9-11 years","12-16 years"), bty="n", cex=0.7, y.intersp=0.8)
    legend(xLeg, y2Leg, fill=c(colRGBle11,colRGBg11), border=c(colRGBle11,colRGBg11), legend=c("\n95% pointwise CIs",""), bty="n", cex=0.7, y.intersp=0.3)
    legend(xLeg, y3Leg, fill=c(colRGB1, colRGB2), border=c(colRGB1, colRGB2), legend=c("Vaccine 9-11 years","Vaccine 12-16 years"), bty="n", cex=0.7, y.intersp=0.8)
  }  
  
  epoint <- switch(markerName, Min="primary", AUC="primary", S1="matched serotype 1", S2="matched serotype 2", S3="matched serotype 3", S4="matched serotype 4")
  if (legOut){ cat("\\vspace*{-3mm}") }
  cat("\\caption[Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint comparing subgroups of 9--11 versus 12--16-year-olds, with the 95\\% pointwise confidence bands]{Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint comparing subgroups of 9--11 versus 12--16-year-olds, with the 95\\% pointwise confidence bands based on $10^3$ bootstrap iterations and the two phase-adjusted histograms of ",markerName3," measured in the vaccine group at month 13 in the CYD15 trial}\\label{Fig: VE(s1) ",markerName2,"_ageStrata}", sep="")
}


plotVEcurveOverlay3A <- function(dataFull, dataFullI, markerName, VEcurveFileAle5, VEcurveFileAg5le11, VEcurveFileAg11, bVEcurveFileA, VEcurveMin=0.8, legOut=FALSE, saveDir){
  markerName2 <- switch(markerName, AUC="AUCMB", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, AUC="the AUC-MB", S1="$\\log_{10}$ serotype 1 titers", S2="$\\log_{10}$ serotype 2 titers", S3="$\\log_{10}$ serotype 3 titers", S4="$\\log_{10}$ serotype 4 titers")
  markerVar <- paste0("IMPSTLOG.", markerName2)
  fstatusVar <- paste0(switch(markerName, AUC="o", S1="s1", S2="s2", S3="s3", S4="s4"),"fstatus_m13") 
  markerX <- switch(markerName, AUC="AUC-MB", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  
  # initiate the plot
  load(file.path(saveDir, VEcurveFileAle5))
  xLim <- c(min(dataFullI[,markerVar], na.rm=TRUE),4.5)
  par(mar=c(4+ifelse(legOut,3.5,0),5,3,4.2), las=1, cex.axis=0.9, cex.lab=1)
  plot(oList$markerVals[oList$markerVals>=VEcurveMin], oList$VE[oList$markerVals>=VEcurveMin], type="n",
       xlab=paste0("Vaccine-Induced ",markerX," at Month 13"), ylab="", ylim=c(-1.2,1), xaxt="n", yaxt="n", 
       xlim=xLim)
  axis(side=1, at=seq(1,4.5,by=0.5))
  axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
  mtext(expression(paste("Estimated VE(",s[1],")"), sep=""), side=2, las=0, line=3.6, cex=1.1)
  
  # interval estimates of VE(s_1)
  load(file.path(saveDir, bVEcurveFileA))  
  results <- resultsle5  
  if (markerName=="S3"){
    bVE <- NULL
    for (i in 1:length(results)){
      if (!is.null(names(results[[i]]))){
        bVE <- rbind(bVE, results[[i]]$bVE)
      }
    }
  } else {
    bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
  }
  s <- results[[1]]$markerVals
  bVE <- bVE[,s>=VEcurveMin]
  s <- s[s>=VEcurveMin]
  xCoord <- c(s, rev(s))
  yCoord <- c(apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE), rev(apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)))
  if (is.na(yCoord[36])){ yCoord[36] <- yCoord[35] + ((yCoord[38]-yCoord[35])/0.3)*0.1 }
  if (is.na(yCoord[37])){ yCoord[37] <- yCoord[35] + ((yCoord[38]-yCoord[35])/0.3)*0.2 }
  # for S4, in this age subgroup, the lower bound of the CI is not available for log10 titers >4
  if (markerName=="S4"){
    yCoord <- yCoord[xCoord<=4]
    xCoord <- xCoord[xCoord<=4]
  }
  colRGBle5 <- c(col2rgb("darkgoldenrod2"))
  colRGBle5 <- rgb(colRGBle5[1], colRGBle5[2], colRGBle5[3], alpha=255*0.55, maxColorValue=255)
  polygon(xCoord, yCoord, col=colRGBle5, border=NA)
  
  results <- resultsbet5_11
  if (markerName=="S3"){
    bVE <- NULL
    for (i in 1:length(results)){
      if (!is.null(names(results[[i]]))){
        bVE <- rbind(bVE, results[[i]]$bVE)
      }
    }
  } else {
    bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
  }
  s <- results[[1]]$markerVals
  bVE <- bVE[,s>=VEcurveMin]
  s <- s[s>=VEcurveMin]
  xCoord <- c(s, rev(s))
  yCoord <- c(apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE), rev(apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)))
  colRGBg5le11 <- c(col2rgb("red1"))
  colRGBg5le11 <- rgb(colRGBg5le11[1], colRGBg5le11[2], colRGBg5le11[3], alpha=255*0.45, maxColorValue=255)
  polygon(xCoord, yCoord, col=colRGBg5le11, border=NA)
  
  results <- resultsgt11
  if (markerName=="S3"){
    bVE <- NULL
    for (i in 1:length(results)){
      if (!is.null(names(results[[i]]))){
        bVE <- rbind(bVE, results[[i]]$bVE)
      }
    }
  } else {
    bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
  }
  s <- results[[1]]$markerVals
  bVE <- bVE[,s>=VEcurveMin]
  s <- s[s>=VEcurveMin]
  xCoord <- c(s, rev(s))
  yCoord <- c(apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE), rev(apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)))
  colRGBg11 <- c(col2rgb("gray40"))
  colRGBg11 <- rgb(colRGBg11[1], colRGBg11[2], colRGBg11[3], alpha=255*0.5, maxColorValue=255)
  polygon(xCoord, yCoord, col=colRGBg11, border=NA)
  
  # point estimates of VE(s_1)
  load(file.path(saveDir, VEcurveFileAle5))
  lines(oList$markerVals[oList$markerVals>=VEcurveMin], oList$VE[oList$markerVals>=VEcurveMin], col="darkgoldenrod2", lwd=2.5, lty="dotdash")
  load(file.path(saveDir, VEcurveFileAg5le11))
  lines(oList$markerVals[oList$markerVals>=VEcurveMin], oList$VE[oList$markerVals>=VEcurveMin], col="red3", lwd=2.5, lty="longdash")
  load(file.path(saveDir, VEcurveFileAg11))
  lines(oList$markerVals[oList$markerVals>=VEcurveMin], oList$VE[oList$markerVals>=VEcurveMin], col="black", lwd=2.5)
  
  # sampling weights in hist
  data <- subset(dataFull, AGEYRS<=5)
  dataI <- subset(dataFullI, AGEYRS<=5)  
  
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  ylim <- c(0, max(hist(S.V[S.V<=4.5], plot=FALSE, breaks=14)$density)*2.5)
  colRGB1 <- c(col2rgb("darkgoldenrod2"))
  colRGB1 <- rgb(colRGB1[1], colRGB1[2], colRGB1[3], alpha=255*0.55, maxColorValue=255)
  par(new=TRUE)
  hist(S.V[S.V<=4.5], col=colRGB1, axes=FALSE, labels=FALSE, 
       main="", xlab="", ylab="", border=0, freq=FALSE, xlim=xLim, ylim=ylim, breaks=14)
  axis(side=4, at=seq(0,1,by=0.25))
  mtext("Density                              ", side=4, las=0, line=3, cex=1.1)
  
  data <- subset(dataFull, AGEYRS>5 & AGEYRS<=11)
  dataI <- subset(dataFullI, AGEYRS>5 & AGEYRS<=11)
  
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  colRGB2 <- c(col2rgb("red1"))
  colRGB2 <- rgb(colRGB2[1], colRGB2[2], colRGB2[3], alpha=255*0.45, maxColorValue=255)
  par(new=TRUE)
  hist(S.V[S.V<=4.5], col=colRGB2, axes=FALSE, labels=FALSE, 
       main="", xlab="", ylab="", border=0, xlim=xLim, ylim=ylim, freq=FALSE, breaks=14)
  
  data <- subset(dataFull, AGEYRS>11)
  dataI <- subset(dataFullI, AGEYRS>11)  
  
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  colRGB3 <- c(col2rgb("gray40"))
  colRGB3 <- rgb(colRGB3[1], colRGB3[2], colRGB3[3], alpha=255*0.5, maxColorValue=255)
  par(new=TRUE)
  hist(S.V[S.V<=4.5], col=colRGB3, axes=FALSE, labels=FALSE, 
       main="", xlab="", ylab="", border=0, freq=FALSE, xlim=xLim, ylim=ylim, breaks=14)
  
  if (legOut){
    legend("bottomleft", lty=c("dotdash","longdash","solid"), lwd=2.5, col=c("darkgoldenrod2","red3","black"), legend=c("2-5 years","6-11 years","12-14 years"), y.intersp=0.9, inset=c(0.07,-0.55), xpd=NA, bty="n", cex=0.7)
    legend("bottom", fill=c(colRGBle5,colRGBg5le11,colRGBg11), border=c(colRGBle5,colRGBg5le11,colRGBg11), legend=c("\n95% pointwise CIs  ",""), bty="n", y.intersp=0.3, inset=c(0,-0.533), xpd=NA, cex=0.7)
    legend("bottomright", fill=c(colRGB1,colRGB2,colRGB3), border=c(colRGB1,colRGB2,colRGB3), legend=c("Vaccine 2-5 years","Vaccine 6-11 years","Vaccine 12-14 years"), bty="n", y.intersp=0.9, inset=c(0.02,-0.55), xpd=NA, cex=0.7)
  } else {
    xLeg1 <- 2
    xLeg2 <- 3
    xLeg3 <- 3
    yLeg1 <- 1.8
    yLeg2 <- 1.8
    yLeg3 <- 1
    legend(xLeg1, yLeg1, lty=c("dotdash","longdash","solid"), lwd=2.5, col=c("darkgoldenrod2","red3","black"), legend=c("2-5 years","6-11 years","12-14 years"), bty="n", cex=0.7, y.intersp=0.8)
    legend(xLeg2, yLeg2, fill=c(colRGBle5,colRGBg5le11,colRGBg11), border=c(colRGBle5,colRGBg5le11,colRGBg11), legend=c("","95% pointwise CIs",""), bty="n", cex=0.7, y.intersp=0.6)
    legend(xLeg3, yLeg3, fill=c(colRGB1,colRGB2,colRGB3), border=c(colRGB1,colRGB2,colRGB3), legend=c("Vaccine 2-5 years","Vaccine 6-11 years","Vaccine 12-14 years"), bty="n", cex=0.7, y.intersp=0.8)
  }  
  
  epoint <- switch(markerName, AUC="overall", S1="serotype 1", S2="serotype 2", S3="serotype 3", S4="serotype 4")
  if (legOut){ cat("\\vspace*{-3mm}") }
  cat("\\caption{Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint comparing subgroups of 2--5, 6--11, and 12--14-year-olds, with the 95\\% pointwise confidence bands based on $10^3$ bootstrap iterations and the two phase-adjusted histograms of ",markerName3," measured in the vaccine group at month 13 in the CYD14 trial}", sep="")
}

# 'plaRiskCurvePointEst' and 'txRiskCurvePointEst' are vectors
# 'plaRiskCurveBootEst' and 'txRiskCurveBootEst' are length(markerVals)-by-1000 matrices
plotLogRRcurve <- function(markerVals, plaRiskCurvePointEst, txRiskCurvePointEst, plaRiskCurveBootEst=NULL, txRiskCurveBootEst=NULL, title=NULL, hingePoint=NULL, plotLegend=TRUE){
  MCEPcurvePointEst <- log(txRiskCurvePointEst/plaRiskCurvePointEst)
  
  if (!is.null(plaRiskCurveBootEst) && !is.null(txRiskCurveBootEst)){
    # transformed MCEP curve (identity)
    tMCEP <- MCEPcurvePointEst
    # transformed bootstrapped MCEP curves (identity)
    # assuming the matrices have the same dimensions
    tbMCEP <- log(txRiskCurveBootEst/plaRiskCurveBootEst)
    
    # bootstrap SE of tMCEP estimates
    bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)
    
    # pointwise confidence bounds for MCEP(s1)
    ptLB.MCEP <- tMCEP - qnorm(0.975) * bSE
    ptUB.MCEP <- tMCEP + qnorm(0.975) * bSE
    
    supAbsZ <- NULL
    for (j in 1:NCOL(tbMCEP)){
      Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
      supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
    }
    qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
    
    smLB.MCEP <- tMCEP - qSupAbsZ * bSE
    smUB.MCEP <- tMCEP + qSupAbsZ * bSE
  } else {
    ptLB.MCEP <- ptUB.MCEP <- smLB.MCEP <- smUB.MCEP <- NULL
  }
  
  cexTitle <- 1.7
  cexLab <- 1.4
  cexAxis <- 1.3
  cexLegend <- 1.2
  
  par(mar=c(5,5,1.5,5), cex.lab=cexLab, cex.axis=cexAxis, las=1)
  
  plot(markerVals, MCEPcurvePointEst, type="l", xlab="Month 13 Average Titer of Vaccinees", ylab="", xlim=range(markerVals), 
       ylim=range(c(MCEPcurvePointEst, ptLB.MCEP, ptUB.MCEP, smLB.MCEP, smUB.MCEP), na.rm=TRUE), lwd=3.5, xaxt="n", yaxt="n")
  axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
  axis(side=2, at=seq(-5,0,by=0.5), cex.axis=cexAxis)
  mtext("Log Relative Risk", side=2, las=0, line=3.4, cex=cexLab)
  axis(side=4, at=log(1-c(seq(0,0.9,by=0.1),0.95,0.99)), labels=c(seq(0,0.9,by=0.1),0.95,0.99)*100, cex.axis=cexAxis)
  mtext("Vaccine Efficacy (%)", side=4, las=0, line=3.2, cex=cexLab)
  if (!is.null(title)){ mtext(title, side=3, cex=cexTitle, line=0) }
  
  if (!is.null(plaRiskCurveBootEst) && !is.null(txRiskCurveBootEst)){
    lines(markerVals, ptLB.MCEP, lty="dashed", lwd=3)
    lines(markerVals, ptUB.MCEP, lty="dashed", lwd=3)
    lines(markerVals, smLB.MCEP, lty="dotdash", lwd=3)
    lines(markerVals, smUB.MCEP, lty="dotdash", lwd=3)
  }
  
  if (plotLegend){ legend("bottomleft", lty=c("dashed","dotdash"), lwd=3, legend=c("Pointwise 95% CI","Simultaneous 95% CI"), cex=cexLegend, bty="n") }
  if (!is.null(hingePoint)){ legend("topright", paste0("Hinge Point = ", hingePoint,"  "), cex=cexLegend, bty="n") }
}

# 'plaRiskCurvePointEst' and 'txRiskCurvePointEst' are vectors
# 'plaRiskCurveBootEst' and 'txRiskCurveBootEst' are length(markerVals)-by-1000 matrices
plotRiskDiffCurve <- function(markerVals, plaRiskCurvePointEst, txRiskCurvePointEst, plaRiskCurveBootEst=NULL, txRiskCurveBootEst=NULL, title=NULL, hingePoint=NULL, plotLegend=TRUE){
  MCEPcurvePointEst <- plaRiskCurvePointEst - txRiskCurvePointEst
  
  if (!is.null(plaRiskCurveBootEst) && !is.null(txRiskCurveBootEst)){
    # transformed MCEP curve (identity)
    tMCEP <- MCEPcurvePointEst
    # transformed bootstrapped MCEP curves (identity)
    # assuming the matrices have the same dimensions
    tbMCEP <- plaRiskCurveBootEst - txRiskCurveBootEst
    
    # bootstrap SE of tMCEP estimates
    bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)
    
    # pointwise confidence bounds for MCEP(s1)
    ptLB.MCEP <- tMCEP - qnorm(0.975) * bSE
    ptUB.MCEP <- tMCEP + qnorm(0.975) * bSE
    
    supAbsZ <- NULL
    for (j in 1:NCOL(tbMCEP)){
      Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
      supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
    }
    qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
    
    smLB.MCEP <- tMCEP - qSupAbsZ * bSE
    smUB.MCEP <- tMCEP + qSupAbsZ * bSE
  } else {
    ptLB.MCEP <- ptUB.MCEP <- smLB.MCEP <- smUB.MCEP <- NULL
  }
  
  cexTitle <- 1.7
  cexLab <- 1.4
  cexAxis <- 1.3
  cexLegend <- 1.2
  
  par(mar=c(5,5,1.5,5), cex.lab=cexLab, cex.axis=cexAxis, las=1)
  
  plot(markerVals, MCEPcurvePointEst, type="l", xlab="Month 13 Average Titer of Vaccinees", ylab="", xlim=range(markerVals), 
       ylim=range(c(MCEPcurvePointEst, ptLB.MCEP, ptUB.MCEP, smLB.MCEP, smUB.MCEP), na.rm=TRUE), lwd=3.5, xaxt="n", yaxt="n")
  axis(side=1, at=c(log10(5),1:5), labels=expression("<10   ",10,100,10^3,10^4,10^5))
  axis(side=2, at=seq(0.01,0.04,by=0.002), labels=FALSE, cex.axis=cexAxis)
  axis(side=2, at=seq(0.01,0.04,by=0.002), cex.axis=cexAxis, line=-0.3, tick=FALSE)
  mtext("Risk Difference (Placebo - Vaccine)", side=2, las=0, line=3.8, cex=cexLab)
  if (!is.null(title)){ mtext(title, side=3, cex=cexTitle, line=0) }
  
  if (!is.null(plaRiskCurveBootEst) && !is.null(txRiskCurveBootEst)){
    lines(markerVals, ptLB.MCEP, lty="dashed", lwd=3)
    lines(markerVals, ptUB.MCEP, lty="dashed", lwd=3)
    lines(markerVals, smLB.MCEP, lty="dotdash", lwd=3)
    lines(markerVals, smUB.MCEP, lty="dotdash", lwd=3)
  }
  
  if (plotLegend){ legend("bottomleft", lty=c("dashed","dotdash"), lwd=3, legend=c("Pointwise 95% CI","Simultaneous 95% CI"), cex=cexLegend, bty="n") }
  if (!is.null(hingePoint)){ legend("topright", paste0("Hinge Point = ", hingePoint,"  "), cex=cexLegend, bty="n") }
}

plotVEcurveBS <- function(data, dataI, markerName, VEcurveFile, bVEcurveFile, saveDir){
  markerName2 <- switch(markerName, AUC="AUCMB", S1="Sero1", S2="Sero2", S3="Sero3", S4="Sero4")
  markerName3 <- switch(markerName, AUC="AUC-MB", S1="$\\log_{10}$ serotype 1 titers", S2="$\\log_{10}$ serotype 2 titers", S3="$\\log_{10}$ serotype 3 titers", S4="$\\log_{10}$ serotype 4 titers")
  markerVar <- paste0("IMPSTLOG.", markerName2)
  bmarkerVar <- paste0("b", markerName)
  fstatusVar <- paste0(switch(markerName, AUC="o", S1="s1", S2="s2", S3="s3", S4="s4"),"fstatus_m13") 
  markerX <- switch(markerName, AUC="AUC-MB", S1="Log10 Serotype 1 Titer", S2="Log10 Serotype 2 Titer", S3="Log10 Serotype 3 Titer", S4="Log10 Serotype 4 Titer")
  
  load(file.path(saveDir, VEcurveFile))
  xLim <- c(min(dataI[,markerVar]-dataI[,bmarkerVar], na.rm=TRUE), max(dataI[,markerVar]-dataI[,bmarkerVar], na.rm=TRUE))
  par(mar=c(4,5,3,4.2), las=1, cex.axis=0.9, cex.lab=1)
  plot(oList$markerVals, oList$VE, type="l", 
       xlab=paste0("Baseline-Subtracted Vaccine-Induced ",markerX," at Month 13"), ylab="", xlim=xLim, ylim=c(-1,1), yaxt="n", 
       col="red3", lwd=2.5)
  #axis(side=1, at=seq(1,4.5,by=0.5))
  axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
  mtext(expression(paste("Estimated VE(",s[1],")"), sep=""), side=2, las=0, line=3.6, cex=1.1)
  
  if (!is.null(bVEcurveFile)){
    load(file.path(saveDir, bVEcurveFile))
    bVE <- NULL
    for (i in 1:length(results)){
      if (!is.null(names(results[[i]]))){
        bVE <- rbind(bVE, results[[i]]$bVE)
      }
    }
    #bVE <- do.call(rbind, lapply(results, "[[", "bVE"))     
    
    s <- results[[1]]$markerVals
    lines(s, apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE), lty="dashed", lwd=2, col="black")
    lines(s, apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE), lty="dashed", lwd=2, col="black")
  }
  
  colHist <- c(col2rgb("olivedrab3"))
  colHist <- rgb(colHist[1], colHist[2], colHist[3], alpha=255*0.4, maxColorValue=255)
  xLeg <- 1.2
  yLeg1 <- 1
  yLeg2 <- 0.9
  legend(xLeg, yLeg1, lty="dashed", lwd=2, col="black", legend="95% pointwise CI", bty="n", cex=0.7)
  legend(xLeg, yLeg2, fill=colHist, border=colHist, legend="Vaccine", bty="n", cex=0.7)  
  
  # sampling weights in hist
  data$fstatus <- data[,fstatusVar]
  data$marker <- data[,markerVar]
  data$bmarker <- data[,bmarkerVar]
  data$marker <- data$marker - data$bmarker
  wt.PControls <- NROW(subset(data, VACC==0 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==0))
  wt.VControls <- NROW(subset(data, VACC==1 & fstatus==0))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==0))
  wt.PCases <- NROW(subset(data, VACC==0 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==0 & fstatus==1))
  wt.VCases <- NROW(subset(data, VACC==1 & fstatus==1))/NROW(subset(data, !is.na(marker) & VACC==1 & fstatus==1))
  
  dataI$fstatus <- dataI[,fstatusVar]
  dataI$marker <- dataI[,markerVar]
  dataI$bmarker <- dataI[,bmarkerVar]
  dataI$marker <- dataI$marker - dataI$bmarker
  S.PControls <- subset(dataI, VACC==0 & fstatus==0)$marker
  S.VControls <- subset(dataI, VACC==1 & fstatus==0)$marker
  S.PCases <- subset(dataI, VACC==0 & fstatus==1)$marker
  S.VCases <- subset(dataI, VACC==1 & fstatus==1)$marker
  
  S.P <- c(rep(S.PControls, round(wt.PControls)), rep(S.PCases, round(wt.PCases)))
  S.V <- c(rep(S.VControls, round(wt.VControls)), rep(S.VCases, round(wt.VCases)))
  
  yLim <- c(0, max(hist(S.V, plot=FALSE, breaks=14)$density)*2.5)
  par(new=TRUE)
  hist(S.V, col=colHist, axes=FALSE, labels=FALSE, 
       main="", xlab="", ylab="", border=0, freq=FALSE, xlim=xLim, ylim=yLim, breaks=16)
  axis(side=4, at=seq(0,0.75,by=0.25))
  #par(new=TRUE)
  #hist(subset(dataI, VACC==0 & IMPSTLOG.AUCMB<=4)$IMPSTLOG.AUCMB, col=rgb(0,0,1,alpha=0.22), axes=FALSE, labels=FALSE, main="",xlab="", ylab="", border=0, freq=TRUE, ylim=ylim, breaks=10)
  mtext("Density                          ", side=4, las=0, line=3, cex=1.1)
  
  epoint <- switch(markerName, AUC="primary", S1="matched serotype 1", S2="matched serotype 2", S3="matched serotype 3", S4="matched serotype 4")  
  cat("\\caption[Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint with the 95\\% pointwise confidence band by baseline-subtracted AUC-MB]{Estimated $VE(s_1)$ against the ",epoint, " dengue endpoint with the 95\\% pointwise confidence band based on $10^3$ bootstrap iterations and the two phase-adjusted histogram of the baseline-subtracted ",markerName3," measured in the vaccine group at month 13 in the CYD15 trial}\\label{Fig: VE(s1) ",markerName2,"_baselineSubtracted}", sep="")
}

# 'applyBoot' calculates the bootstrap confidence bands as pointwise quantiles from the large number of VE curves
# the purpose of this function is to speed up plotting and report generation and to avoid running out of memory
# this function creates a new .RData file that can be loaded in a plotting function
applyBoot <- function(loadFile, saveFile, saveDir){
  load(file.path(saveDir, loadFile))
  bVE <- NULL
  for (i in 1:length(results)){
    if (!is.null(names(results[[i]]))){
      bVE <- rbind(bVE, results[[i]]$bVE)
    }
  }    
  #bVE <- do.call(rbind, lapply(results, "[[", "bVE"))
  
  s <- results[[1]]$markerVals
  UB <- apply(bVE, 2, quantile, probs=0.975, na.rm=TRUE)
  LB <- apply(bVE, 2, quantile, probs=0.025, na.rm=TRUE)
  bList <- list(s=s, UB=UB, LB=LB)
  save(bList, file=file.path(saveDir,saveFile))
}
