library(MASS)
library(mvtnorm)
library(osDesign)
library(np)
library(splines)

# determine beta_0 to achieve a prespecified marginal probability of infection in the placebo group
# 'scenario' is one of "A1" and "A2"
getBeta0 <- function(beta2, meanS0, varS0, prob, scenario="A1", threshold=NULL){
  if (scenario=="A1"){
    f <- function(beta0, beta2, meanS0, varS0, prob){ integrate(function(s0, beta0, beta2, meanS0, varS0){ pnorm(beta0+beta2*s0)*dnorm(s0, mean=meanS0, sd=sqrt(varS0)) }, -Inf, Inf, beta0=beta0, beta2=beta2, meanS0=meanS0, varS0=varS0)$value - prob }  
  }
  
  if (scenario=="A2"){
    f <- function(beta0, beta2, meanS0, varS0, prob){ pnorm(beta0 + beta2*threshold)*pnorm(threshold, mean=meanS0, sd=sqrt(varS0)) + integrate(function(s0, beta0, beta2, meanS0, varS0){ pnorm(beta0+beta2*s0)*dnorm(s0, mean=meanS0, sd=sqrt(varS0)) }, threshold, Inf, beta0=beta0, beta2=beta2, meanS0=meanS0, varS0=varS0)$value - prob }  
  }
  
  return(uniroot(f, c(-10,10), beta2=beta2, meanS0=meanS0, varS0=varS0, prob=prob)$root)
}

# determine the starting value of beta_1 such that p_1(0,0)/p_0(0,0) = rr
startBeta1 <- function(beta0, rr){
  f <- function(beta1, beta0, rr){ pnorm(beta0+beta1)/pnorm(beta0) - rr }
  return(uniroot(f, c(-10,10), beta0=beta0, rr=rr)$root)
}

# determine beta_3 such that p_1(s_{1,0.9})/p_1(s_{1,0.1}) = rr, where s_{1,p} is the quantile of the distribution of S(1) at probability p
getBeta3 <- function(beta0, beta1, beta2, meanS1, varS1, rr){
  # 'q' is assumed to be a numeric vector of length 2
  f <- function(beta3, beta0, beta1, beta2, q, rr){ pnorm(beta0+beta1+beta3*q[2])/pnorm(beta0+beta1+beta3*q[1]) - rr }
  q <- qnorm(c(0.1,0.9), mean=meanS1, sd=sqrt(varS1))
  return(uniroot(f, c(-10,10), beta0=beta0, beta1=beta1, beta2=beta2, q=q, rr=rr)$root)
}

# update beta_1 to achieve a prespecified marginal probability of infection in the vaccine group
# 'scenario' is one of "A1" and "A2"
getBeta1 <- function(beta0, beta2, beta3, meanS1, varS1, prob, scenario="A1", threshold=NULL){
  if (scenario=="A1"){
    f <- function(beta1, beta0, beta2, beta3, meanS1, varS1, prob){ 
      integrate(function(s1, beta0, beta1, beta2, beta3, meanS1, varS1){ pnorm(beta0+beta1+beta3*s1)*dnorm(s1, mean=meanS1, sd=sqrt(varS1)) }, -Inf, Inf, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, meanS1=meanS1, varS1=varS1)$value - prob
    }  
  }
  
  if (scenario=="A2"){
    f <- function(beta1, beta0, beta2, beta3, meanS1, varS1, prob){ 
      pnorm(beta0+beta1+beta3*threshold)*pnorm(threshold, mean=meanS1, sd=sqrt(varS1)) + integrate(function(s1, beta0, beta1, beta2, beta3, meanS1, varS1){pnorm(beta0+beta1+beta3*s1)*dnorm(s1, mean=meanS1, sd=sqrt(varS1)) }, threshold, Inf, beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, meanS1=meanS1, varS1=varS1)$value - prob
    }  
  }
  
  return(uniroot(f, c(-10,10), beta0=beta0, beta2=beta2, beta3=beta3, meanS1=meanS1, varS1=varS1, prob=prob)$root)
}

# 'dmvnorm.s1vector' calculates f(s0,s1) where f is bivariate normal density
# 's0' is a scalar
# 's1' is a numeric vector
dmvnorm.s1vector <- function(s0, s1, meanS0, varS0, meanS1, varS1, covS0S1){
  dmnormVector <- sapply(s1, function(s1val, s0, meanS0, varS0, meanS1, varS1, covS0S1){ 
    dmvnorm(c(s0,s1val), mean=c(meanS0, meanS1), sigma=matrix(c(varS0, covS0S1, covS0S1, varS1),2,2)) 
  }, s0=s0, meanS0=meanS0, varS0=varS0, meanS1=meanS1, varS1=varS1, covS0S1=covS0S1)
  return(dmnormVector)
}

# 'dmvnorm.s0vector' calculates f(s0,s1) where f is bivariate normal density
# 's0' is a numeric vector
# 's1' is a scalar
dmvnorm.s0vector <- function(s0, s1, meanS0, varS0, meanS1, varS1, covS0S1){
  dmnormVector <- sapply(s0, function(s0val, s1, meanS0, varS0, meanS1, varS1, covS0S1){ 
    dmvnorm(c(s0val,s1), mean=c(meanS0, meanS1), sigma=matrix(c(varS0, covS0S1, covS0S1, varS1),2,2)) 
  }, s1=s1, meanS0=meanS0, varS0=varS0, meanS1=meanS1, varS1=varS1, covS0S1=covS0S1)
  return(dmnormVector)
}

# 's0' is a scalar
pS0 <- function(s0, meanS0, varS0, meanS1, varS1, covS0S1, threshold){ 
  return(integrate(function(s1, s0, meanS0, varS0, meanS1, varS1, covS0S1){ dmvnorm.s1vector(s0, s1, meanS0, varS0, meanS1, varS1, covS0S1) }, 
                   -Inf, threshold, s0=s0, meanS0=meanS0, varS0=varS0, meanS1=meanS1, varS1=varS1, covS0S1=covS0S1)$value / pnorm(threshold, mean=meanS1, sd=sqrt(varS1)))
}

# 's0' is a numeric vector
wRiskS0 <- function(s0, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1, threshold){
  wRisk <- sapply(s0, function(s0val, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1){ 
    pnorm(beta0+beta2*s0val)*pS0(s0val, meanS0, varS0, meanS1, varS1, covS0S1, threshold) 
    }, beta0=beta0, beta2=beta2, meanS0=meanS0, varS0=varS0, meanS1=meanS1, varS1=varS1, covS0S1=covS0S1)
  return(wRisk)
}

# 's1' is a scalar
p0S <- function(s1, meanS0, varS0, meanS1, varS1, covS0S1, threshold){ 
  return(integrate(function(s0, s1, meanS0, varS0, meanS1, varS1, covS0S1){ dmvnorm.s0vector(s0, s1, meanS0, varS0, meanS1, varS1, covS0S1) }, 
                   -Inf, threshold, s1=s1, meanS0=meanS0, varS0=varS0, meanS1=meanS1, varS1=varS1, covS0S1=covS0S1)$value / dnorm(s1, mean=meanS1, sd=sqrt(varS1)))
}

# 's0' is a numeric vector
# 's1' is a scalar
wRisk0S <- function(s0, s1, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1){
  return(pnorm(beta0+beta2*s0) * dmvnorm.s0vector(s0, s1, meanS0, varS0, meanS1, varS1, covS0S1) / dnorm(s1, mean=meanS1, sd=sqrt(varS1)))
}

# calculate P(Y(0)=1 | S(1)=s1)
# 's1' is assumed to be a scalar
# 'scenario' is one of "A1" and "A2"
risk0 <- function(s1, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1, scenario="A1", threshold=NULL){
  if (scenario=="A1"){
    p <- integrate(function(s0, s1, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1){ 
      pnorm(beta0+beta2*s0)*dnorm(s0, mean=meanS0 + covS0S1*(s1-meanS1)/varS1, sd=sqrt(varS0 - (covS0S1^2)/varS1)) 
      }, -Inf, Inf, s1=s1, beta0=beta0, beta2=beta2, meanS0=meanS0, varS0=varS0, meanS1=meanS1, varS1=varS1, covS0S1=covS0S1
      )$value
  }
  
  if (scenario=="A2"){
    
    if (s1==threshold){
      
      p00 <- pmvnorm(lower=rep(-Inf,2), upper=rep(threshold,2), mean=c(meanS0,meanS1), sigma=matrix(c(varS0, covS0S1, covS0S1, varS1),2,2)) / pnorm(threshold, mean=3, sd=sqrt(varS1))
      
      p <- pnorm(beta0+beta2*threshold)*p00 + integrate(function(s0, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1, threshold){ 
        wRiskS0(s0, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1, threshold) }, threshold, Inf, 
        beta0=beta0, beta2=beta2, meanS0=meanS0, varS0=varS0, meanS1=meanS1, varS1=varS1, covS0S1=covS0S1, threshold=threshold)$value
      
    } else { # i.e., s1 > threshold
      
      p <- pnorm(beta0+beta2*threshold)*p0S(s1, meanS0, varS0, meanS1, varS1, covS0S1, threshold) + 
        integrate(function(s0, s1, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1){ 
          wRisk0S(s0, s1, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1)
          }, threshold, Inf, s1=s1, beta0=beta0, beta2=beta2, meanS0=meanS0, varS0=varS0, meanS1=meanS1, varS1=varS1, covS0S1=covS0S1)$value
    }
  }
  
  return(p)
}

# the true VE curve assuming the probit model for the probability of infection
# 's1' is a vector
# 'scenario' is one of "A1" and "A2"
trueVE <- function(s1, beta0, beta1, beta2, beta3, meanS0, varS0, meanS1, varS1, covS0S1, scenario="A1", threshold=NULL){
  1 - pnorm(beta0+beta1+beta3*s1)/sapply(s1, function(s1Value){ risk0(s1Value, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1, scenario=scenario, threshold=threshold) }) 
}

# 'n' is the total sample size
# 'beta' is a vector of probit risk model coefficients
# 'pi' is a sampling probability for sampling a subcohort of controls with biomarker measurements
# 'truncateMarker' is TRUE indicates truncated marker distributions
getData <- function(n, beta, pi, truncateMarker, seed){
  set.seed(seed)
  
  Z <- rep(0:1, each=n/2)
  S <- mvrnorm(n, mu=c(2,2,3), Sigma=matrix(c(1,0.9,0.7,0.9,1,0.7,0.7,0.7,1), nrow=3))
  if (truncateMarker){
    S <- ifelse(S<1.5, 1.5, S)  
  }
  p <- pnorm(drop(cbind(1,Z,(1-Z)*S[,2],Z*S[,3]) %*% beta))
  Y <- sapply(p, function(risk){ rbinom(1,1,risk) })
  
  # delete S(1) in placebo recipients
  S[Z==0,3] <- NA
  
  # delete S(0) in vaccine recipients
  S[Z==1,2] <- NA
  
  # generate the indicator of being sampled into the Phase 2 subcohort
  phase2 <- rbinom(n,1,pi)
  S[Y==0 & phase2==0,] <- c(NA,NA,NA)
  
  # delete Sb for cases not included in the Phase 2 subcohort
  S[Y==1 & phase2==0,1] <- NA
  
  data <- data.frame(Z,S,Y)
  colnames(data) <- c("Z","Sb","S0","S1","Y")
  
  return(data)
}

# 'tpsPredict' returns predicted values from a model fitted by tps
# columns of newMatrix in the same order as the coefficient vector from tps
tpsPredict <- function(fit, newMatrix){
  linPred <- newMatrix %*% fit$coef
  return(drop(1/(1+exp(-linPred))))
}

# 'hNum' returns function values at s0 of the integrand in the numerator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1 is a scalar
hNum <- function(s0, s1, tpsFit, npcdensFit1, npcdensFit2){
  phat.s0 <- tpsPredict(tpsFit, cbind(1, s0))
  fhat.s0 <- predict(npcdensFit1, newdata=data.frame(Sb=s0, S1=s1))
  ghat.s0 <- predict(npcdensFit2, newdata=data.frame(S0=s0))
  return(phat.s0*fhat.s0*ghat.s0)
}

# 'phNum' returns function values at s0 of the integrand in the numerator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1 is a scalar
# 'phNum' uses parametric (Gaussian) density estimation, whereas 'hNum' uses nonparametric kernel density estimation
phNum <- function(s0, s1, tpsFit, lmFit1, lmFit2){
  phat.s0 <- tpsPredict(tpsFit, cbind(1, s0))
  fMean.s0 <- predict(lmFit1, newdata=data.frame(Sb=s0))
  fSD.s0 <- summary(lmFit1)$sigma
  fhat.s0 <- dnorm(s1, mean=fMean.s0, sd=fSD.s0)
  gMean.s0 <- predict(lmFit2, newdata=data.frame(1))
  gSD.s0 <- summary(lmFit2)$sigma
  ghat.s0 <- dnorm(s0, mean=gMean.s0, sd=gSD.s0)
  return(phat.s0*fhat.s0*ghat.s0)
}

# 'hDen' returns function values at s0 of the integrand in the denominator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1, x.male, x.age, x.country are scalars
hDen <- function(s0, s1, npcdensFit1, npcdensFit2){
  fhat.s0 <- predict(npcdensFit1, newdata=data.frame(Sb=s0, S1=s1))
  ghat.s0 <- predict(npcdensFit2, newdata=data.frame(S0=s0))
  return(fhat.s0*ghat.s0)
}

# 'phDen' returns function values at s0 of the integrand in the denominator of risk_{(0)}(s_1)
# s0 is a numeric vector, whereas s1 is a scalar
# 'phDen' uses parametric (Gaussian) density estimation, whereas 'hDen' uses nonparametric kernel density estimation
phDen <- function(s0, s1, lmFit1, lmFit2){
  fMean.s0 <- predict(lmFit1, newdata=data.frame(Sb=s0))
  fSD.s0 <- summary(lmFit1)$sigma
  fhat.s0 <- dnorm(s1, mean=fMean.s0, sd=fSD.s0)
  gMean.s0 <- predict(lmFit2, newdata=data.frame(1))
  gSD.s0 <- summary(lmFit2)$sigma
  ghat.s0 <- dnorm(s0, mean=gMean.s0, sd=gSD.s0)
  return(fhat.s0*ghat.s0)
}

# 'riskP' returns the value of risk_{(0)}(s1)
# s1 is a scalar
riskP <- function(s1, data, tpsFit, npcdensFit1, npcdensFit2){
  UL <- max(data$S0, na.rm=TRUE) + 0.2 # if integration over (0,Inf) fails, use (0,UL)
  
  num <- try(integrate(hNum, 0, Inf, s1=s1, tpsFit=tpsFit, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000)$value, silent=TRUE)
  if (inherits(num, 'try-error')){
    num <- try(integrate(hNum, 0, UL, s1=s1, tpsFit=tpsFit, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000)$value, silent=TRUE)
  }
  
  den <- try(integrate(hDen, 0, UL, s1=s1, npcdensFit1=npcdensFit1, npcdensFit2=npcdensFit2, subdivisions=2000, rel.tol=30*.Machine$double.eps^0.25)$value, silent=TRUE)
  
  out <- NULL
  if ((!inherits(num, 'try-error')) & (!inherits(den, 'try-error'))){ out <- num/den }
  
  return(out)
}

# 'riskP' returns the value of risk_{(0)}(s1)
# s1 is a scalar
# 'pRiskP' uses parametric (Gaussian) density estimation, whereas 'riskP' uses nonparametric kernel density estimation
pRiskP <- function(s1, data, tpsFit, lmFit1, lmFit2){
  UL <- max(data$S0, na.rm=TRUE) + 0.2 # if integration over (0,Inf) fails, use (0,UL)
  
  num <- try(integrate(phNum, 0, Inf, s1=s1, tpsFit=tpsFit, lmFit1=lmFit1, lmFit2=lmFit2, subdivisions=2000)$value, silent=TRUE)
  if (inherits(num, 'try-error')){
    num <- try(integrate(phNum, 0, UL, s1=s1, tpsFit=tpsFit, lmFit1=lmFit1, lmFit2=lmFit2, subdivisions=2000)$value, silent=TRUE)
  }
  
  den <- try(integrate(phDen, 0, UL, s1=s1, lmFit1=lmFit1, lmFit2=lmFit2, subdivisions=2000, rel.tol=30*.Machine$double.eps^0.25)$value, silent=TRUE)
  
  out <- NULL
  if ((!inherits(num, 'try-error')) & (!inherits(den, 'try-error'))){ out <- num/den }
  
  return(out)
}

# 'riskV' returns the value of risk_{(1)}(s1)
# s1 is a scalar
riskV <- function(s1, data, dataI){
  nVControls <- NROW(subset(data, Z==1 & Y==0))    
  nVCases <- NROW(subset(data, Z==1 & Y==1))
  group <- rep(1, NROW(subset(dataI, Z==1 & !is.na(Y))))
  fit <- tps(Y ~ S1, data=subset(dataI, Z==1 & !is.na(Y)), nn0=nVControls, nn1=nVCases, group=group, method="PL", cohort=TRUE)
  return(tpsPredict(fit, cbind(1, s1)))
}

# 'estVE' returns the value of VE(s1)
# s1 is a scalar
estVE <- function(s1, data, dataI, tpsFit, npcdensFit1, npcdensFit2){ 
  riskPvalue <- riskP(s1, data, tpsFit, npcdensFit1, npcdensFit2)
  
  VE <- NA
  if (!is.null(riskPvalue)){ VE <- 1 - riskV(s1, data, dataI)/riskPvalue }
  
  return(VE)
}

# 'pEstVE' returns the value of VE(s1)
# s1 is a scalar
pEstVE <- function(s1, data, dataI, tpsFit, lmFit1, lmFit2){
  riskPvalue <- pRiskP(s1, data, tpsFit, lmFit1, lmFit2)
  
  VE <- NA
  if (!is.null(riskPvalue)){ VE <- 1 - riskV(s1, data, dataI)/riskPvalue }
  
  return(VE)
}

# 'VEcurve' returns the estimated VE(s1) curve evaluated on the grid of the 's1grid' values
# 'data' is a data frame with variables Z, Sb, S0, S1, and Y
VEcurve <- function(data, s1grid){
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(S0) | !is.na(S1))
  
  # extract subsets of controls ('dataControls') and cases ('dataCases') to be used for resampling
  # in addition, within each treatment group in the immunogenicity set, delete cases to recover
  # the case:control ratio in the Phase 1 population
  dataControls <- subset(data, Y==0)
  nPControlsI <- NROW(dataPControlsI <- subset(dataI, Z==0 & Y==0))
  nVControlsI <- NROW(dataVControlsI <- subset(dataI, Z==1 & Y==0))  
  
  dataCases <- subset(data, Y==1)  
  nPCasesI <- NROW(dataPCasesI <- subset(dataI, Z==0 & Y==1))
  nVCasesI <- NROW(dataVCasesI <- subset(dataI, Z==1 & Y==1))
  
  nPControls <- NROW(subset(dataControls, Z==0))
  nVControls <- NROW(subset(dataControls, Z==1))
  nPCases <- NROW(subset(dataCases, Z==0))
  nVCases <- NROW(subset(dataCases, Z==1))
  
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
  fbw <- npcdensbw(S1 ~ Sb, data=dataVIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  gbw <- npudensbw(~ S0, data=dataPIcorrectRatio, ckertype="epanechnikov")
  
  group <- rep(1, NROW(subset(dataI, Z==0)))
  
  # weighted logistic regression model using the placebo group in the immunogenicity set
  fit1 <- tps(Y ~ S0, data=subset(dataI, Z==0), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
  
  # kernel density estimator for f(s1|Sb=sb) using the vaccine group in the immunogenicity set
  fhat <- npcdens(fbw)
  
  # kernel density estimator for g(s0) using the placebo group in the immunogenicity set
  ghat <- npudens(gbw)
  
  VEcurvePointEst <- sapply(s1grid, function(s1val){ estVE(s1val, data, dataI, fit1, fhat, ghat) })
  
  return(VEcurvePointEst)
}

# 'pVEcurve' is a "parametric" version of 'VEcurve' with parametric (Gaussian) estimates of conditional densities
# 'pVEcurve' returns the estimated VE(s1) curve evaluated on the grid of the 's1grid' values
# 'data' is a data frame with variables Z, Sb, S0, S1, and Y
pVEcurve <- function(data, s1grid){
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(S0) | !is.na(S1))
  
  # calculate the sampling weights
  nPControls <- NROW(subset(data, Z==0 & Y==0))
  wtPControls <- NROW(subset(data, Z==0 & Y==0))/NROW(subset(dataI, Z==0 & Y==0))
  wtVControls <- NROW(subset(data, Z==1 & Y==0))/NROW(subset(dataI, Z==1 & Y==0))
  
  nPCases <- NROW(subset(data, Z==0 & Y==1))
  wtPCases <- NROW(subset(data, Z==0 & Y==1))/NROW(subset(dataI, Z==0 & Y==1))
  wtVCases <- NROW(subset(data, Z==1 & Y==1))/NROW(subset(dataI, Z==1 & Y==1))
  group <- rep(1, NROW(subset(dataI, Z==0)))
  
  # weighted logistic regression model using the placebo group in the immunogenicity set
  fit1 <- tps(Y ~ S0, data=subset(dataI, Z==0), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
  
  # normal density estimator for f(s1|Sb=sb) using the vaccine group in the immunogenicity set with baseline markers
  # sampling weights are incorporated
  dataB <- subset(dataI, Z==1 & !is.na(Sb))
  fLM <- lm(S1 ~ ns(Sb, df=2), data=dataB, weights=ifelse(dataB$Y==1, wtVCases, wtVControls))
  
  # normal density estimator for g(s0) using the placebo group in the immunogenicity set
  dataB <- subset(dataI, Z==0)
  gLM <- lm(S0 ~ 1, data=dataB, weights=ifelse(dataB$Y==1, wtPCases, wtPControls))
  
  # a single VE(s1) curve
  VEcurvePointEst <- sapply(s1grid, function(s1val){ pEstVE(s1val, data, dataI, fit1, fLM, gLM) })
  
  return(VEcurvePointEst)
}

# 'coverVEcurve' returns a vector of 0s and 1s indicating whether the true VE(s1) values on the 's1grid' are covered by pointwise Wald-type bootstrap CIs;
# the last value in the vector indicates coverage of the whole VE(s1) curve by the simultaneous Wald-type bootstrap CI
# 'data' is a data frame with variables Z, Sb, S0, S1, and Y
# 'nBoot' is the number of bootstrap iterations
coverVEcurve <- function(data, s1grid, trueVEcurve, nBoot){
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(S0) | !is.na(S1))
  
  # extract subsets of controls ('dataControls') and cases ('dataCases') to be used for resampling
  # in addition, within each treatment group in the immunogenicity set, delete cases to recover
  # the case:control ratio in the Phase 1 population
  dataControls <- subset(data, Y==0)
  nPControlsI <- NROW(dataPControlsI <- subset(dataI, Z==0 & Y==0))
  nVControlsI <- NROW(dataVControlsI <- subset(dataI, Z==1 & Y==0))  
  
  dataCases <- subset(data, Y==1)  
  nPCasesI <- NROW(dataPCasesI <- subset(dataI, Z==0 & Y==1))
  nVCasesI <- NROW(dataVCasesI <- subset(dataI, Z==1 & Y==1))
  
  nPControls <- NROW(subset(dataControls, Z==0))
  nVControls <- NROW(subset(dataControls, Z==1))
  nPCases <- NROW(subset(dataCases, Z==0))
  nVCases <- NROW(subset(dataCases, Z==1))
  
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
  fbw <- npcdensbw(S1 ~ Sb, data=dataVIcorrectRatio, cxkertype="epanechnikov", cykertype="epanechnikov")
  gbw <- npudensbw(~ S0, data=dataPIcorrectRatio, ckertype="epanechnikov")  
  
  group <- rep(1, NROW(subset(dataI, Z==0)))
  
  # weighted logistic regression model using the placebo group in the immunogenicity set
  fit1 <- tps(Y ~ S0, data=subset(dataI, Z==0), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
  
  # kernel density estimator for f(s1|Sb=sb) using the vaccine group in the immunogenicity set
  fhat <- npcdens(fbw)
  
  # kernel density estimator for g(s0) using the placebo group in the immunogenicity set
  ghat <- npudens(gbw)
  
  VEcurvePointEst <- sapply(s1grid, function(s1val){ estVE(s1val, data, dataI, fit1, fhat, ghat) })
  
  bSampleControls <- matrix(sample(1:nControls, nControls*nBoot, replace=TRUE), nrow=nControls, ncol=nBoot)
  bSampleCases <- matrix(sample(1:nCases, nCases*nBoot, replace=TRUE), nrow=nCases, ncol=nBoot)
  
  # 'bVEcurves' is a matrix with 'nBoot' columns each of which is a vector of bootstrap estimates of the VE curve on 's1grid'
  bVEcurves <- sapply(1:nBoot, function(i){
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,i],], dataCases[bSampleCases[,i],])
    # extract the bootstrapped immunogenicity set
    bdataI <- subset(bdata, !is.na(S0) | !is.na(S1))
    
    bdataControls <- subset(bdata, Y==0)
    nPControlsI <- NROW(bdataPControlsI <- subset(bdataI, Z==0 & Y==0))
    nVControlsI <- NROW(bdataVControlsI <- subset(bdataI, Z==1 & Y==0))
    
    bdataCases <- subset(bdata, Y==1)    
    nPCasesI <- NROW(bdataPCasesI <- subset(bdataI, Z==0 & Y==1))
    nVCasesI <- NROW(bdataVCasesI <- subset(bdataI, Z==1 & Y==1))
    
    nPControls <- NROW(subset(bdataControls, Z==0))
    nVControls <- NROW(subset(bdataControls, Z==1))
    nPCases <- NROW(subset(bdataCases, Z==0))
    nVCases <- NROW(subset(bdataCases, Z==1)) 
    
    # within each treatment group, calculate the number of cases in the bootstrapped immunogenicity set 
    # needed to achieve the correct case:control ratio
    nPCasesInew <- nPCases * nPControlsI / nPControls
    nVCasesInew <- nVCases * nVControlsI / nVControls
    
    # within each treatment group, sample as many cases in the bootstrapped immunogenicity set as needed 
    # to achieve the correct case:control ratio
    bdataPIcorrectRatio <- rbind(bdataPControlsI, bdataPCasesI[sample(1:nPCasesI, nPCasesInew),])
    bdataVIcorrectRatio <- rbind(bdataVControlsI, bdataVCasesI[sample(1:nVCasesI, nVCasesInew),])
    rm(bdataPControlsI); rm(bdataPCasesI); rm(bdataVControlsI); rm(bdataVCasesI)
    
    group <- rep(1, NROW(subset(bdataI, Z==0)))
    
    # weighted logistic regression model using the placebo group in the bootstrapped immunogenicity set
    fit1 <- tps(Y ~ S0, data=subset(bdataI, Z==0), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
    
    # kernel density estimator for f(s1|Sb=sb) using the vaccine group in the bootstrapped immunogenicity set
    bfbw <- npcdensbw(S1 ~ Sb, data=bdataVIcorrectRatio, bws=fbw, bandwidth.compute=FALSE)
    fhat <- npcdens(bfbw)
    
    # kernel density estimator for g(s0) using the placebo group in the bootstrapped immunogenicity set
    bgbw <- npudensbw(~ S0, data=bdataPIcorrectRatio, bws=gbw, bandwidth.compute=FALSE)    
    ghat <- npudens(bgbw)
    
    # a single bootstrap VE(s1) curve
    VEcurveBootEst <- sapply(s1grid, function(s1val){ estVE(s1val, bdata, bdataI, fit1, fhat, ghat) })
    return(VEcurveBootEst)
  })
  
  logRR <- log(1-VEcurvePointEst)
  bLogRRs <- log(1-bVEcurves)
  
  # bootstrap SE of log RR estimates
  bSE <- apply(bLogRRs, 1, sd, na.rm=TRUE)
  
  # pointwise confidence bounds for VE(s1)
  LB.VE <- 1 - exp(logRR + qnorm(0.975) * bSE)
  UB.VE <- 1 - exp(logRR - qnorm(0.975) * bSE)
  
  # indicator of the truth on 's1grid' being covered by pointwise CIs
  cover <- as.numeric(LB.VE<trueVEcurve & UB.VE>trueVEcurve)
  
  supAbsZ <- NULL
  for (j in 1:NCOL(bLogRRs)){
    Zstat <- abs((bLogRRs[,j]-logRR)/bSE)
    supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
  }
  qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
  
  LB.VE <- 1 - exp(logRR + qSupAbsZ * bSE)
  UB.VE <- 1 - exp(logRR - qSupAbsZ * bSE)
  
  # indicator of the truth on 's1grid' being covered by the simultaneous CI
  smCover <- as.numeric(all(LB.VE<trueVEcurve) && all(UB.VE>trueVEcurve))
  
  # the last value of 'cover' pertains to the simultaneous coverage
  cover <- c(cover, smCover)
  
  return(cover)
}

# 'pCoverVEcurve' returns a vector of 0s and 1s indicating whether the true VE(s1) values on the 's1grid' are covered by pointwise Wald-type bootstrap CIs;
# the last value in the vector indicates coverage of the whole VE(s1) curve by the simultaneous Wald-type bootstrap CI
# 'data' is a data frame with variables Z, Sb, S0, S1, and Y
# 'nBoot' is the number of bootstrap iterations
# parametric Gaussian density estimation is employed
pCoverVEcurve <- function(data, s1grid, trueVEcurve, nBoot){
  # extract the immunogenicity set
  dataI <- subset(data, !is.na(S0) | !is.na(S1))
  
  # calculate the sampling weights
  dataControls <- subset(data, Y==0)
  nPControls <- NROW(subset(dataControls, Z==0))
  wtPControls <- NROW(subset(dataControls, Z==0))/NROW(subset(dataI, Z==0 & Y==0))
  wtVControls <- NROW(subset(dataControls, Z==1))/NROW(subset(dataI, Z==1 & Y==0))
  
  dataCases <- subset(data, Y==1)
  nPCases <- NROW(subset(dataCases, Z==0))
  wtPCases <- NROW(subset(dataCases, Z==0))/NROW(subset(dataI, Z==0 & Y==1))
  wtVCases <- NROW(subset(dataCases, Z==1))/NROW(subset(dataI, Z==1 & Y==1))
  group <- rep(1, NROW(subset(dataI, Z==0)))
  
  # the overall numbers of controls and cases for resampling
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
  # weighted logistic regression model using the placebo group in the immunogenicity set
  fit1 <- tps(Y ~ S0, data=subset(dataI, Z==0), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
  
  # normal density estimator for f(s1|Sb=sb) using the vaccine group in the immunogenicity set with baseline markers
  # sampling weights are incorporated
  dataB <- subset(dataI, Z==1 & !is.na(Sb))
  fLM <- lm(S1 ~ ns(Sb, df=2), data=dataB, weights=ifelse(dataB$Y==1, wtVCases, wtVControls))
  
  # normal density estimator for g(s0) using the placebo group in the immunogenicity set
  dataB <- subset(dataI, Z==0)
  gLM <- lm(S0 ~ 1, data=dataB, weights=ifelse(dataB$Y==1, wtPCases, wtPControls))
  
  VEcurvePointEst <- sapply(s1grid, function(s1val){ pEstVE(s1val, data, dataI, fit1, fLM, gLM) })
  
  bSampleControls <- matrix(sample(1:nControls, nControls*nBoot, replace=TRUE), nrow=nControls, ncol=nBoot)
  bSampleCases <- matrix(sample(1:nCases, nCases*nBoot, replace=TRUE), nrow=nCases, ncol=nBoot)
  
  # 'bVEcurves' is a matrix with 'nBoot' columns each of which is a vector of bootstrap estimates of the VE curve on 's1grid'
  bVEcurves <- sapply(1:nBoot, function(i){
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,i],], dataCases[bSampleCases[,i],])
    # extract the immunogenicity set
    bdataI <- subset(bdata, !is.na(S0) | !is.na(S1))
    
    # calculate the sampling weights
    bdataControls <- subset(bdata, Y==0)
    nPControls <- NROW(subset(bdataControls, Z==0))
    wtPControls <- NROW(subset(bdataControls, Z==0))/NROW(subset(bdataI, Z==0 & Y==0))
    wtVControls <- NROW(subset(bdataControls, Z==1))/NROW(subset(bdataI, Z==1 & Y==0))
    
    bdataCases <- subset(bdata, Y==1)
    nPCases <- NROW(subset(bdataCases, Z==0))
    wtPCases <- NROW(subset(bdataCases, Z==0))/NROW(subset(bdataI, Z==0 & Y==1))
    wtVCases <- NROW(subset(bdataCases, Z==1))/NROW(subset(bdataI, Z==1 & Y==1))
    group <- rep(1, NROW(subset(bdataI, Z==0)))
    
    # weighted logistic regression model using the placebo group in the immunogenicity set
    fit1 <- tps(Y ~ S0, data=subset(bdataI, Z==0), nn0=nPControls, nn1=nPCases, group=group, method="PL", cohort=TRUE)
    
    # normal density estimator for f(s1|Sb=sb) using the vaccine group in the immunogenicity set with baseline markers
    # sampling weights are incorporated
    bdataB <- subset(bdataI, Z==1 & !is.na(Sb))
    fLM <- lm(S1 ~ ns(Sb, df=2), data=bdataB, weights=ifelse(bdataB$Y==1, wtVCases, wtVControls))
    
    # normal density estimator for g(s0) using the placebo group in the immunogenicity set
    bdataB <- subset(bdataI, Z==0)
    gLM <- lm(S0 ~ 1, data=bdataB, weights=ifelse(bdataB$Y==1, wtPCases, wtPControls))
    
    VEcurveBootEst <- sapply(s1grid, function(s1val){ pEstVE(s1val, bdata, bdataI, fit1, fLM, gLM) })
    return(VEcurveBootEst)
  })
  
  logRR <- log(1-VEcurvePointEst)
  bLogRRs <- log(1-bVEcurves)
  
  # bootstrap SE of log RR estimates
  bSE <- apply(bLogRRs, 1, sd, na.rm=TRUE)
  
  # pointwise confidence bounds for VE(s1)
  LB.VE <- 1 - exp(logRR + qnorm(0.975) * bSE)
  UB.VE <- 1 - exp(logRR - qnorm(0.975) * bSE)
  
  # indicator of the truth on 's1grid' being covered by pointwise CIs
  cover <- as.numeric(LB.VE<trueVEcurve & UB.VE>trueVEcurve)
  
  supAbsZ <- NULL
  for (j in 1:NCOL(bLogRRs)){
    Zstat <- abs((bLogRRs[,j]-logRR)/bSE)
    supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
  }
  qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
  
  LB.VE <- 1 - exp(logRR + qSupAbsZ * bSE)
  UB.VE <- 1 - exp(logRR - qSupAbsZ * bSE)
  
  # indicator of the truth on 's1grid' being covered by the simultaneous CI
  smCover <- as.numeric(all(LB.VE<trueVEcurve) && all(UB.VE>trueVEcurve))
  
  # the last value of 'cover' pertains to the simultaneous coverage
  cover <- c(cover, smCover)
  
  return(cover)
}

# 'getEstVE' performs 1 MC iteration, i.e., it generates the data-set and estimates the VE(s1) curve
getEstVE <- function(s1grid, n, beta, pi, truncateMarker, seed){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  return(VEcurve(data=data, s1grid=s1grid))
}

# 'getPestVE' performs 1 MC iteration, i.e., it generates the data-set and estimates the VE(s1) curve
# parametric Gaussian density estimation is employed
getPestVE <- function(s1grid, n, beta, pi, truncateMarker, seed){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  return(pVEcurve(data=data, s1grid=s1grid))
}

# 'getCoverVE' generates a data-set representing 1 MC iteration, computes the bootstrap SE based on 'nBoot'
# bootstrap iterations, and returns a vector of 0s and 1s indicating whether the truth is covered by
# the bootstrap Wald-type CI;
# the last value of the output vector pertains to the simultaneous CI
getCoverVE <- function(s1grid, trueVEcurve, n, beta, pi, truncateMarker, seed, nBoot){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  return(coverVEcurve(data=data, s1grid=s1grid, trueVEcurve=trueVEcurve, nBoot=nBoot))
}

# 'getPcoverVE' generates a data-set representing 1 MC iteration, computes the bootstrap SE based on 'nBoot'
# bootstrap iterations, and returns a vector of 0s and 1s indicating whether the truth is covered by
# the bootstrap Wald-type CI;
# the last value of the output vector pertains to the simultaneous CI;
# parametric Gaussian density estimation is employed
getPcoverVE <- function(s1grid, trueVEcurve, n, beta, pi, truncateMarker, seed, nBoot){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  return(pCoverVEcurve(data=data, s1grid=s1grid, trueVEcurve=trueVEcurve, nBoot=nBoot))
}