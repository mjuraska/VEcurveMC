library(MASS)
library(mvtnorm)
library(osDesign)
library(np)
library(splines)

logit <- function(p){
  return(log(p/(1-p)))
}

expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

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

# calculate the values of the probit risk model coefficients
getBeta <- function(beta2, beta3=NULL, meanS0, varS0, meanS1, varS1, probP, probV, rr, rr0, scenario, threshold){
  beta0 <- getBeta0(beta2, meanS0=meanS0, varS0=varS0, prob=probP, scenario=scenario, threshold=threshold)
  beta1start <- startBeta1(beta0, rr0)
  if (is.null(beta3)){
    beta3 <- getBeta3(beta0, beta1start, beta2, meanS1=meanS1, varS1=varS1, rr)  
  }
  beta1 <- getBeta1(beta0, beta2, beta3, meanS1=meanS1, varS1=varS1, probV, scenario=scenario, threshold=threshold)
  return(c(beta0, beta1, beta2, beta3))
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

# the true "risk curves" assuming the probit model for the probability of infection
# 's1' is a vector
# 'scenario' is one of "A1" and "A2"
trueRisk <- function(s1, beta0, beta1, beta2, beta3, meanS0, varS0, meanS1, varS1, covS0S1, scenario="A1", threshold=NULL){
  txRisk <- pnorm(beta0+beta1+beta3*s1)
  plaRisk <- sapply(s1, function(s1Value){ risk0(s1Value, beta0, beta2, meanS0, varS0, meanS1, varS1, covS0S1, scenario=scenario, threshold=threshold) })
  return(list(plaRisk=plaRisk, txRisk=txRisk))
}

# 'n' is the total sample size
# 'beta' is a vector of probit risk model coefficients
# 'pi' is a sampling probability for sampling a subcohort of controls with biomarker measurements
# 'truncateMarker' is TRUE indicates truncated marker distributions
getData <- function(n, beta, pi, sigma12, truncateMarker, seed){
  set.seed(seed)
  
  Z <- rep(0:1, each=n/2)
  S <- mvrnorm(n, mu=c(2,2,3), Sigma=matrix(c(1,sigma12,0.7,sigma12,1,0.7,0.7,0.7,1), nrow=3))
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

# 'pEstRisk' is a version of 'pEstVE' and returns a numeric vector with the estimated P(Y(0)=1|S(1)=s1) and P(Y(1)=1|S(1)=s1) in this order
# s1 is a scalar
pEstRisk <- function(s1, data, dataI, tpsFit, lmFit1, lmFit2){
  riskPvalue <- pRiskP(s1, data, tpsFit, lmFit1, lmFit2)
  
  if (is.null(riskPvalue)){ 
    return(c(NA,NA))
  } else {
    return(c(riskPvalue, riskV(s1, data, dataI)))
  }
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
  fLM <- lm(S1 ~ ns(Sb, df=3), data=dataB, weights=ifelse(dataB$Y==1, wtVCases, wtVControls))
  
  # normal density estimator for g(s0) using the placebo group in the immunogenicity set
  dataB <- subset(dataI, Z==0)
  gLM <- lm(S0 ~ 1, data=dataB, weights=ifelse(dataB$Y==1, wtPCases, wtPControls))
  
  # a single VE(s1) curve
  VEcurvePointEst <- sapply(s1grid, function(s1val){ pEstVE(s1val, data, dataI, fit1, fLM, gLM) })
  
  return(VEcurvePointEst)
}

# 'pRiskCurve' is a version of 'pVEcurve' outputting the "risk curve" on 's1grid' in each treatment arm instead of the VE curve
# it allows evaluation of an arbitrary mCEP curve estimand, i.e., an arbitrary contrast function
# 'data' is a data frame with variables Z, Sb, S0, S1, and Y
# the output is a list with 2 components being the 2 risk curves
pRiskCurve <- function(data, s1grid){
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
  fLM <- lm(S1 ~ ns(Sb, df=3), data=dataB, weights=ifelse(dataB$Y==1, wtVCases, wtVControls))
  
  # normal density estimator for g(s0) using the placebo group in the immunogenicity set
  dataB <- subset(dataI, Z==0)
  gLM <- lm(S0 ~ 1, data=dataB, weights=ifelse(dataB$Y==1, wtPCases, wtPControls))
  
  # a matrix with the rows being the "risk curves" on 's1grid' in the placebo and treatment arm
  riskCurvePointEst <- sapply(s1grid, function(s1val){ pEstRisk(s1val, data, dataI, fit1, fLM, gLM) })
  
  return(list(plaRiskCurvePointEst=riskCurvePointEst[1,], txRiskCurvePointEst=riskCurvePointEst[2,]))
}

# 'VEcurve2' is a version of 'VEcurve' that, besides the estimated VE(s1) curve evaluated on the grid of the 's1grid' values, returns other
# quantities needed for the bootstrap
# 'data' is a data frame with variables Z, Sb, S0, S1, and Y
VEcurve2 <- function(data, s1grid){
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
  
  return(list(VEcurvePointEst=VEcurvePointEst, fbw=fbw, gbw=gbw))
}

bVEcurve2 <- function(data, s1grid, nBoot, fbw, gbw){
  dataControls <- subset(data, Y==0)
  dataCases <- subset(data, Y==1)
  
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
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
  
  return(bVEcurves)
}

pbVEcurve <- function(data, s1grid, nBoot){
  dataControls <- subset(data, Y==0)
  dataCases <- subset(data, Y==1)
  
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
  bSampleControls <- matrix(sample(1:nControls, nControls*nBoot, replace=TRUE), nrow=nControls, ncol=nBoot)
  bSampleCases <- matrix(sample(1:nCases, nCases*nBoot, replace=TRUE), nrow=nCases, ncol=nBoot)
  
  # 'bVEcurves' is a matrix with 'nBoot' columns each of which is a vector of bootstrap estimates of the VE curve on 's1grid'
  bVEcurves <- sapply(1:nBoot, function(i){
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,i],], dataCases[bSampleCases[,i],])
    # extract the bootstrapped immunogenicity set
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
    fLM <- lm(S1 ~ ns(Sb, df=3), data=bdataB, weights=ifelse(bdataB$Y==1, wtVCases, wtVControls))
    
    # normal density estimator for g(s0) using the placebo group in the immunogenicity set
    bdataB <- subset(bdataI, Z==0)
    gLM <- lm(S0 ~ 1, data=bdataB, weights=ifelse(bdataB$Y==1, wtPCases, wtPControls))
    
    VEcurveBootEst <- sapply(s1grid, function(s1val){ pEstVE(s1val, bdata, bdataI, fit1, fLM, gLM) })
    return(VEcurveBootEst)
  })
  
  return(bVEcurves)
}

# a version of 'pbVEcurve' outputting bootstrapped "risk curves"
pbRiskCurves <- function(data, s1grid, nBoot){
  dataControls <- subset(data, Y==0)
  dataCases <- subset(data, Y==1)
  
  nControls <- NROW(dataControls)
  nCases <- NROW(dataCases)
  
  bSampleControls <- matrix(sample(1:nControls, nControls*nBoot, replace=TRUE), nrow=nControls, ncol=nBoot)
  bSampleCases <- matrix(sample(1:nCases, nCases*nBoot, replace=TRUE), nrow=nCases, ncol=nBoot)
  
  # 'bVEcurves' is a matrix with 'nBoot' columns each of which is a vector of bootstrap estimates of the VE curve on 's1grid'
  bRiskCurves <- lapply(1:nBoot, function(i){
    # create a bootstrap sample
    bdata <- rbind(dataControls[bSampleControls[,i],], dataCases[bSampleCases[,i],])
    # extract the bootstrapped immunogenicity set
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
    fLM <- lm(S1 ~ ns(Sb, df=3), data=bdataB, weights=ifelse(bdataB$Y==1, wtVCases, wtVControls))
    
    # normal density estimator for g(s0) using the placebo group in the immunogenicity set
    bdataB <- subset(bdataI, Z==0)
    gLM <- lm(S0 ~ 1, data=bdataB, weights=ifelse(bdataB$Y==1, wtPCases, wtPControls))
    
    riskCurvesBootEst <- sapply(s1grid, function(s1val){ pEstRisk(s1val, bdata, bdataI, fit1, fLM, gLM) })
    return(list(plaRiskCurveBootEst=riskCurvesBootEst[1,], txRiskCurveBootEst=riskCurvesBootEst[2,]))
  })
  
  return(bRiskCurves)
}


# 'coverVEcurve' returns indicators of coverage of the true VE(s1) curve by both pointwise and simultaneous CIs based on 1 MC iteration
# assumes that all 3 input arguments are evaluated on the same grid of s1 values
coverVEcurve <- function(VEcurvePointEst, bVEcurves, trueVEcurve){
  logRR <- log(1-VEcurvePointEst)
  bLogRRs <- log(1-bVEcurves)
  
  # bootstrap SE of log RR estimates
  bSE <- apply(bLogRRs, 1, sd, na.rm=TRUE)
  
  # pointwise confidence bounds for VE(s1)
  LB.VE <- 1 - exp(logRR + qnorm(0.975) * bSE)
  UB.VE <- 1 - exp(logRR - qnorm(0.975) * bSE)
  
  # indicator of the truth on 's1grid' being covered by pointwise CIs
  coverInd <- as.numeric(LB.VE<trueVEcurve & UB.VE>trueVEcurve)
  
  supAbsZ <- NULL
  for (j in 1:NCOL(bLogRRs)){
    Zstat <- abs((bLogRRs[,j]-logRR)/bSE)
    supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
  }
  qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
  
  LB.VE <- 1 - exp(logRR + qSupAbsZ * bSE)
  UB.VE <- 1 - exp(logRR - qSupAbsZ * bSE)
  
  # indicator of the truth on 's1grid' being covered by the simultaneous CI
  smCoverInd <- as.numeric(all(LB.VE<trueVEcurve) && all(UB.VE>trueVEcurve))
  
  return(list(coverInd=coverInd, smCoverInd=smCoverInd))
}

# 'coverMCEPcurve' returns indicators of coverage of the true MCEP(s1) curve by both pointwise and simultaneous CIs based on 1 MC iteration
# assumes that all 3 input arguments are evaluated on the same grid of s1 values
# 'normTransform' is one of "identity" and "logit"
coverMCEPcurve <- function(riskCurvePointEst, bRiskCurves, trueRiskCurves, normTransform){
  # true MCEP curve
  trueMCEP <- trueRiskCurves$plaRisk - trueRiskCurves$txRisk
  
  # cbind all bootstrap estimates to make easier to transform
  plaRiskCurveBootEst <- sapply(bRiskCurves, "[[", "plaRiskCurveBootEst")
  txRiskCurveBootEst <- sapply(bRiskCurves, "[[", "txRiskCurveBootEst")
  
  if (normTransform=="identity"){
    # transformed estimated MCEP curve
    tMCEP <- riskCurvePointEst$plaRiskCurvePointEst - riskCurvePointEst$txRiskCurvePointEst
    # transformed bootstrapped MCEP curves
    tbMCEP <- plaRiskCurveBootEst - txRiskCurveBootEst
    
    # bootstrap SE of tMCEP estimates
    bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)
    
    # pointwise confidence bounds for MCEP(s1)
    LB.MCEP <- tMCEP - qnorm(0.975) * bSE
    UB.MCEP <- tMCEP + qnorm(0.975) * bSE
    
    # indicator of the truth on 's1grid' being covered by pointwise CIs
    coverInd <- as.numeric(LB.MCEP<trueMCEP & UB.MCEP>trueMCEP)
    
    supAbsZ <- NULL
    for (j in 1:NCOL(tbMCEP)){
      Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
      supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
    }
    qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
    
    LB.MCEP <- tMCEP - qSupAbsZ * bSE
    UB.MCEP <- tMCEP + qSupAbsZ * bSE
    
    # indicator of the truth on 's1grid' being covered by the simultaneous CI
    smCoverInd <- as.numeric(all(LB.MCEP<trueMCEP) && all(UB.MCEP>trueMCEP))
  }
  
  if (normTransform=="logit"){
    # transformed MCEP curve
    tMCEP <- logit(((riskCurvePointEst$plaRiskCurvePointEst - riskCurvePointEst$txRiskCurvePointEst) + 1)/2)
    # transformed bootstrapped MCEP curves
    tbMCEP <- logit(((plaRiskCurveBootEst - txRiskCurveBootEst) + 1)/2)
    
    # bootstrap SE of tMCEP estimates
    bSE <- apply(tbMCEP, 1, sd, na.rm=TRUE)
    
    # pointwise confidence bounds for MCEP(s1)
    LB.MCEP <- 2*expit(tMCEP - qnorm(0.975) * bSE) - 1
    UB.MCEP <- 2*expit(tMCEP + qnorm(0.975) * bSE) - 1
    
    # indicator of the truth on 's1grid' being covered by pointwise CIs
    coverInd <- as.numeric(LB.MCEP<trueMCEP & UB.MCEP>trueMCEP)
    
    supAbsZ <- NULL
    for (j in 1:NCOL(tbMCEP)){
      Zstat <- abs((tbMCEP[,j]-tMCEP)/bSE)
      supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat))))
    }
    qSupAbsZ <- quantile(supAbsZ, probs=0.95, na.rm=TRUE)
    
    LB.MCEP <- 2*expit(tMCEP - qSupAbsZ * bSE) - 1
    UB.MCEP <- 2*expit(tMCEP + qSupAbsZ * bSE) - 1
    
    # indicator of the truth on 's1grid' being covered by the simultaneous CI
    smCoverInd <- as.numeric(all(LB.MCEP<trueMCEP) && all(UB.MCEP>trueMCEP))
  }
  
  return(list(coverInd=coverInd, smCoverInd=smCoverInd))
}

# returns a 2-sided p-value of {H0: VE(s1)=VE for all s1} against {H1: non-H0}
testConstancy <- function(VEcurve, bVEcurves){
  logRR <- log(1-VEcurve)
  bLogRRs <- log(1-bVEcurves)
  
  # bootstrap SE of log RR estimates
  bSE <- apply(bLogRRs, 1, sd, na.rm=TRUE)
  
  # calculate the supremum statistic for each bootstrap sample
  supAbsZ <- NULL
  for (i in 1:NCOL(bLogRRs)){ 
    Zstat <- abs((bLogRRs[,i]-logRR)/bSE)
    supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat)))) 
  }
  
  # 2-sided p-value
  p <- uniroot(f, 0:1, logRR=logRR, bSE=bSE, supAbsZ=supAbsZ)$root
  return(p)
}

f <- function(alpha, logRR, bSE, supAbsZ){
  qSupAbsZ <- ifelse(alpha==1, 0, ifelse(alpha==0, Inf, quantile(supAbsZ, probs=1-alpha, na.rm=TRUE)))
  minUB <- min(logRR + qSupAbsZ * bSE, na.rm=TRUE)
  maxLB <- max(logRR - qSupAbsZ * bSE, na.rm=TRUE)
  return(minUB - maxLB)
}

# returns an indicator of rejecting {H0: VE(s1)=VE for all s1} in favor of {H1: non-H0}
# 'overallVE' is the point estimate of 1 - P(Y=1|Z=1)/P(Y=1|Z=0)
testConstancy2 <- function(VEcurve, bVEcurves, overallVE, alpha=0.05){
  oLogRR <- log(1-overallVE)
  logRR <- log(1-VEcurve)
  bLogRRs <- log(1-bVEcurves)
  
  # bootstrap SE of log RR estimates
  bSE <- apply(bLogRRs, 1, sd, na.rm=TRUE)
  
  # calculate the supremum statistic for each bootstrap sample
  supAbsZ <- NULL
  for (i in 1:NCOL(bLogRRs)){ 
    Zstat <- abs((bLogRRs[,i]-logRR)/bSE)
    supAbsZ <- c(supAbsZ, max(Zstat, na.rm=!all(is.na(Zstat)))) 
  }
  # critical value
  qSupAbsZ <- quantile(supAbsZ, probs=1-alpha, na.rm=TRUE)
  
  testStat <- max(abs(logRR-oLogRR)/bSE, na.rm=TRUE)
  
  return(as.numeric(testStat > qSupAbsZ))
}

# returns a 2-sided p-value of {H0: VE_x(s1)=VE_y(s1) for all s1 in [s_l,s_u]} against {H1: non-H0}, 
# where VE_x and VE_y represent two independent populations (e.g., two trials in different geographic regions);
# the code assumes that all VE curves were estimated on a common grid of marker values
testTrialEquality <- function(VEcurve1, VEcurve2, bVEcurves1, bVEcurves2, s1grid=NULL, markerInterval=NULL){
  if (!is.null(markerInterval)){
    if (!is.null(s1grid)){
      VEcurve1 <- VEcurve1[s1grid>=markerInterval[1] & s1grid<=markerInterval[2]]
      VEcurve2 <- VEcurve2[s1grid>=markerInterval[1] & s1grid<=markerInterval[2]]
      bVEcurves1 <- bVEcurves1[s1grid>=markerInterval[1] & s1grid<=markerInterval[2],]
      bVEcurves2 <- bVEcurves2[s1grid>=markerInterval[1] & s1grid<=markerInterval[2],]  
    } else {
      stop("The argument 's1grid' is missing.")
    }
  }
  
  # difference in the log relative risk estimates
  dlogRR <- log(1-VEcurve1) - log(1-VEcurve2)
  dbLogRR <- log(1-bVEcurves1) - log(1-bVEcurves2)
  
  # bootstrap SE of the difference in the log RR estimates calculated, due to independence, as the square root of the sum of the population-specific sample variances
  bSE <- sqrt(apply(log(1-bVEcurves1), 1, var, na.rm=TRUE) + apply(log(1-bVEcurves2), 1, var, na.rm=TRUE))
  
  # calculate the supremum statistic for each bootstrap sample
  supAbsZ <- NULL
  for (i in 1:NCOL(dbLogRR)){ 
    Z <- abs((dbLogRR[,i]-dlogRR)/bSE)
    supAbsZ <- c(supAbsZ, max(Z, na.rm=!all(is.na(Z)))) 
  }
  
  # one of these two equations may not have a solution but at least one of them will always have a solution
  alpha1 <- try(uniroot(fInf, 0:1, dlogRR=dlogRR, bSE=bSE, supAbsZ=supAbsZ)$root, silent=TRUE)
  alpha2 <- try(uniroot(fSup, 0:1, dlogRR=dlogRR, bSE=bSE, supAbsZ=supAbsZ)$root, silent=TRUE)
  
  if (inherits(alpha1, 'try-error')){
    # then 'alpha1' does not exist
    alpha1 <- NA
  } else {
    # now we can plug 'alpha1' into 'fSup'
    alpha1 <- ifelse(fSup(alpha1, dlogRR, bSE, supAbsZ)<=0, alpha1, NA)  
  }
  
  if (inherits(alpha2, 'try-error')){
    # then 'alpha2' does not exist
    alpha2 <- NA
  } else {
    # now we can plug 'alpha2' into 'fInf'
    alpha2 <- ifelse(fInf(alpha2, dlogRR, bSE, supAbsZ)>=0, alpha2, NA)
  }
  
  # 2-sided p-value
  return(min(alpha1, alpha2, na.rm=TRUE))
}

fInf <- function(alpha, dlogRR, bSE, supAbsZ){
  qSupAbsZ <- ifelse(alpha==1, 0, ifelse(alpha==0, Inf, quantile(supAbsZ, probs=1-alpha, na.rm=TRUE)))
  return(min(dlogRR + qSupAbsZ * bSE, na.rm=TRUE))
}

fSup <- function(alpha, dlogRR, bSE, supAbsZ){
  qSupAbsZ <- ifelse(alpha==1, 0, ifelse(alpha==0, Inf, quantile(supAbsZ, probs=1-alpha, na.rm=TRUE)))
  return(max(dlogRR - qSupAbsZ * bSE, na.rm=TRUE))
}

testTrialEquality2 <- function(VEcurve1, VEcurve2, bVEcurves1, bVEcurves2, s1grid=NULL, markerInterval=NULL, alpha=0.05){
  if (!is.null(markerInterval)){
    if (!is.null(s1grid)){
      VEcurve1 <- VEcurve1[s1grid>=markerInterval[1] & s1grid<=markerInterval[2]]
      VEcurve2 <- VEcurve2[s1grid>=markerInterval[1] & s1grid<=markerInterval[2]]
      bVEcurves1 <- bVEcurves1[s1grid>=markerInterval[1] & s1grid<=markerInterval[2],]
      bVEcurves2 <- bVEcurves2[s1grid>=markerInterval[1] & s1grid<=markerInterval[2],]  
    } else {
      stop("The argument 's1grid' is missing.")
    }
  }
  
  # difference in the log relative risk estimates
  dlogRR <- log(1-VEcurve1) - log(1-VEcurve2)
  dbLogRR <- log(1-bVEcurves1) - log(1-bVEcurves2)
  
  # bootstrap SE of the difference in the log RR estimates calculated, due to independence, as the square root of the sum of the population-specific sample variances
  bSE <- sqrt(apply(log(1-bVEcurves1), 1, var, na.rm=TRUE) + apply(log(1-bVEcurves2), 1, var, na.rm=TRUE))
  
  # calculate the supremum statistic for each bootstrap sample
  supAbsZ <- NULL
  for (i in 1:NCOL(dbLogRR)){ 
    Z <- abs((dbLogRR[,i]-dlogRR)/bSE)
    supAbsZ <- c(supAbsZ, max(Z, na.rm=!all(is.na(Z)))) 
  }
  # critical value
  qSupAbsZ <- quantile(supAbsZ, probs=1-alpha, na.rm=TRUE)
  
  testStat <- max(abs(dlogRR)/bSE, na.rm=TRUE)
  
  return(as.numeric(testStat > qSupAbsZ))
}


# 'inferenceVEcurve' returns a list with the following components:
# coverInd                - a vector of 0s and 1s indicating whether the true VE(s1) values on the 's1grid' are covered by the pointwise bootstrap Wald-type CI
# smCoverInd              - a 0 or 1 indicating coverage of the entire VE(s1) curve on the 's1grid' by the simultaneous bootstrap Wald-type CI
# pSizeTestConstancy      - a 2-sided p-value of {H0: VE(s1)=VE for all s1} against {H1: non-H0} in a scenario that satisfies H0
# pPowerTestConstancy     - a 2-sided p-value of {H0: VE(s1)=VE for all s1} against {H1: non-H0} in the same scenario as that used for coverage
# pSizeTestTrialEquality  - a 2-sided p-value of {H0: VE_x(s1)=VE_y(s1) for all s1 in [s_l,s_u]} against {H1: non-H0} in a scenario that satisfies H0
# pPowerTestTrialEquality - a 2-sided p-value of {H0: VE_x(s1)=VE_y(s1) for all s1 in [s_l,s_u]} against {H1: non-H0} in a scenario that satisfies H1
inferenceVEcurve <- function(data, dataSizeT1, dataSizeT3, dataPowerT3, s1grid, trueVEcurve, nBoot){
  VEcurveEst <- VEcurve2(data, s1grid)
  
  bVEcurves <- bVEcurve2(data, s1grid, nBoot, VEcurveEst$fbw, VEcurveEst$gbw)
  
  # a list with 'coverInd' and 'smCoverInd'
  cover <- coverVEcurve(VEcurveEst$VEcurvePointEst, bVEcurves, trueVEcurve)
  
  # compute the p-value from the test of constancy under H0
  VEcurveEstSizeT1 <- VEcurve2(dataSizeT1, s1grid)
  bVEcurvesSizeT1 <- bVEcurve2(dataSizeT1, s1grid, nBoot, VEcurveEstSizeT1$fbw, VEcurveEstSizeT1$gbw)
  pSizeTestConstancy <- testConstancy(VEcurveEstSizeT1$VEcurvePointEst, bVEcurvesSizeT1)
  
  # compute the p-value from the test of constancy under H1
  pPowerTestConstancy <- testConstancy(VEcurveEst$VEcurvePointEst, bVEcurves)
  
  # compute the p-value from the test of equality in two populations under H0
  VEcurveEstSizeT3 <- VEcurve2(dataSizeT3, s1grid)
  bVEcurvesSizeT3 <- bVEcurve2(dataSizeT3, s1grid, nBoot, VEcurveEstSizeT3$fbw, VEcurveEstSizeT3$gbw)
  pSizeTestTrialEquality <- testTrialEquality(VEcurveEst$VEcurvePointEst, VEcurveEstSizeT3$VEcurvePointEst, bVEcurves, bVEcurvesSizeT3)
  
  # compute the p-value from the test of equality in two populations under H1
  VEcurveEstPowerT3 <- VEcurve2(dataPowerT3, s1grid)
  bVEcurvesPowerT3 <- bVEcurve2(dataPowerT3, s1grid, nBoot, VEcurveEstPowerT3$fbw, VEcurveEstPowerT3$gbw)
  pPowerTestTrialEquality <- testTrialEquality(VEcurveEst$VEcurvePointEst, VEcurveEstPowerT3$VEcurvePointEst, bVEcurves, bVEcurvesPowerT3)
  
  return(list(coverInd=cover$coverInd, smCoverInd=cover$smCoverInd, pSizeTestConstancy=pSizeTestConstancy, pPowerTestConstancy=pPowerTestConstancy,
              pSizeTestTrialEquality=pSizeTestTrialEquality, pPowerTestTrialEquality=pPowerTestTrialEquality))
}

inferenceVEcurve2 <- function(data, dataSizeT1, dataSizeT3, dataPowerT3, s1grid, trueVEcurve, nBoot){
  VEcurveEst <- VEcurve2(data, s1grid)
  
  bVEcurves <- bVEcurve2(data, s1grid, nBoot, VEcurveEst$fbw, VEcurveEst$gbw)
  
  # a list with 'coverInd' and 'smCoverInd'
  cover <- coverVEcurve(VEcurveEst$VEcurvePointEst, bVEcurves, trueVEcurve)
  
  # compute the indicator of rejection of the null hypothesis in the test of H^1_0 under validity of the null
  VEcurveEstSizeT1 <- VEcurve2(dataSizeT1, s1grid)
  bVEcurvesSizeT1 <- bVEcurve2(dataSizeT1, s1grid, nBoot, VEcurveEstSizeT1$fbw, VEcurveEstSizeT1$gbw)
  
  fit <- glm(Y ~ Z, data=data, family=binomial)
  prob <- predict(fit, newdata=data.frame(Z=0:1), type="response")
  overallVE <- 1 - prob[2]/prob[1]
  rejectIndSizeTestH10 <- testConstancy2(VEcurveEstSizeT1$VEcurvePointEst, bVEcurvesSizeT1, overallVE)
  
  # compute the indicator of rejection of the null hypothesis in the test of H^2_0 under validity of the null
  rejectIndSizeTestH20 <- testConstancy2(VEcurveEstSizeT1$VEcurvePointEst, bVEcurvesSizeT1, 0.5)
  
  # compute the indicator of rejection of the null hypothesis in the test of H^1_0 under an alternative
  rejectIndPowerTestH10 <- testConstancy2(VEcurveEst$VEcurvePointEst, bVEcurves, overallVE)
  
  # compute the indicator of rejection of the null hypothesis in the test of H^2_0 under an alternative
  rejectIndPowerTestH20 <- testConstancy2(VEcurveEst$VEcurvePointEst, bVEcurves, 0.5)
  
  # compute the indicator of rejection of the null hypothesis in the test of H^4_0 under validity of the null
  VEcurveEstSizeT3 <- VEcurve2(dataSizeT3, s1grid)
  bVEcurvesSizeT3 <- bVEcurve2(dataSizeT3, s1grid, nBoot, VEcurveEstSizeT3$fbw, VEcurveEstSizeT3$gbw)
  rejectIndSizeTestH40 <- testTrialEquality2(VEcurveEst$VEcurvePointEst, VEcurveEstSizeT3$VEcurvePointEst, bVEcurves, bVEcurvesSizeT3)
  
  # compute the indicator of rejection of the null hypothesis in the test of H^4_0 under an alternative
  VEcurveEstPowerT3 <- VEcurve2(dataPowerT3, s1grid)
  bVEcurvesPowerT3 <- bVEcurve2(dataPowerT3, s1grid, nBoot, VEcurveEstPowerT3$fbw, VEcurveEstPowerT3$gbw)
  rejectIndPowerTestH40 <- testTrialEquality2(VEcurveEst$VEcurvePointEst, VEcurveEstPowerT3$VEcurvePointEst, bVEcurves, bVEcurvesPowerT3)
  
  return(list(coverInd=cover$coverInd, smCoverInd=cover$smCoverInd, rejectIndSizeTestH10=rejectIndSizeTestH10, rejectIndPowerTestH10=rejectIndPowerTestH10, 
              rejectIndSizeTestH20=rejectIndSizeTestH20, rejectIndPowerTestH20=rejectIndPowerTestH20, rejectIndSizeTestH40=rejectIndSizeTestH40, 
              rejectIndPowerTestH40=rejectIndPowerTestH40))
}

# parametric Gaussian density estimation is employed
pInferenceVEcurve <- function(data, dataSizeT1, dataSizeT3, dataPowerT3, s1grid, trueVEcurve, nBoot){
  VEcurveEst <- pVEcurve(data, s1grid)
  
  bVEcurves <- pbVEcurve(data, s1grid, nBoot)
  
  # a list with 'coverInd' and 'smCoverInd'
  cover <- coverVEcurve(VEcurveEst, bVEcurves, trueVEcurve)
  
  # compute the p-value from the test of constancy under H0
  VEcurveEstSizeT1 <- pVEcurve(dataSizeT1, s1grid)
  bVEcurvesSizeT1 <- pbVEcurve(dataSizeT1, s1grid, nBoot)
  pSizeTestConstancy <- testConstancy(VEcurveEstSizeT1, bVEcurvesSizeT1)
  
  # compute the p-value from the test of constancy under H1
  pPowerTestConstancy <- testConstancy(VEcurveEst, bVEcurves)
  
  # compute the p-value from the test of equality in two populations under H0
  VEcurveEstSizeT3 <- pVEcurve(dataSizeT3, s1grid)
  bVEcurvesSizeT3 <- pbVEcurve(dataSizeT3, s1grid, nBoot)
  pSizeTestTrialEquality <- testTrialEquality(VEcurveEst, VEcurveEstSizeT3, bVEcurves, bVEcurvesSizeT3)
  
  # compute the p-value from the test of equality in two populations under H1
  VEcurveEstPowerT3 <- pVEcurve(dataPowerT3, s1grid)
  bVEcurvesPowerT3 <- pbVEcurve(dataPowerT3, s1grid, nBoot)
  pPowerTestTrialEquality <- testTrialEquality(VEcurveEst, VEcurveEstPowerT3, bVEcurves, bVEcurvesPowerT3)
  
  return(list(coverInd=cover$coverInd, smCoverInd=cover$smCoverInd, pSizeTestConstancy=pSizeTestConstancy, pPowerTestConstancy=pPowerTestConstancy,
              pSizeTestTrialEquality=pSizeTestTrialEquality, pPowerTestTrialEquality=pPowerTestTrialEquality))
}

# parametric Gaussian density estimation is employed
pInferenceVEcurve2 <- function(data, dataSizeT1, dataSizeT3, dataPowerT3, s1grid, trueVEcurve, nBoot){
  VEcurveEst <- pVEcurve(data, s1grid)
  
  bVEcurves <- pbVEcurve(data, s1grid, nBoot)
  
  # a list with 'coverInd' and 'smCoverInd'
  cover <- coverVEcurve(VEcurveEst, bVEcurves, trueVEcurve)
  
  # compute the p-value from the test of constancy under H0
  VEcurveEstSizeT1 <- pVEcurve(dataSizeT1, s1grid)
  bVEcurvesSizeT1 <- pbVEcurve(dataSizeT1, s1grid, nBoot)
  fit <- glm(Y ~ Z, data=data, family=binomial)
  prob <- predict(fit, newdata=data.frame(Z=0:1), type="response")
  #overallVE <- 1 - prob[2]/prob[1]
  overallVE <- 0.5
  rejectIndSizeTestConstancy <- testConstancy2(VEcurveEstSizeT1, bVEcurvesSizeT1, overallVE)
  
  # compute the p-value from the test of constancy under H1
  rejectIndPowerTestConstancy <- testConstancy2(VEcurveEst, bVEcurves, overallVE)
  
  # compute the p-value from the test of equality in two populations under H0
  VEcurveEstSizeT3 <- pVEcurve(dataSizeT3, s1grid)
  bVEcurvesSizeT3 <- pbVEcurve(dataSizeT3, s1grid, nBoot)
  rejectIndSizeTestTrialEquality <- testTrialEquality2(VEcurveEst, VEcurveEstSizeT3, bVEcurves, bVEcurvesSizeT3)
  
  # compute the p-value from the test of equality in two populations under H1
  VEcurveEstPowerT3 <- pVEcurve(dataPowerT3, s1grid)
  bVEcurvesPowerT3 <- pbVEcurve(dataPowerT3, s1grid, nBoot)
  rejectIndPowerTestTrialEquality <- testTrialEquality2(VEcurveEst, VEcurveEstPowerT3, bVEcurves, bVEcurvesPowerT3)
  
  return(list(coverInd=cover$coverInd, smCoverInd=cover$smCoverInd, rejectIndSizeTestConstancy=rejectIndSizeTestConstancy, rejectIndPowerTestConstancy=rejectIndPowerTestConstancy,
              rejectIndSizeTestTrialEquality=rejectIndSizeTestTrialEquality, rejectIndPowerTestTrialEquality=rejectIndPowerTestTrialEquality))
}

# a version of 'pInferenceVEcurve2' for different contrasts (currently implemented for additive difference only)
# two normalizing transformations are considered: identity and logit((x+1)/2)
pInferenceMCEPcurve2 <- function(data, dataSizeT1, dataSizeT3, dataPowerT3, s1grid, trueRiskCurves, nBoot, contrast="additive"){
  # a list with components 'plaRiskCurvePointEst' and 'txRiskCurvePointEst'
  riskCurveEst <- pRiskCurve(data, s1grid)
  
  # a list of length 'nBoot' each component of which is a list with components 'plaRiskCurveBootEst' and 'txRiskCurveBootEst'
  bRiskCurves <- pbRiskCurves(data, s1grid, nBoot)
  
  # a list with 'coverInd' and 'smCoverInd'
  coverIdTransform <- coverMCEPcurve(riskCurveEst, bRiskCurves, trueRiskCurves, normTransform="identity")
  coverLogitTransform <- coverMCEPcurve(riskCurveEst, bRiskCurves, trueRiskCurves, normTransform="logit")
  
  # # compute the p-value from the test of constancy under H0
  # VEcurveEstSizeT1 <- pVEcurve(dataSizeT1, s1grid)
  # bVEcurvesSizeT1 <- pbVEcurve(dataSizeT1, s1grid, nBoot)
  # fit <- glm(Y ~ Z, data=data, family=binomial)
  # prob <- predict(fit, newdata=data.frame(Z=0:1), type="response")
  # #overallVE <- 1 - prob[2]/prob[1]
  # overallVE <- 0.5
  # rejectIndSizeTestConstancy <- testConstancy2(VEcurveEstSizeT1, bVEcurvesSizeT1, overallVE)
  # 
  # # compute the p-value from the test of constancy under H1
  # rejectIndPowerTestConstancy <- testConstancy2(VEcurveEst, bVEcurves, overallVE)
  # 
  # # compute the p-value from the test of equality in two populations under H0
  # VEcurveEstSizeT3 <- pVEcurve(dataSizeT3, s1grid)
  # bVEcurvesSizeT3 <- pbVEcurve(dataSizeT3, s1grid, nBoot)
  # rejectIndSizeTestTrialEquality <- testTrialEquality2(VEcurveEst, VEcurveEstSizeT3, bVEcurves, bVEcurvesSizeT3)
  # 
  # # compute the p-value from the test of equality in two populations under H1
  # VEcurveEstPowerT3 <- pVEcurve(dataPowerT3, s1grid)
  # bVEcurvesPowerT3 <- pbVEcurve(dataPowerT3, s1grid, nBoot)
  # rejectIndPowerTestTrialEquality <- testTrialEquality2(VEcurveEst, VEcurveEstPowerT3, bVEcurves, bVEcurvesPowerT3)
  # 
  # return(list(coverInd=cover$coverInd, smCoverInd=cover$smCoverInd, rejectIndSizeTestConstancy=rejectIndSizeTestConstancy, rejectIndPowerTestConstancy=rejectIndPowerTestConstancy,
  #             rejectIndSizeTestTrialEquality=rejectIndSizeTestTrialEquality, rejectIndPowerTestTrialEquality=rejectIndPowerTestTrialEquality))
  return(list(coverIndIdTransform=coverIdTransform$coverInd, smCoverIndIdTransform=coverIdTransform$smCoverInd, coverIndLogitTransform=coverLogitTransform$coverInd, smCoverIndLogitTransform=coverLogitTransform$smCoverInd))
}

# 'getEstVE' performs 1 MC iteration, i.e., it generates the data-set and estimates the VE(s1) curve
getEstVE <- function(s1grid, n, beta, pi, sigma12, truncateMarker, seed){
  data <- getData(n=n, beta=beta, pi=pi, sigma12=sigma12, truncateMarker=truncateMarker, seed=seed)
  return(VEcurve(data=data, s1grid=s1grid))
}

# 'getPestVE' performs 1 MC iteration, i.e., it generates the data-set and estimates the VE(s1) curve
# parametric Gaussian density estimation is employed
getPestVE <- function(s1grid, n, beta, pi, truncateMarker, seed){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  return(pVEcurve(data=data, s1grid=s1grid))
}

# 'getPestRisk' performs 1 MC iteration, i.e., it generates the data-set and estimates the risk curve in each treatment group
# parametric Gaussian density estimation is employed
getPestRisk <- function(s1grid, n, beta, pi, truncateMarker, seed){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  return(pRiskCurve(data=data, s1grid=s1grid))
}

# 'getInferenceVE' generates data representing 1 MC iteration, computes the bootstrap SE based on 'nBoot'
# bootstrap iterations, and returns a list with the following components:
# coverInd                - a vector of 0s and 1s indicating whether the truth is covered by the bootstrap Wald-type CI
# smCoverInd              - a 0 or 1 indicating coverage of the entire VE(s1) curve on the 's1grid' by the simultaneous bootstrap Wald-type CI
# pSizeTestConstancy      - a 2-sided p-value of {H0: VE(s1)=VE for all s1} against {H1: non-H0} in a scenario that satisfies H0
# pPowerTestConstancy     - a 2-sided p-value of {H0: VE(s1)=VE for all s1} against {H1: non-H0} in the same scenario as that used for coverage
# pSizeTestTrialEquality  - a 2-sided p-value of {H0: VE_x(s1)=VE_y(s1) for all s1 in [s_l,s_u]} against {H1: non-H0} in a scenario that satisfies H0
# pPowerTestTrialEquality - a 2-sided p-value of {H0: VE_x(s1)=VE_y(s1) for all s1 in [s_l,s_u]} against {H1: non-H0} in a scenario that satisfies H1
getInferenceVE <- function(s1grid, trueVEcurve, n, beta, betaSizeT1, betaPowerT3, pi, truncateMarker, seed, nBoot){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT1 <- getData(n=n, beta=betaSizeT1, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT3 <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  dataPowerT3 <- getData(n=n, beta=betaPowerT3, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  return(inferenceVEcurve(data=data, dataSizeT1=dataSizeT1, dataSizeT3=dataSizeT3, dataPowerT3=dataPowerT3, s1grid=s1grid, trueVEcurve=trueVEcurve, nBoot=nBoot))
}

getInferenceVE2 <- function(s1grid, trueVEcurve, n, beta, betaSizeT1, betaPowerT3, pi, truncateMarker, seed, nBoot){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT1 <- getData(n=n, beta=betaSizeT1, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT3 <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  dataPowerT3 <- getData(n=n, beta=betaPowerT3, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  return(inferenceVEcurve2(data=data, dataSizeT1=dataSizeT1, dataSizeT3=dataSizeT3, dataPowerT3=dataPowerT3, s1grid=s1grid, trueVEcurve=trueVEcurve, nBoot=nBoot))
}

# 'getPinferenceVE' generates data representing 1 MC iteration, computes the bootstrap SE based on 'nBoot'
# bootstrap iterations, and returns a list with the following components:
# coverInd                - a vector of 0s and 1s indicating whether the truth is covered by the bootstrap Wald-type CI
# smCoverInd              - a 0 or 1 indicating coverage of the entire VE(s1) curve on the 's1grid' by the simultaneous bootstrap Wald-type CI
# pSizeTestConstancy      - a 2-sided p-value of {H0: VE(s1)=VE for all s1} against {H1: non-H0} in a scenario that satisfies H0
# pPowerTestConstancy     - a 2-sided p-value of {H0: VE(s1)=VE for all s1} against {H1: non-H0} in the same scenario as that used for coverage
# pSizeTestTrialEquality  - a 2-sided p-value of {H0: VE_x(s1)=VE_y(s1) for all s1 in [s_l,s_u]} against {H1: non-H0} in a scenario that satisfies H0
# pPowerTestTrialEquality - a 2-sided p-value of {H0: VE_x(s1)=VE_y(s1) for all s1 in [s_l,s_u]} against {H1: non-H0} in a scenario that satisfies H1;
# parametric Gaussian density estimation is employed
getPinferenceVE <- function(s1grid, trueVEcurve, n, beta, betaSizeT1, betaPowerT3, pi, truncateMarker, seed, nBoot){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT1 <- getData(n=n, beta=betaSizeT1, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT3 <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  dataPowerT3 <- getData(n=n, beta=betaPowerT3, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  return(pInferenceVEcurve(data=data, dataSizeT1=dataSizeT1, dataSizeT3=dataSizeT3, dataPowerT3=dataPowerT3, s1grid=s1grid, trueVEcurve=trueVEcurve, nBoot=nBoot))
}

getPinferenceVE2 <- function(s1grid, trueVEcurve, n, beta, betaSizeT1, betaPowerT3, pi, truncateMarker, seed, nBoot){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT1 <- getData(n=n, beta=betaSizeT1, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT3 <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  dataPowerT3 <- getData(n=n, beta=betaPowerT3, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  return(pInferenceVEcurve2(data=data, dataSizeT1=dataSizeT1, dataSizeT3=dataSizeT3, dataPowerT3=dataPowerT3, s1grid=s1grid, trueVEcurve=trueVEcurve, nBoot=nBoot))
}

# 'contrast' could be one of "additive", "logRR", "ve" (currently only "additive" is implemented)
# 'normTransform' is the normalizing transformation; currently only one of "identity" and "logit"
getPinferenceMCEPcurve2 <- function(s1grid, trueRiskCurves, n, beta, betaSizeT1, betaPowerT3, pi, truncateMarker, seed, nBoot, contrast="additive"){
  data <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT1 <- getData(n=n, beta=betaSizeT1, pi=pi, truncateMarker=truncateMarker, seed=seed)
  dataSizeT3 <- getData(n=n, beta=beta, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  dataPowerT3 <- getData(n=n, beta=betaPowerT3, pi=pi, truncateMarker=truncateMarker, seed=seed+100000)
  return(pInferenceMCEPcurve2(data=data, dataSizeT1=dataSizeT1, dataSizeT3=dataSizeT3, dataPowerT3=dataPowerT3, s1grid=s1grid, trueRiskCurves=trueRiskCurves, nBoot=nBoot, 
                              contrast=contrast))
}
