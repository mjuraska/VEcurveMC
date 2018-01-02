# this code computes MC point estimates of the VE(s1) curve and saves them in a file for evaluation of bias and MSE
# this estimation procedure considers parametric Gaussian density estimation

rm(list=ls(all=TRUE))

# Ted, please specify the path for VEcurveMC_myFunctions.R in your local repo
source("VEcurveMC_myFunctions.R")

# the output files are large, therefore I prefer that we store them in this shared folder rather than push them to the remote GitHub repo
outDir <- "h:/SCHARP/Sanofi/manuscript_VEcurveMethod/Routput"

# the total sample size
n <- 5000

# the indicator of whether marker distributions are truncated
truncateMarker <- TRUE

# grid of biomarker values on which the performance of the estimators is evaluated
if (truncateMarker){
  # calculate the values of the probit risk model coefficients
  beta2 <- -0.02
  # the marginal probability of infection in the placebo group set to 0.1
  beta0 <- getBeta0(beta2, meanS0=2, varS0=1, 0.1, scenario="A2", threshold=1.5)
  beta1start <- startBeta1(beta0, 0.75)
  beta3 <- getBeta3(beta0, beta1start, beta2, meanS1=3, varS1=1, 0.25)
  # the marginal probability of infection in the vaccine group set to 0.05
  beta1 <- getBeta1(beta0, beta2, beta3, meanS1=3, varS1=1, 0.05, scenario="A2", threshold=1.5)
  beta <- c(beta0, beta1, beta2, beta3)
  
  # the grid of marker values
  s1 <- seq(1.5, qnorm(0.95, mean=3, sd=1), by=0.05)
} else {
  # calculate the values of the probit risk model coefficients
  beta2 <- -0.02
  # the marginal probability of infection in the placebo group set to 0.1
  beta0 <- getBeta0(beta2, meanS0=2, varS0=1, 0.1)
  beta1start <- startBeta1(beta0, 0.75)
  beta3 <- getBeta3(beta0, beta1start, beta2, meanS1=3, varS1=1, 0.25)
  # the marginal probability of infection in the vaccine group set to 0.05
  beta1 <- getBeta1(beta0, beta2, beta3, meanS1=3, varS1=1, 0.05)
  beta <- c(beta0, beta1, beta2, beta3)
  
  # the grid of marker values
  q <- qnorm(c(0.05,0.95), mean=3, sd=1)
  s1 <- seq(q[1], q[2], by=0.05)
}

# the number of MC iterations
nMC <- 10500

# a sampling probability for sampling a Phase 2 subcohort with biomarker measurements
pi <- c(0.1, 0.25, 0.5)

for (i in 1:length(pi)){
  VE.MCestimates <- sapply(1:nMC, function(seed, s1grid, n, beta, pi, truncateMarker){ getPestVE(s1grid, n, beta, pi, truncateMarker, seed) }, s1grid=s1, n=n, beta=beta, pi=pi[i], truncateMarker=truncateMarker)
  write.table(VE.MCestimates, file=file.path(outDir, paste0("pVE_MCestimates_nMC=",nMC,"_N=",n,"_truncateMarker=",truncateMarker,"_pi=",pi[i],".txt")), quote=FALSE, row.names=FALSE, col.names=FALSE)
}
