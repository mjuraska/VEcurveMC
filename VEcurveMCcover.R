# for each MC iteration, this code computes indicators of covering the true VE(s1) values by the CIs and saves them in a file for evaluation of CI coverage

rm(list=ls(all=TRUE))

source("t:/vaccine/sanofipasteur/dengue/manuscript_VEcurveMethod/code/VEcurveMC_myFunctions.R")

outDir <- "t:/vaccine/sanofipasteur/dengue/manuscript_VEcurveMethod/Routput"

# calculate the values of the probit risk model coefficients
beta2 <- -0.02
# the marginal probability of infection in the placebo group set to 0.1
beta0 <- getBeta0(beta2, meanS0=2, varS0=1, 0.1, scenario="A2", threshold=1.5)
beta1start <- startBeta1(beta0, 0.75)
beta3 <- getBeta3(beta0, beta1start, beta2, meanS1=3, varS1=1, 0.25)
# the marginal probability of infection in the vaccine group set to 0.05
beta1 <- getBeta1(beta0, beta2, beta3, meanS1=3, varS1=1, 0.05, scenario="A2", threshold=1.5)
beta <- c(beta0, beta1, beta2, beta3)

# grid of biomarker values on which the performance of the estimators is evaluated
s1 <- seq(1.5, qnorm(0.95, mean=3, sd=1), by=0.05)

# VE(s1) curve evaluated on the grid
trueVEcurve <- trueVE(s1, beta0, beta1, beta2, beta3, 2, 1, 3, 1, 0.7, scenario="A2", threshold=1.5)

nMC <- 10500
nBoot <- 500
pi <- c(0.25, 0.5, 1)

for (i in 1:length(pi)){
  # the last row of 'coverVE' is for the simultaneous CI
  coverVE <- sapply(1:nMC, function(seed, s1grid, trueVEcurve, beta, pi, nBoot){ getCoverVE(s1grid, trueVEcurve, beta, pi, seed, nBoot) }, s1grid=s1, trueVEcurve=trueVEcurve, beta=beta, pi=pi[i], nBoot=nBoot)
  write.table(coverVE, file=file.path(outDir, paste0("VE_cover_nMC=",nMC,"_nBoot=",nBoot,"_pi=",pi[i],".txt")), quote=FALSE, row.names=FALSE, col.names=FALSE)
}
