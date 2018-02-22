# this code computes bootstrap estimates of the "risk curve" in each treatment arm in the pooled analysis of 9-16 year-olds in CYD14+15
# these bootstrap estimates are used to construct CIs and hypothesis tests reported in the application section of the treatment effect modification manuscript
# bootstrap estimates of the risk curve in each treatment arm are computed to allow inference about two mCEP curves: log RR and additive risk difference
# nonparametric kernel density estimation is employed

rm(list=ls(all=TRUE))

# Ted, please specify the path for CYD14and15_myFunctions.R in your local repo
source("CYD14and15_myFunctions.R")

# the output files are large, therefore I prefer that we store them in this shared folder rather than push them to the remote GitHub repo
outDir <- "t:/vaccine/sanofipasteur/dengue/manuscript_VEcurveMethod/Routput"

# the data are shared here so that the many people working on them can access them from the same location
# however, feel free to create a local copy elsewhere if it's easier for the process
dataDir <- "t:/vaccine/sanofipasteur/dengue/CYD14and15/Mon13CoR/data"

# target set: ITT set at-risk at month 13 with no prior infection
data <- read.csv(file.path(dataDir,"cyd14and15m13CoRdata.csv"), header=TRUE)
data <- subset(data, AGEYRS>=9)

bRiskCurve(data, markerName="AUC", iter=1000, saveFile="bRiskCurves_CYD14and15_9to16_AUCMB_hingePoint.RData", saveDir=outDir)
