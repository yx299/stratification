library(BART)
library(PReMiuMar)
library(mclust)
library(mcclust)
library(infotheo)
library(dplyr)
library(aricode)

# step 1: predict potential outcomes using BART 

# set the number of iterations in BART
nburnin <- 1000
npost <- 5000 

# data: training data
# data.test: testing data
# first column: treatment A
# second column: outcome Y
# third to last columns: covariates X
model <- wbart(data[,-2],data[,2],data.test,nskip=nburnin,
               ndpost=npost)

# predicted potential outcomes (posterior mean)
yhat.1 <- model$yhat.test.mean[1:n]
yhat.2 <- model$yhat.test.mean[n+(1:n)]
yhat.3 <- model$yhat.test.mean[2*n+(1:n)]

dat <- data.frame(cbind(cov.all,yhat.1,yhat.2,yhat.3))

# step 2: profile regression for clustering

# names of covariates
cov.names <- c("x1","x2","x3")

# run profile regression model
runInfoObj <- profRegr(yModel="MVN",xModel="Normal",nSweeps=1000,  
                       nBurn=1000,data=dat,covNames=cov.names,
                       outcome=c("yhat.1","yhat.2","yhat.3"))

# calculate the dissimilarity matrix
dissimObj <- calcDissimilarityMatrix(runInfoObj)

# calculate the optimal clustering
clusObj <- calcOptimalClustering(dissimObj)

# clustering allocation of the optimal partition
clust.nolim <- clusObj$clustering

# number of clusters in the optimal clustering
numclust <- clusObj$nClusters

# calculate adjusted rand index
# Z is the true clustering (ground truth)
ari <- adjustedRandIndex(clust.nolim,Z)

# calculate homogeneity and completeness
homo <- v.measure(Z,clust.nolim)[1]
complete <- v.measure(Z,clust.nolim)[2]
# calculate average risks and profiles
riskProfileObj <- calcAvgRiskAndProfile(clusObj)

# plot covariate and outcome profiles
clusterOrderObj <- plotRiskProfile(riskProfileObj,"profile.png")