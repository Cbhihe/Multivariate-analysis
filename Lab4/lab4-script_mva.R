##  MIRI:     MVA
##  LAB #4:   Clustering Analysis
##  Authors:  Cedric Bhihe <cedric.bhihe@gmail.com>
##            Santi Calvo <s.calvo93@gmail.com>  
##  Delivery: before 2018.04.23 - 23:55

##  Script name: lab4-script_mva.R

rm(list=ls(all=TRUE))

library(MASS)         # to use 'fractions'
require(graphics)
library(FactoMineR)   # to use PCA method
library("mice")
# library(chemometrics) # to invoke the canned algorithm NIPALS: function nipals()
# requires package "rparts" or "FactoMineR"
library(ggplot2)      # to enhance graph plotting
library(ggrepel)      # to plot with well behaved labeling


options(scipen=6) # R switches to sci notation above 5 digits on plot axes
ccolors=c("red","green","blue","orange","cyan","tan1","darkred","honeydew2","violetred",
          "palegreen3","peachpuff4","lavenderblush3","lightgray","lightsalmon","wheat2")

set.seed(932178)
datestamp <- format(Sys.time(),"%Y%m%d-%H%M%S"); 

#setwd("C:/Users/calvo/Desktop/UPC/Courses/Third_semester/MVA/Exercises/Labs/")
setwd("~/Documents/Study/UPC/MIRI-subjects/MVA_multivariate-analysis/Labs/")

#############################################
# 1: Calculate or load imputed Russet data from Lab 1
#############################################


# ----------------------------------------------
# Calculate imputation
list_files = as.vector(grep("^russet1964.*\\.csv$",
                            list.files(path="Data",full.names=FALSE),
                            ignore.case = TRUE,
                            perl = TRUE,
                            fixed = FALSE,
                            inv=FALSE,
                            value=TRUE)
)
russet_data = read.csv(paste0("Data/",list_files[1]),
                       header = TRUE, 
                       quote = "\"", 
                       dec = ".",
                       sep=" ", 
                       check.names=TRUE)
data_imp_mice <- mice(russet_data[2:10],
                      where=is.na(russet_data[2:10]),
                      method="pmm",
                      m = 8,
                      maxint=10);
X <- complete(data_imp_mice)
rownames(X) <- as.matrix(russet_data[1])

# ----------------------------------------------
# Load imputation
list_files = as.vector(grep("^russet.*\\mpp.csv$",
                            list.files(path="Data",full.names=FALSE),
                            ignore.case = TRUE,
                            perl = TRUE,
                            fixed = FALSE,
                            inv=FALSE,
                            value=TRUE)
)

X = read.csv(paste0("Data/",list_files[1]),
             header = TRUE, 
             quote = "\"", 
             dec = ".",
             sep=",", 
             check.names=TRUE)
# ----------------------------------------------

# save categorical variable in object 'Xdemo'
Xdemo <- X$demo
# get rid of categorical variable 'demo' in main data set
# X <- X[,-9]

Nobs=nrow(X) # number of observations / individuals

#############################################
# 2: Conduct PCA, relying on Factominer
#############################################

# -------- Perform PCA manually (to check FactoMineR results) ----------
# Uniform weight matrix
# N <- diag(rep(1,Nobs)/Nobs,Nobs,Nobs)  

# Load weights for Cuba-centered PCA (Cuba's weight is set to 0)
obs_weights <- read.csv(file=sprintf("Data/%s_arbitrary-wparam.csv","20180319-131250"))
W <- obs_weights[,2]
N <- diag(W/sum(W),Nobs,Nobs)  # build diagonal matrix with normalized (non-uniform) weights
rm(obs_weights)

# Initialize matrices (centered and standardized obs matrices):
X_ctd <- matrix(0,Nobs,ncol(X[,-9]))
X_std <- X_ctd

# Centroid G of individuals.
X_wgt <-  N %*% as.matrix(X[,-9])  # weighed observations
centroid <- apply(X_wgt, 2, sum)  # use 'sum' because N is made of normalized ponderation factors
rm(X_wgt)

# Centered X matrix
X_ctd <- as.matrix(X[,-9] - matrix(rep(t(centroid), Nobs), ncol=ncol(X[,-9]), byrow=T))
colnames(X_ctd) <- colnames(X[,-9])
# Covariance matrix, X_ctd 
covX <- t(X_ctd) %*% N  %*% X_ctd

# Standardized matrix, X_std (taking into acccount the non-uniform obs weight distribution)
X_std <- sweep(X_ctd,2,sqrt(diag(covX)),"/") # Is there a multiplicative correction of Nobs/(Nobs-1) on covX ?
colnames(X_std) <- colnames(X[,-9])
# Correlation matrix (on X_std)
corX <- t(X_std) %*% N  %*% X_std
# corX-cor(sqrt(N) %*% as.matrix(X_std))  is N or sqrt(N), obs ponderation matrix, justified here ?? check low residuals

# X_std2 <- scale(N %*% as.matrix(X[,-9]),T,T)
# corX2 <- t(X_std2) %*% N  %*% X_std2
# corX2-cor(N %*% as.matrix(X[,-9]))
# 
# X_std3 <- sweep(X_ctd,2,sqrt(diag(covX *Nobs/(Nobs-1))),"/")
# corX3 <- t(X_std3) %*% N  %*% X_std3
# corX3-cor(N %*% as.matrix(X[,-9]))

eigX <- eigen(corX) 
evalsMNL <- eigX$values
evecsMNL <- eigX$vectors


# -------- Perform PCA with FactoMineR ------------------
par(mfrow = c(1,3))
# pcaX <- PCA(X_std,scale.unit=T,row.w=W,quali.sup=c(9)) # crashes
pcaX <- PCA(X,scale=T, ind.sup = c(11), quali.sup = c(9))
pcaX$eig   # eigenvalues, the percentage of variance and the cumulative percentage of variance
evalsFMR <- pcaX$eig[,1]



#############################################
# 3: Number of significant dimensions
#############################################

# choose all evals >= 1 and calculate number of significant diemensions, nd
nd <- length(pcaX$eig[pcaX$eig[,1]>=1,1])

# check that the cumulative variance explanations is greater than 70%
try(if(pcaX$eig[nd,3] < 70) 
  stop("WARNING: insufficient variance representation with current significant dimensions.")) 

plottitle = sprintf("PC explanatory power\nfor uniform observation weights\n(with \"Cuba\" as outlier)")
plot(seq(1:length(evalsFMR)),evalsFMR,
     pch=15, 
     cex=1,
     col="blue",
     type="b",
     main=plottitle,
     sub="(Red labels show cumulative variance explanatory power)",
     xlab="Index (sorted)",
     ylab="Eigenvalues"
)
text(x=1:length(evalsFMR), y=evalsFMR, 
     labels=as.character(round(pcaX$eig[,3],2)),
     cex=0.75,
     pos=1,
     col="red")  # add labels
abline(h=mean(evalsFMR),col="gray")
text(x=6, y=mean(evalsFMR)-0.004, 
     labels=paste0("Mean eigenvalue: ",as.character(round(mean(evalsFMR),2))),
     cex=0.75,
     pos=3,
     col="red")
#screeplot(prcomp(covX),npcs=min(ncol(X),length(evals)), type="l", log="y")
grid()

par(mfrow = c(1,1))


#############################################
# 4: Hierarchical clustering
#    Clustering consolidation
#############################################

