##  MIRI:     MVA
##  LAB #2:   PCA
##  Authors:  Cedric Bhihe <cedric.bhihe@gmail.com>
##            Santi Calvo <s.calvo93@gmail.com>  
##  Delivery: before 2018.03.19 - 23:55


rm(list=ls(all=TRUE))
#library("mice")
#library("DMwR")
#library("VIM")    # exec `install.packages("VIM")` in R shell first
#library("dplyr")  # exec `install.packages("dplyr")` in R shell first
                  # split-apply-combine, specialized to df
#library(chemometrics)  # exec `install.packages("chemometrics")` in R shell first
                  # to have access to MCD outliers based on Mahalanobis distance
require(graphics)
#require(ggplot2)  # exec `install.packages("ggplot2")` in R shell first

options(scipen=6) # R switches to sci notation above 5 digits on plot axes
ccolors=c("red","green","blue","orange","cyan","tan1","darkred","honeydew2","violetred",
          "palegreen3","peachpuff4","lavenderblush3","lightgray","lightsalmon","wheat2")

set.seed(932178)

#setwd("C:/Users/calvo/Desktop/UPC/Courses/Third_semester/MVA/Exercises/Lab1v2/")
setwd("~/Documents/Study/UPC/MIRI-subjects/MVA_multivariate-analysis/Labs/")

#############################################
# Load imputed Russet data from Lab 1
#############################################
list_files = as.vector(grep("^russet.*\\mpp.csv$",
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
                       sep=",", 
                       check.names=TRUE)
dim(russet_data)  # list nbr of observations followed by nbr of vars
X <- russet_data

#############################################
# function for PCA analysis
#############################################

pcaF <- function(X,wflag,wparam, ...) {

  # "wflag" is in c("random","uniform","arbitrary")
  # if arg "wflag" is:
  #   <> "random"
  #      default weight distribution is random uniform 
  #      over interval [0,1], normalized to 1
  #   <> "uniform"
  #      arg "wparam" is ignored 
  #   <> "arbitrary"
  #      weight distribution is given by vector "wparam"
  #      e.g. c(1,2,3,4,5,6,123), normalized to 1.

  if (wflag == "random") {  
    # generate random weight for each individual
    W <- runif(nrow(X))
    # normalize weights to 1
    W <- W/sum(W)
    N <- diag(W,nrow(X),nrow(X))
    sum(diag(N))  # check that matrix trace = 1
    
  } else if (wflag == "uniform") {
    # generate uniform weights distribution for each individual
    N <- diag(1/nrow(X),nrow(X),nrow(X))
    sum(diag(N))  # check that matrix trace = 1
    
  } else if (wflag == "arbitrary") {
   try(if(length(wparam) != nrow(X) | !is.numeric(wparam)) 
     stop("WARNING: invalid individual weights given to function \"pcaF\".\n
          Program abort.\n"))
    N <- diag(wparam/sum(wparam),nrow(X),nrow(X))
    
  } else {
        stop("WARNING: invalid parameter \"wflag\" given to function \"pcaF\".\n
             Program abort.\n")
      }

  # centroid G of individuals.
  (centroid <- colMeans((X)))
  
  # centered X matrix.
  X_ctd <- X - matrix(rep(centroid, nrow(X)), ncol=9, byrow=T)
  rownames(X_ctd) <- rownames(X)
    # standardized X matrix
  for (cc in 1:ncol(X)) {
    X_std[,cc] <- X_ctd[,cc]/sd(X[,cc])  
    }
  rownames(X_std) <- rownames(X)
  
    # covariance matrix (on) X_ctd) 

  
    # correlation matrix (on X_std)
  
  
  
  
  return
} #  function closure

#X_cov <- cov(X)
X_cor <- cor(X)
X_cor[is.na(X_cor)] <- 0
p <- eigen(X_cor)$vectors
d <- diag(eigen(X_cor)$values)
p %*% d %*% solve(p)

mm <- data.matrix(X_cor)
mm == p

t(X) %% N %% X_std
