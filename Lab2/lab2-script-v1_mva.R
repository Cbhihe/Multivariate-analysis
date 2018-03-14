##  MIRI:     MVA
##  LAB #2:   PCA
##  Authors:  Cedric Bhihe <cedric.bhihe@gmail.com>
##            Santi Calvo <s.calvo93@gmail.com>  
##  Delivery: before 2018.03.19 - 23:55


rm(list=ls(all=TRUE))
library(MASS)    #  exec `install.packages("MASS")` to use 'fractions'
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

russet_imp_data = read.csv(paste0("Data/",list_files[1]),
                       header = TRUE, 
                       quote = "\"", 
                       dec = ".",
                       sep=",", 
                       check.names=TRUE)
dim(russet_imp_data)  # list nbr of observations followed by nbr of vars

# initialize matrices
X <- russet_imp_data ; rm(russet_imp_data)


#############################################
# function for PCA analysis
#############################################

pcaF <- function(X,wflag,wparam,...) {
  X_ctd <- matrix(0,nrow(X),ncol(X))
  X_std <- X_ctd

  # "wflag" is in c("random","uniform","arbitrary")
  # if arg "wflag" is:
  #   <> "random"
  #      default weight distribution is random uniform 
  #      over interval [0,1], normalized to 1
  #      arg "wparam" is ignored 
  #   <> "uniform"
  #      arg "wparam" is ignored 
  #   <> "arbitrary"
  #      weight distribution is given by vector "wparam"
  #      e.g. c(1,2,3,4,5,6,123), normalized to 1.

  if (wflag == "random") { 
    print("ok 'random'")
    # generate random weight for each individual
    W <- runif(nrow(X))
     
  } else if (wflag == "uniform") {
    print("ok 'uniform'")
    # generate uniform weights distribution for each individual
    W=rep(1,nrow(X))

  } else if (wflag == "arbitrary") {
    try(if(length(wparam) != nrow(X) | !is.numeric(wparam)) 
        stop("WARNING: invalid individual weights given to function 'pcaF'. 
        Program abort."))
    print("ok 'arbitrary'")
    W <- wparam
    
  } else {
    stop("WARNING: invalid parameter 'wflag' given to function 'pcaF'.
         Program abort.")
  }
  
  N <- diag(W/sum(W),nrow(X),nrow(X))  # build diagonal matrix with normalized weights
  rm(W)
  try(if(abs(sum(diag(N))-1) >= 10**-6) 
    stop("WARNING: invalid normalization of individual weights"))  # check that matrix trace = 1
  
  # centroid G of individuals.
  centroid <- colMeans((X))
  
  # centered X matrix.
  X_ctd <- as.matrix(X - matrix(rep(t(centroid), nrow(X)), ncol=9, byrow=T))
  rownames(X_ctd) <- rownames(X)
  colnames(X_ctd) <- colnames(X)
  # covariance matrix (on) X_ctd) 
  covX <- t(X_ctd) %*% N  %*% X_ctd
  # compare with cov(X)
  print("ok 'covX'")
  
  # standardized X matrix
  for (cc in 1:ncol(X)) {
    X_std[,cc] <- as.matrix(X_ctd[,cc]/sd(X_ctd[,cc]))
  }
  rownames(X_std) <- rownames(X)
  colnames(X_std) <- colnames(X)
  # correlation matrix (on X_std)
  corX <- t(X_std) %*% N  %*% X_std
  # compare with cor(X)
  print("ok 'corX'")

  evals <- round(eigen(covX)$values,3)
  cat("Eigenvalues: ",evals,"\n")
  cat("Eigenvectors:","\n"); (evecs <- eigen(covX)$vectors)
  cat("Rank of observation matrix: ",min(ncol(X),ceiling(sum(diag(corX)))),"\n")
  
  # screen eigenvalues to ascertain the nbr of significant dimensions
  #   -> cumulative variance explanatory power of principal directions 
  #     (evecs) associated with eigenvalues
  cum_varexp=matrix(0,length(evals)+1,1)
  for (ii in (1:length(evals))) {
    cum_varexp[ii+1] <- round(cum_varexp[ii] + 100*evals[ii]/sum(evals),3)
    cat(ii," ",evals[ii]/sum(evals)," ",cum_varexp[ii+1],"\n")
    #cum_varexp <- c(cum_varexp,100*(cum_varexp+evals[ii]/sum(evals))) 
  }
  
  # save plot in pdf file
  barplotfile = sprintf("Lab2/Report/screeplot_%s.pdf",wflag)
  pdf(file = barplotfile)    # open pdf file
  plottitle = sprintf("Screeplot for %s observation weights",wflag)
  plot(seq(1:length(evals)),evals, 
       pch=15, 
       cex=1,
       col="blue",
       type="b",
       main=plottitle,
       sub="(labels show cumulative variance explanatory power)",
       xlab="Index",
       ylab="Eigenvalues",
       log="y")
  text(x=1:length(evals), y=evals, 
       labels=as.character(round(cum_varexp[-1,1],2)),
       cex=0.75,
       pos=1,
       col="red")  # add labels
  #screeplot(prcomp(covX),npcs=min(ncol(X),length(evals)), type="l", log="y")
  grid()
  dev.off()    # close pdf file

    return()
} #  function closure


pcaF(X,wflag="random")

pcaF(X,wflag="uniform")

weights <- rep(1:10,5); length(wparam) <- nrow(X) 
pcaF(X,wflag="arbitrary",wparam=weights)

