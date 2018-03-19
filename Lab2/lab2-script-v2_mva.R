##  MIRI:     MVA
##  LAB #2:   PCA
##  Authors:  Cedric Bhihe <cedric.bhihe@gmail.com>
##            Santi Calvo <s.calvo93@gmail.com>  
##  Delivery: before 2018.03.19 - 23:55


rm(list=ls(all=TRUE))
library(MASS)    #  exec `install.packages("MASS")` to use 'fractions'
#library(FactoMineR)
#library("mice")
#library("DMwR")
#library("VIM")    # exec `install.packages("VIM")` in R shell first
#library("dplyr")  # exec `install.packages("dplyr")` in R shell first
                  # split-apply-combine, specialized to df
#library(chemometrics)  # exec `install.packages("chemometrics")` in R shell first
                  # to have access to MCD outliers based on Mahalanobis distance
require(graphics)
require(ggplot2)  # exec `install.packages("ggplot2")` in R shell first
library(ggrepel)

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

pcaF <- function(X,datestamp,wflag,wparam,...) {
  #attach(X)
  Xdemo <- X$demo
  cat("X$demo:",Xdemo)
  X <- X[,-9]

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
  try(if(abs(sum(diag(N))-1) >= 1e-4) 
    stop("WARNING: invalid normalization of individual weights"))  # check that matrix trace = 1
  
  # centroid G of individuals.
  X_wgt <-  N %*% as.matrix(X)  # weighed observations
  centroid <- apply(X_wgt, 2, sum)  # use 'sum' because N is made of normalized ponderation factors
  rm(X_wgt)
  
  # centered X matrix
  X_ctd <- as.matrix(X - matrix(rep(t(centroid), nrow(X)), ncol=ncol(X), byrow=T))
  colnames(X_ctd) <- colnames(X)
  # covariance matrix (on) X_ctd) 
  covX <- t(X_ctd) %*% N  %*% X_ctd
  # compare with cov(X) or cov(X_ctd)  ###########################################################
  print("ok 'covX'")
  
  # standardized X matrix
  colsd <- c()
  #X_sd <- 1/sqrt(diag(covX))
  for (cc in 1:ncol(X)) { colsd <- c(colsd,1/sd(X[,cc])) } # build vector of variables' sd(X[,.])
  X_std <- X_ctd %*% diag(colsd,ncol(X),ncol(X)) # build standardized observation matrix
  # for (cc in 1:ncol(X)) { X_std[,cc] <- as.matrix(X_ctd[,cc]/sd(X[,cc])) } # build standardized observation matrix
  colnames(X_std) <- colnames(X)
  # correlation matrix (on X_std)
  corX <- t(X_std) %*% N  %*% X_std
  # compare with cor(X)
  print("ok 'corX'")
  
  evals <- eigen(corX)$values
  cat("Eigenvalues (obs): ",round(evals,4),"\n")
  evecs <- eigen(corX)$vectors
  cat("Eigenvectors (obs):","\n",evecs)
  cat("Rank of observation matrix: ",min(ncol(X),ceiling(sum(diag(corX)))),"\n\n")
  
  # screen eigenvalues -> ascertain nbr of significant dimensions
  #   -> cumulative variance explanatory power of principal directions 
  #     (evecs) associated with eigenvalues
  cum_exp_pow=matrix(0,length(evals)+1,1)
  cat("Rank ","eval_power ", "eval_cumul_power\n") 
  for (ii in (1:length(evals))) {
    cum_exp_pow[ii+1] <- cum_exp_pow[ii] + 100*evals[ii]/sum(evals)
    cat(ii," ",round(100*evals[ii]/sum(evals),2)," ",round(cum_exp_pow[ii+1],2),"\n")
    #cum_varexp <- c(cum_varexp,100*(cum_varexp+evals[ii]/sum(evals))) 
  }
  
  # save plot in pdf file
  plotfile = sprintf("Lab2/Report/%s_screen-evals_%s.pdf",
                     datestamp,
                     substr(wflag,1,4))
  pdf(file = plotfile)    # open pdf file
  plottitle = sprintf("PC explanatory power for %s observation weights",wflag)
  plot(seq(1:length(evals)),evals,
       pch=15, 
       cex=1,
       col="blue",
       type="b",
       main=plottitle,
       sub="(Red labels show cumulative variance explanatory power)",
       xlab="Index (sorted)",
       ylab="Eigenvalues"
  )
  text(x=1:length(evals), y=evals, 
       labels=as.character(round(cum_exp_pow[-1,1],2)),
       cex=0.75,
       pos=1,
       col="red")  # add labels
  abline(h=mean(evals),col="gray")
  text(x=6, y=mean(evals)-0.004, 
       labels=paste0("Mean eigenvalue: ",as.character(round(mean(evals),3))),
       cex=0.75,
       pos=3,
       col="red")
  #screeplot(prcomp(covX),npcs=min(ncol(X),length(evals)), type="l", log="y")
  grid()
  dev.off()    # close off pdf file socket
  
  
  
  # projections of centered individuals in the EV's direction, psi (n=47 by p=9 matrix)
  psi <- X_std %*% evecs
  colnames(psi) <- paste0("PC",1:ncol(psi))
  cat("Individuals' projections on principal directions: ok","\n")
  
  # 1st check that roundoff is contained (compare with eigenvalues, evals).
  #    remember: sum of eigenvalues = rank, p, of multivariate distribution
  #              = trace of cor matrix computed on standardized data in R^p
  try(if(sum(abs(diag(t(psi)%*%N%*%psi) - evals)) >= 1e-4) 
    stop("WARNING: invalid roundoff error in eigenvector and/or eigenvalue computations."))   
  
  # 2nd check that roundoff is contained
  # toto <- 0; for (ii in (1:ncol(X))) toto=toto+var(X[,ii])
  # (1-sum(evals)/toto)
  
  # 3rd check that roundoff is contained
  # cat ("Eigenvalues (eigen):\n",evals,"\n")
  # cat ("Eigenvalues (var(psi columns)): ")
  # for (ii in (1:ncol(X))) {cat(var(psi[,ii])," ")}
  
  # projections of centered individuals in the PC1-2 factorial plane
  demo_col=as.factor(Xdemo)
  levels(demo_col) <- c("Stable", "Unstable", "Dictatorship")
  
  # plotfile = sprintf("Lab2/Report/%s_indiv-proj12_%s.pdf",
  #                    datestamp,
  #                    substr(wflag,1,4))
  # pdf(file = plotfile)    # open pdf file
  # plot(psi[,1],psi[,2],
  #      pch=18,
  #      cex=1,
  #      col=demo_col,
  #      type="p",
  #      main=plottitle,
  #      sub=plotsubtitle,
  #      xlab="PC_1",
  #      ylab="PC_2")
  # text(x=psi[,1], y=psi[,2],
  #      labels=rownames(X),
  #      cex=0.75,
  #      pos=3,
  #      col="black")  # add labels
  # abline(h=0,v=0, col="gray")
  # grid()
  # dev.off()    # close pdf file
  
  # In 1st PC plane, PC1 x PC2, compute represented fraction of individuals. 
  represented12 <- c()
  label12 <- c()
  for (ii in 1:nrow(psi)) {
    represented12 <- c(represented12,
                       round(100*norm(as.matrix(psi[ii,1:2]), type="F")/norm(as.matrix(psi[ii,]),type="F"),2))
    label12 <- c(label12,
                 paste0(rownames(psi)[ii],
                        " (",
                        round(100*norm(as.matrix(psi[ii,1:2]), type="F")/norm(as.matrix(psi[ii,]),type="F"),0),
                        "%)")
    )
  }
  cat("Represented inertia of individuals in PC1-PC2 plane:\n",label12) 
  
  # plot obs projection in PC1-2 factorial plane
  cat("plot obs projections in PC1-2 factorial plane\n")
  plottitle=sprintf("Individuals\' projection in PC1-2 factorial plane (%s weights)",wflag)
  plotdata <- data.frame(PC1=psi[,1],PC2=psi[,2],z=label12)
  plotfile <- sprintf("Lab2/Report/%s_indiv-proj12_%s.pdf",
                      datestamp,
                      substr(wflag,1,4))
  #pdf(file = plotfile)
  ggplot(data = plotdata) + 
    theme_bw() +
    geom_vline(xintercept = 0, col="gray") +
    geom_hline(yintercept = 0, col="gray") +
    geom_text_repel(aes(PC1,PC2,label = z),
                    size=3,
                    point.padding = 0.5,
                    box.padding = unit(0.55, "lines"),
                    segment.size = 0.3,
                    segment.color = 'grey') +
    geom_point(aes(PC1,PC2,col = factor(demo_col)), size = 2) +
    scale_color_discrete(name = 'Political\nregime') +
    labs(title = plottitle)
  #dev.off()
  ggsave(plotfile)

  
  # In 2nd PC plane, PC2 x PC3, compute represented fraction of individuals. 
  represented23 <- c()
  label23 <- c()
  for (ii in 1:nrow(psi)) {
    represented23 <- c(represented23,
                       round(100*norm(as.matrix(psi[ii,c(2,3)]), type="F")/norm(as.matrix(psi[ii,]),type="F"),2))
    label23 <- c(label23,
                 paste0(rownames(psi)[ii],
                        " (",
                        round(100*norm(as.matrix(psi[ii,2:3]), type="F")/norm(as.matrix(psi[ii,]),type="F"),0),
                        "%)")
    )
  }
  cat("Represented inertia of individuals in PC2-PC3 plane:\n",label23) 
  
  # plot in 2nd PC plane, PC2 x PC3
  cat("plot obs projections in PC2-3 factorial plane\n")
  plottitle=sprintf("Individuals\' projection in PC2-3 factorial plane (%s weights)",wflag)
  plotdata <- data.frame(PC2=psi[,2],PC3=psi[,3],z=label23)
  plotfile <- sprintf("Lab2/Report/%s_indiv-proj23_%s.pdf",
                      datestamp,
                      substr(wflag,1,4))
  #pdf(file = plotfile)
  ggplot(data = plotdata) + 
    theme_bw() +
    geom_vline(xintercept = 0, col="gray") +
    geom_hline(yintercept = 0, col="gray") +
    geom_text_repel(aes(PC2,PC3,label = z),
                    size=3,
                    point.padding = 0.5,
                    box.padding = unit(0.55, "lines"),
                    segment.size = 0.3,
                    segment.color = 'grey') +
    geom_point(aes(PC2,PC3,col = factor(demo_col)), size = 2) +
    scale_color_discrete(name = 'Political\nregime') +
    labs(title = plottitle)
  ggsave(plotfile)
  #dev.off()

  
  # In 3rd PC plane, PC1 x PC3, compute represented fraction of individuals. 
  represented13 <- c()
  label13 <- c()
  for (ii in 1:nrow(psi)) {
    represented13 <- c(represented13,
                       round(100*norm(as.matrix(psi[ii,c(1,3)]), type="F")/norm(as.matrix(psi[ii,]),type="F"),2))
    label13 <- c(label13,
                 paste0(rownames(psi)[ii],
                        " (",
                        round(100*norm(as.matrix(psi[ii,c(1,3)]), type="F")/norm(as.matrix(psi[ii,]),type="F"),0),
                        "%)")
    )
  }
  cat("Represented inertia of individuals in PC1-PC3 plane:\n",label13) 
  
  # plot in 3rd PC plane, PC1 x PC3
  cat("plot obs projections in PC1-3 factorial plane\n")
  plottitle=sprintf("Individuals\' projection in PC1-3 factorial plane (%s weights)",wflag)
  plotdata <- data.frame(PC1=psi[,1],PC3=psi[,3],z=label13)
  plotfile <- sprintf("Lab2/Report/%s_indiv-proj13_%s.pdf",
                      datestamp,
                      substr(wflag,1,4))
  #pdf(file = plotfile)
  ggplot(data = plotdata) + 
    theme_bw() +
    geom_vline(xintercept = 0, col="gray") +
    geom_hline(yintercept = 0, col="gray") +
    geom_text_repel(aes(PC1,PC3,label = z),
                    size=3,
                    point.padding = 0.5,
                    box.padding = unit(0.55, "lines"),
                    segment.size = 0.3,
                    segment.color = 'grey') +
    geom_point(aes(PC1,PC3,col = factor(demo_col)), size = 2) +
    scale_color_discrete(name = 'Political\nregime') +
    labs(title = plottitle)
  ggsave(plotfile)
  #dev.off()

  
  # Write representativeness of individual projection in three factorial planes to disk
  representPC123 <- data.frame(cbind(Countries=rownames(X),
                                     Obs_weights=round(diag(N),3),
                                     PC12=represented12,
                                     PC23=represented23,
                                     PC13=represented13))
  colnames(representPC123) <- c("Countries", 
                                "Observation weights", 
                                "PC1-PC2 inertia (%)", 
                                "PC2-PC3 inertia (%)",
                                "PC1-PC3 inertia (%)")
  filename <- sprintf("Lab2/Report/%s_indiv-proj-inert_%s.csv",
                      datestamp,
                      substr(wflag,1,4))
  write.table(representPC123,file=filename,append=F,quote=F,
              sep=",",eol="\n",row.names=F,col.names=T)
  
  
  # projections of variables in the EV's direction, phi (n=47 by p=9 matrix)
  cat("Variables' projections on principal directions","\n")
  fooX <- sqrt(N) %*% as.matrix(X_std) %*% as.matrix(t(X_std)) %*% sqrt(N) # n x n matrix, with n=47 
  
  # compute eigenvalues in R^n
  evals_var <- eigen(fooX)$values
  # check that vectors are normed, using norm(as.matrix(evecs_var[,x]),type="F") for 1<=x<=8
  #norm(as.matrix(evecs_var[,5]),type="F")
  
  # display relative explanatory power of eigenvalues in R^n
  cum_exp_pow_evals_var=matrix(0,length(evals_var)+1,1)
  cat("Rank ","eval_var_power (%) ", "eval_var_cumul_power (%)\n") 
  for (ii in (1:length(evals_var))) {
    cum_exp_pow_evals_var[ii+1] <- cum_exp_pow_evals_var[ii] + 100*evals_var[ii]/sum(evals_var)
    #cum_varexp <- c(cum_varexp,100*(cum_varexp+evals[ii]/sum(evals))) 
    if (round(cum_exp_pow_evals_var[ii+1],1) < 100 ) {
      cat(ii," ",round(100*evals_var[ii]/sum(evals_var),2)," ",round(cum_exp_pow_evals_var[ii+1],2),"\n")
      jj=ii+1 # last index to display
      }
  }
  cat(jj," ",round(100*evals_var[jj]/sum(evals_var),2)," ",round(cum_exp_pow_evals_var[jj+1],2),"\n")
  
  # compute eigenvalues in R^n
  evecs_var <- eigen(fooX)$vectors
  cat("Rank of observation matrix: ",min(nrow(X),ceiling(sum(diag(fooX)))),"\n\n")
  
  phi <- as.matrix(t(X_std)) %*% evecs_var[,1:8]
  colnames(phi) <- rownames(X)
  # row names are variables names
  
  cat("plot variables' projection in PC1-2 factorial plane\n")
  plotfile <- sprintf("Lab2/Report/%s_var-proj12_%s.pdf",
                      datestamp,
                      substr(wflag,1,4))
  plottitle = sprintf("Variables\' projection in PC1-2 factorial plane (%s obs. weights)", wflag)
  
  # pdf(file = plotfile)
  # plot(phi[,1],phi[,2],
  #      pch=15, 
  #      cex=1,
  #      col="blue",
  #      type="p",
  #      main=plottitle,
  #      # sub=sprintf("(%s obs. weights)", wflag),
  #      xlab="PC_1_var",
  #      ylab="PC_2_var")
  # text(x=phi[,1], y=phi[,2], 
  #      labels=rownames(phi),
  #      cex=0.75,
  #      pos=1,
  #      col="red")  # add labels
  # grid()
  # dev.off()
  # plot in 3rd PC plane, PC1 x PC3

  plotdata <- data.frame(PC1=phi[,1],PC2=phi[,2],z=rownames(phi))
  ggplot(data = plotdata) + 
    theme_bw() +
    geom_vline(xintercept = 0, col="gray") +
    geom_hline(yintercept = 0, col="gray") +
    geom_text_repel(aes(PC1,PC2,label = z))+
    # geom_text_repel(aes(PC1,PC2,label = z),
    #                 size=3,
    #                 point.padding = 0.5,
    #                 box.padding = unit(0.55, "lines"),
    #                 segment.size = 0.3,
    #                 segment.color = 'grey') +
    geom_point(aes(PC1,PC2),col = "blue", size = 1) +
    geom_segment(data = plotdata, 
                 mapping = aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 color="green",
                 arrow=arrow(length=unit(4,"mm")),
                 alpha=0.80) +
    labs(title = plottitle)
  ggsave(plotfile)

  
} #  function closure
  


datestamp <- format(Sys.time(),"%Y%m%d-%H%M%S"); pcaF(X,datestamp,wflag="random")

datestamp <- format(Sys.time(),"%Y%m%d-%H%M%S"); pcaF(X,datestamp,wflag="uniform")

weights <- rep(1:10,5); length(weights) <- nrow(X) 
datestamp <- format(Sys.time(),"%Y%m%d-%H%M%S"); pcaF(X,datestamp,wflag="arbitrary",wparam=weights)

