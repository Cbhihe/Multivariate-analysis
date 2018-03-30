##  MIRI:     MVA
##  LAB #3:   Beyond PCA
##  Authors:  Cedric Bhihe <cedric.bhihe@gmail.com>
##            Santi Calvo <s.calvo93@gmail.com>  
##  Delivery: before 2018.04.02 - 23:55

##  Script name: lab3-script_mva.R

rm(list=ls(all=TRUE))

library(MASS)         # to use 'fractions'
require(graphics)
library(FactoMineR)   # to use PCA method
library(chemometrics) # to invoke the canned algorithm NIPALS: function nipals()
                      # requires package "rparts" or "FactoMineR"
library(ggplot2)      # to enhance graph plotting
library(ggrepel)      # to plot with well behaved labeling


options(scipen=6) # R switches to sci notation above 5 digits on plot axes
ccolors=c("red","green","blue","orange","cyan","tan1","darkred","honeydew2","violetred",
          "palegreen3","peachpuff4","lavenderblush3","lightgray","lightsalmon","wheat2")

set.seed(932178)

#setwd("C:/Users/calvo/Desktop/UPC/Courses/Third_semester/MVA/Exercises/Labs/")
setwd("~/Documents/Study/UPC/MIRI-subjects/MVA_multivariate-analysis/Labs/")

#############################################
# 1: Load imputed Russet data from Lab 1
#############################################
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

Xdemo <- X$demo   # save categorical variable in object 'Xdemo' 
X_std <- scale(X[,-9],T,T) # get rid of categorical variable in main data-set
N <- diag(rep(1,nrow(X))/nrow(X),nrow(X_std),nrow(X_std))  # uniform weight matrix

#############################################
# 2: PC determination using NIPALS
#############################################

# NIPALS = Nonlinear Iterative Partial Least Squares
# Per question 4 of Lab #2, we retain the observation "Cuba" on the ground that 
# it is significant.

nd  <- 3  # consider 3 significant dimensions
lbd <- c()

psi <- c()   # scores     from X = psi %*% t(p)
ldg <- c()   # loadings as weights for each original variable when calculating PCs
#psi <- as.data.frame(matrix(nrow=nrow(X), ncol=nd)) 
#p <- as.data.frame(matrix(nrow=ncol(X), ncol=nd))

X_tmp <- X_std
threshold = 1e-5  # convergence criterion on eigenvalues

# calculate norm2 (Frobenius) of each column in matrix Y
maxresF <- function(Y) {
  maxres=c()
  for (cc in 1:ncol(Y)) { maxres <- cbind(maxres,norm(as.matrix(Y[,cc]),type="F")) }
  return(maxres)
}

for(dd in 1:ncol(X_std)) {
  maxres <- maxresF(X_tmp) 
  # psi_tmp <- rowmeans(X_tmp)                    # initial scores (iter=0)
  psi_tmp <- X_tmp[,which(maxres == max(maxres))] # initial scores (iter=0)
  lbd_old <- 1e4
  lbd_tmp <- 0
  
  while (abs(lbd_tmp -lbd_old) >= threshold) {
    lbd_old <- lbd_tmp
    ldg_tmp <- t(as.matrix(X_tmp)) %*% as.matrix(psi_tmp) # loadings
    ldg_tmp <- ldg_tmp / sqrt(sum(ldg_tmp * ldg_tmp))     # normalize loadings to 1
    # psi_tmp <- sqrt(N) %*% as.matrix(X_tmp) %*% as.matrix(ldg_tmp)    # calculate new scores
    psi_tmp <- as.matrix(X_tmp) %*% as.matrix(ldg_tmp)    # calculate new scores
    lbd_tmp <- t(psi_tmp) %*% psi_tmp / nrow(X)           # calculate new eigenvalue
  }
  lbd <- rbind(lbd,lbd_tmp)
  psi <- cbind(psi,psi_tmp)  # psi, projections of observations on eigenvectors in R^p
  ldg <- cbind(ldg,ldg_tmp)  # PCA loadings, XX' eigenvectors in individual obs space R^p,
  # ldg is "p x n", equivalent to evecs in Lab#2
  X_tmp <- X_tmp - psi_tmp %*% t(ldg_tmp)
}

colnames(psi) <- paste0("PC",1:ncol(psi)) # psi, projections of observations on eigenvectors in R^p
rownames(psi) <- rownames(X_std)

colnames(ldg) <- paste0("PC",1:ncol(ldg)) # loadings, XX' eigenvectors in individual obs space R^p
rownames(ldg) <- colnames(X_std)

# NIPALS cumulative variance explanatory power in principal directions (per eigenvalues)
cumvarexp <- c(lbd[1]/sum(lbd))
for (ii in 2:nd) {cumvarexp <- c(cumvarexp,cumvarexp[ii-1]+lbd[ii]/sum(lbd))}
cat("Eigenvalues derived from custom NIPALS\n     (3rd decimal round-off):\n",
    round(lbd,2),
    "\nCumulative variance representative power in principal directions (%):\n     (2nd decimal round-off):\n",
    round(100*cumvarexp,1),"\n")

# check that loading vectors verify the norm <=1 condition
#for (rr in 1:nrow(p)) {  cat (diag(sqrt(sum(ldg[rr,]**2))))  }
diag(sqrt(t(ldg[1:nrow(ldg),]) %*% ldg[1:nrow(ldg),] ))

# Compute variable space (R^n) eigenvectors
LBD <- matrix(0,nd,nd) ; diag(LBD) <- lbd[1:nd]
evecs_var <- sqrt(N) %*% psi %*% solve(LBD^(0.5))
diag(sqrt(t(ldg[,1:nd]) %*% ldg[,1:nd]))  # verify that every vector is normed

# ###############################
# Variables' projections on principal directions
#fooX <- sqrt(N) %*% as.matrix(X_std) %*% as.matrix(t(X_std)) %*% sqrt(N) # n x n matrix, with n=47 
# check that independently computed eigenvectors in variable space (R^n) are as evecs_var
#evecs_var2 <- eigen(fooX)$vectors

# ###############################
# check based on correlation matrix, corX ( or on cor(X_std) )
#corX <- t(X_std) %*% N  %*% X_std
#evecs <- eigen(corX)$vectors
#psi2 <- X_std %*% evecs

# ##############################
# check based on canned NIPALS (package: "chemometrics") algorithm
# Use Gram-Schmidt orthogonalization at each iter step to avoid flop error accumulation ?
canned_nipals <- nipals(X_std,nd,it=400,tol=threshold) #,gramschmidt=TRUE)
#canned_nipals$T  # PCA scores  
#canned_nipals$P  # PCA loadings, eigenvectors in individual obs space R^p, equivalent to evecs
colnames(canned_nipals$P) <- paste0("PC",1:nd)
# check that (based on three first components) loading vectors verify the norm <=1 condition
# diag(sqrt(t(canned_nipals$P[1:nrow(canned_nipals$P),]) %*% canned_nipals$P[1:nrow(canned_nipals$P),])) 
diag(sqrt(t(canned_nipals$P[1:nd,]) %*% canned_nipals$P[1:nd,])) 


#############################################
# 3: Biplot in R^p based on canned and custom NIPALS 
#    for X = N^(-1/2) U V'
#############################################

# Interpretation: variables only give direction of growth
par(mfrow=c(1,2))
biplot(canned_nipals$T,canned_nipals$P,
       var.axes=T,
       xlab="PC1",ylab="PC2",
       cex = rep(par("cex"), 0.75),
       arrow.len = 0.1,
       main=bquote("Biplot in R"^"p"*" (uniform obs weights)"), #\nnot excluding Cuba"),
       sub="(From chemometrics::NIPALS)") 
#biplot(psi*sqrt(nrow(X_std)),ldg /sqrt(nrow(X_std)),
biplot(psi[,1:nd],ldg[,1:nd],
       var.axes=T,
       cex = rep(par("cex"), 0.75),
       arrow.len = 0.1,
#     xlim=4,ylim=4,
       main=expression("Biplot in R"^"p"*" (uniform obs weights)"), #\nnot excluding Cuba"),
       sub="(From custom NIPALS)") 


#############################################
# 4: VARIMAX
#############################################

var_rot <- varimax(ldg[,1:nd] %*% sqrt(LBD),normalize=T,eps=1e-5)
phi_rot <- var_rot$loadings

# plot in rotated PC1-2 factorial plane
cat("plot variables' projection in PC1-2 factorial plane\n")
rotvec_labels <- colnames(X_std)
# compute unit radius circle
theta <- seq(-pi, pi, length = 200); circ_data <- data.frame(xc=cos(theta),yc=sin(theta))
# plotfile <- sprintf("Lab2/Report/%s_var-proj12_%s.pdf",
#                     datestamp,
#                     substr(wflag,1,4))
plottitle = sprintf("Rotated Variables\' projection in PC1-2 factorial plane\n(uniform observations' weights)")
plotdata <- data.frame(PC1=phi_rot[,1],PC2=phi_rot[,2],z=rotvec_labels)
varproj_plot <- ggplot(data = plotdata) + 
  theme_bw()+
  geom_vline(xintercept = 0, col="gray") +
  geom_hline(yintercept = 0, col="gray") +
  #geom_text_repel(aes(PC1,PC2,label = z))+
  geom_text_repel(aes(PC1,PC2,label = z),
                  size=4,
                  point.padding = 0.5,
                  box.padding = unit(0.55, "lines"),
                  segment.size = 0.3,
                  segment.color = 'grey') +
  geom_point(aes(PC1,PC2),col = "blue", size = 1) +
  #geom_path(circ_data,aes(xc,yc), inherit.aes =F) +
  #geom_point(aes(x_coord,y_coord),col = "black", size = 0.2) +
  geom_segment(data = plotdata, 
               mapping = aes(x = 0, y = 0, xend = PC1, yend = PC2),
               color="blue",
               arrow=arrow(length=unit(2,"mm")))+
  labs(title = plottitle)+
  coord_fixed()
varproj_plot + geom_path(aes(xc, yc), data = circ_data, col="grey70")


#############################################
# 5: Scores of individuals projected on rotated components
#      using FactoMineR's PCA
#############################################

pcaX <- PCA(X_std,quali.sup=c(9),graph=F)
pcaX_psi <- pcaX$ind$coord[,1:nd]  # projections of individuals on PC axes; scores
pcaX_phi <- pcaX$var$coord[,1:nd]  # correlations with principal directions
pcaX_var_rot <- varimax(pcaX_phi)
pcaX_phi_rot <- pcaX_var_rot$loadings[1:ncol(X_std),]
# evals in newly rotated PC nd-dimensional space
pcaX_lbd_rot = diag(t(pcaX_var_rot$loadings) %*% pcaX_var_rot$loadings)
# check that sum of SS loadings = sum of evals before rotation 
sum(pcaX_lbd_rot);sum(pcaX$eig[1:nd,1])

pcaX_psi_rot <- X_std %*% solve(cor(X_std)) %*% pcaX_phi_rot %*% diag(sqrt(pcaX_lbd_rot))

# place calculated pcaX_psi_rot components in pcaX$ind$coord[,1:nd]
pcaX$ind$coord[,1:nd] <- pcaX_psi_rot
# point out variables (and categories) most characteristic for each principal directions
dimdesc(pcaX,axes=1:nd,proba=0.5)

#############################################
# 6-7: Read PCA_quetaltecaen data and symmetrize the matrix
#############################################

# express joint feeling between CCAA.


#############################################
#  8: Transform similarity matrix into dissimilarity matrix 
#  9: Perform PCA upon dissimilarity matrix
# 10: Plot the first two components
#############################################

# notice that max. similarity = 10
