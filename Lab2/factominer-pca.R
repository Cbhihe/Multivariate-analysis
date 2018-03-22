
# ######################################
#      FACTOMINER
# ######################################

##  MIRI:     MVA
##  LAB #2:   PCA
##  Authors:  Cedric Bhihe <cedric.bhihe@gmail.com>
##            Santi Calvo <s.calvo93@gmail.com>  
##  Delivery: before 2018.03.19 - 23:55

##  Script name: lab2-factominer-pca_mva.R


rm(list=ls(all=TRUE))

#library(MASS)    #  exec `install.packages("MASS")` to use 'fractions'
#library(FactoMineR)
#library("mice")
#library("DMwR")
#library("VIM")    # exec `install.packages("VIM")` in R shell first
#library("dplyr")  # exec `install.packages("dplyr")` in R shell first
                   # split-apply-combine, specialized to df
#library(chemometrics)  # exec `install.packages("chemometrics")` in R shell first
                        # MCD outliers based on Mahalanobis distance
library(FactoMineR) # exec `install.packages("FactoMineR")` in R shell first
                    # PCA analysis
require(graphics)
# require(ggplot2)  # exec `install.packages("ggplot2")` in R shell first
# library(ggrepel)

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

X$demo <- as.factor(X$demo)
par(mfrow = c(1,3))
pcaX <- PCA(X,quali.sup=c(9))
#attributes(pcaX)
pcaX$quali.sup
pcaX$eig  # compare with evals and cumulative percentage of variance obtained manually. 


# SCREENING OF EIGENVALUE
plottitle = sprintf("Screening of eigenvalues (uniform obs. weights)")
plot(pcaX$eig[,1],
     pch=15, 
     cex=1,
     col="blue",
     type="b",
     main=plottitle,
     #sub="(Red labels show cumulative variance explanatory power)",
     xlab="Index (sorted)",
     ylab="Eigenvalues")  # screen eigenvalues
abline(h=mean(pcaX$eig[1:8]),col="gray")
text(x=6, y=mean(pcaX$eig[1:8])-0.004, 
     labels=paste0("Mean eigenvalue: ",as.character(format(round(mean(pcaX$eig[1:8]),3),digits=3))),
     cex=0.75,
     pos=3,
     col="red")
psiFMR <- pcaX$ind$coord # psi matrix

par(mfrow = c(1,1))

# BEST AND WORST COUNTRY REPRESENTATION IN PC-1-PC2

#In 1st PC plane, PC1 x PC2, compute represented inertia for each individual
represented12 <- c()
for (ii in 1:nrow(psiFMR)) {
  represented12 <- c(represented12,
                     round(100*norm(as.matrix(psiFMR[ii,1:2]),type="F")/norm(as.matrix(psiFMR[ii,]),type="F"),2))
  #cat (rownames(psiFMR)[ii]," ",represented12[ii],"% \n")
  }
cat("\n\nSorted individuals' projection's representativeness in PC1-PC2\n")
for (ii in 1:nrow(psiFMR)) {
  cat(rownames(psiFMR)[order(represented12, decreasing=T)][ii], 
      represented12[order(represented12, decreasing=T)][ii],"\n")
      }

# three countries most influencing the formation of the first PC
# Answer:  Cuba, Switzerland, Canada
for (ii in 1:3) {
  cat(rownames(psiFMR)[order(abs(psiFMR[,1]), decreasing=T)][ii], 
      round(psiFMR[order(abs(psiFMR[,1]), decreasing=T)][ii],2),"\n")
}
# Which are the three countries most influencing the formation of the second PC 
# Answer:  India, Yugoslavia, Cuba
psiFMRpc2 <- matrix(psiFMR[,2],ncol=1,byrow=T)
rownames(psiFMRpc2) <- rownames(psiFMR)
for (ii in 1:3) {
  cat(rownames(psiFMRpc2)[order(abs(psiFMRpc2[,1]),decreasing=T)][ii], 
      round(psiFMRpc2[order(abs(psiFMRpc2[,1]),decreasing=T)][ii],2),"\n")
}

# Which is the variable best represented in the first factorial plane?
# Which are the worst represented in the first factorial plane?

phiFMR <- pcaX$var$cor # variable projection in R^n space

represented12 <- c()
for (ii in 1:nrow(phiFMR)) {
  represented12 <- c(represented12,
                     round(100*norm(as.matrix(phiFMR[ii,1:2]), type="F")/
                             norm(as.matrix(phiFMR[ii,]),type="F"),2))
  #cat (rownames(phiFMR)[ii]," ",represented12[ii],"% \n")
}

for (ii in 1:nrow(phiFMR)) {
  cat(rownames(phiFMR)[order(represented12, decreasing=T)][ii], 
      represented12[order(represented12, decreasing=T)][ii],"\n")
}

# Which are the three variables most influencing the formation of the first principal component?
# Answer: farm, Gini, Gnpr
for (ii in 1:3) {
  cat(rownames(phiFMR)[order(abs(phiFMR[,1]), decreasing=T)][ii], 
      round(phiFMR[order(abs(phiFMR[,1]), decreasing=T)][ii],2),"\n")
}

# what are the three variables most influencing the formation of the second principal component?
# Answer: Gini, Instab, farm
phiFMRpc2 <- matrix(phiFMR[,2],ncol=1,byrow=T)
rownames(phiFMRpc2) <- rownames(phiFMR)
for (ii in 1:3) {
  cat(rownames(phiFMRpc2)[order(abs(phiFMRpc2[,1]),decreasing=T)][ii], 
      round(phiFMRpc2[order(abs(phiFMRpc2[,1]),decreasing=T)][ii],2),"\n")
}

# Which modalities of the variable “demo” are significant in the first 2 PCs?
(pcaX$quali.sup)  # look at suplpementary variables
plotellipses(pcaX,9)  # plot confidence ellipses around the categorical 9th variable "demo"
hcpcX <- HCPC(pcaX, nb.clust=0, conso=0, min=3, max=10)
