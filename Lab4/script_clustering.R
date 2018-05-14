# ######################################
##  Topic:    Clustering Analysis
##  Authors:  Cedric Bhihe, Santi Calvo
##  Date:	  2018.04.23 - 23:55
##  Script name: script_clustering.R
# ######################################

rm(list=ls(all=TRUE))

library(MASS)           # to use 'fractions'
library(FactoMineR)     # to use PCA method
# library(chemometrics)   # invoke canned algorithm NIPALS: function nipals()
                        # requires package "rparts" or "FactoMineR"
library("mice")         # for MICE imputation
library("cluster")
require(graphics)       # enhanced graphics
library(ggplot2)        # to enhance graph plotting
library(ggrepel)        # to plot with well behaved labeling

 
options(scipen=6) # R switches to sci notation above 5 digits on plot axes
ccolors=c("red","green","blue","orange","cyan","tan1","darkred","honeydew2","violetred",
          "palegreen3","peachpuff4","lavenderblush3","lightgray","lightsalmon","wheat2")

set.seed(932178)
datestamp <- format(Sys.time(),"%Y%m%d-%H%M%S"); 

#setwd("C:/Users/calvo/Desktop/UPC/Courses/Third_semester/MVA/Exercises/Labs/")
setwd("~/Documents/Academic/UPC/MIRI/Subjects/MVA_multivariate-analysis/Labs/")

# ############################################
# 1: Calculate or load imputed Russet data from Lab 1
# ############################################


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
rm(russet_data)
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

# ############################################
# 2: Conduct PCA
# ############################################

# ##########
# Perform PCA manually (to check later FactoMineR results)
# ##########

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

# ##########
# Perform PCA with FactoMineR
# ##########

# Treat "demo" as supplementary (qualitative) variable
supFac=c(9)

# Treat "Cuba" as outlier, and non-active, supplementary  variable 
supInd=c(11)

# pcaX <- PCA(X_std,scale.unit=T,row.w=W,quali.sup=c(9)) # crashes as 'row.w=W' not implemented yet
par(mfrow=c(1,3),new=F)
pcaX <- PCA(X,scale=T, ind.sup=supInd, quali.sup=supFac)
pcaX$eig   # eigenvalues, the percentage of variance and the cumulative percentage of variance
evalsFMR <- pcaX$eig[,1]

Navar <- ncol(X)-length(c(supi,supf))  # number of active variables

# ############################################
# 3: Number of significant dimensions
# ############################################

# 1/ Choose all evals >= 1 and sort in decreasing order of intertia representativeness.
# 2/ Calculate cumulative intertia representation from each PC's representativeness.
# 3) If total explained inertia < 80% chose next eigenvalue in list of evals in decreasing order. 
# 4/ Extract 'nd' as number of significant dimensions

n_dim <- length(pcaX$eig[pcaX$eig[,1]>=1,1])
if(pcaX$eig[n_dim,3] < 80) {
  for (nd in (n_dim+1):Navar) {
    if(pcaX$eig[nd,3] >= 80) break
  }
}
# check that the cumulative variance explanations is greater than 80%
try(if(pcaX$eig[nd,3] < 80) 
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
     labels=paste0("Eigenvalue mean: ",as.character(round(mean(evalsFMR),2))),
     cex=0.75,
     pos=3,
     col="red")
#screeplot(prcomp(covX),npcs=min(ncol(X),length(evals)), type="l", log="y")
grid()


# ############################################
# 4: Probabilistic clustering with k-means replications
#    Hierarchical clustering, 
#    Clustering consolidation, k-means
# ############################################

# ##########
# Probabilistic k-means approach
# ##########

Nr=10 # nbr of replications
Nk=10 # max probed nbr (index) of clusters

in_over_tot <- matrix(NA,Nr,Nk) # ratio of within-cluster sum of squares over total sum of squares
CH_index<- matrix(NA,Nr,Nk)     # corrected ratio of between-cluster over tot within-cluster sum of squares 
# = Calinsky – Harabasz index

for (r in 1:Nr) {
  for (k in 2:Nk) {
    resC <- kmeans(Psi,k,iter.max = 10, nstart = 1)
    in_over_tot[r,k] <- resC$betweenss/resC$totss
    CH_index[r,k] <- (resC$betweenss/(k-1))/(resC$tot.withinss/(Nobs-k))
  }
}

# index has zero value for 1 cluster
in_over_tot[,1] <- 0  
CH_index[,1] <- 0

par(mfrow = c(1,2))
plot(colMeans(in_over_tot),
     ylab='Normalized within-cluster SS',
     type="b",col="blue")
plot(colMeans(CH_index),
     ylab='Calinsky-Harabasz index',
     type="b",
     col="dark red")

# OR ...
# ... plotting on the same graph, with 2 different vertical axes
par(mfrow = c(1,1), mar=c(5, 4, 4, 6) + 0.1)
# 1st plot
plot(colMeans(in_over_tot),
     pch=16, 
     axes=F,
     xlab='',ylab='',
     ylim=c(0,1),
     type="b",
     col="blue",
     main="Selection of optimal number of clusters (by k-means)")
axis(2,ylim=c(0,1),col='black',col.axis='black',las=1)
mtext('Normalized within-cluster SS',side=2,line=2.5)
# 2nd plot
par(new=T)
plot(colMeans(CH_index),
     pch=18, 
     axes=F,
     ylim=c(0,max(colMeans(CH_index))),
     xlab='',ylab='',
     type="b",
     col="dark red")
axis(4,ylim=c(0,max(colMeans(CH_index))),col='dark red',col.axis='dark red',las=1)
mtext('Calinsky-Harabasz index',col='dark red',side=4,line=2.5)
# abcissa
axis(1,pretty(range(1:Nk),10))
mtext("Index (cluster nbr)",side=1,col="black",line=2.5) 
# dashed line at optimal cluster index
abline(v=which(colMeans(CH_index)==max(colMeans(CH_index))),
       lty=2,
       col='gray')

# # legend
# legend("bottomright",
#        legend=c("Norm within-SS","CH index"),
#        text.col=c("blue","dark red"),
#        pch=c(16,18),
#        col=c("blue","dark red"))

# Use silhouette() method to confirm cluster nbr
par(mfrow = c(1,3),new=F)
for (kk in 2:4) {
  resC <- kmeans(Psi,kk,iter.max = 10, nstart = 1)
  sil <- silhouette(resC$cluster,dist(Psi))
  plot(sil,col=1:kk, main='')
}


# ##########
# Hierarchical clustering approach
# ##########

# Input for the clustering algorithm is the psi matrix even as the PCs form an artificial var space.
Psi <- pcaX$ind$coord[,1:nd]   # yields obs coordinates, for nd significant dim, in PC factorial space, R^min(p,nd)
Nobs <- Nobs-length(supInd)    # nbr of obs - nbr of supplementary individuals

distX <- dist(Psi, method = "euclidean")
treeX <- hclust(distX, method = "ward.D2")

par(mfrow=c(1,2),new=F)
plot(treeX,
     main='Hierarchical Clustering (Ward.D2)',
     xlab='Distance',
     cex=0.6)
abline(h=12,col='red')

barplot(treeX$height,
        xlab='Agglomeration index', 
        main='Clustering heights',
        ylab='Agglomeration criterion\'s value')
abline(h=12,col='red')
text(x=24, y=12-0.04, 
     labels="Dendrogram tree-cut",
     cex=0.75,
     pos=3,
     col="red")

Nclusters = 2
cutX <- cutree(treeX,k=Nclusters)
length(cutX)

# Centroids for 2 classes
centroids <- aggregate(Psi,list(cutX),mean)[,2:(nd+1)]

# Quality index
Bss <- sum(rowSums(centroids^2)*as.numeric(table(cutX)))
Tss <- sum(rowSums(Psi^2))
Ib <- 100*Bss/Tss

# Visualize partition
par(mfrow=c(1,1),new=F)

plot(Psi[,1],Psi[,2],
     xlab='PC1',ylab='PC2',
     pch=16,type="n",
     col=cutX,
     main="Clustering of observations in 2 classes",
     sub="(With \'Cuba\' as outlier in PCA)"
     )
text(Psi[,1],Psi[,2],col=cutX,labels=rownames(X[-11,]),cex = 0.8)
abline(h=0,v=0,col="gray")
# legend("topright",paste0("Class_",1:Nclusters),pch=20,col=c(1:Nclusters))

# OR
plotdata <- data.frame(PC1=Psi[,1],PC2= Psi[,2],z=rownames(X[-11,]))
plottitle <- "Clustering of observations in PC1-2 factorial plane\n(\'Cuba\' as outlier in PCA)"
ggplot(data = plotdata,col=cutX) +
  theme_bw()+
  # geom_vline(xintercept = 0, col="gray") +
  # geom_hline(yintercept = 0, col="gray") +
  geom_text_repel(aes(PC1,PC2,label = z),
                  col=cutX,
                  size=4,
                  point.padding = 0.5,
                  box.padding = unit(0.55, "lines"),
                  segment.size = 0.3,
                  segment.color = 'grey') +
  geom_point(aes(PC1,PC2),col = "blue", size = 1) +
  labs(title = plottitle)+
  coord_fixed()



# ##########
# Consolidation using k-means
# ##########

resC_consol <- kmeans(Psi,centers=centroids)

# Quality index
Bss <- sum(rowSums(resC_consol$centers^2)*resC_consol$size)  # resC_consol$betweenss
Wss <- sum(resC_consol$withinss)                             # resC_consol$tot.withinss
Ib_consol <- 100*Bss/(Bss+Wss)

# Plot
par(mfrow=c(1,1),new=F)
plot(Psi[,1],Psi[,2],
     xlab='PC1',ylab='PC2',
     pch=cutX+16,
     type="p",
     col=cutX+4,
     main="Consolidated clustering of observations in 2 classes",
     sub="(With \'Cuba\' as outlier in PCA)"
     )
points(centroids,pch=18,type='p',col='blue',cex = 1.5)
text(centroids,labels=paste0("G",unique(cutX)),pos=1,cex=1,2)
abline(h=0,v=0,col="gray")

# Use silhouette to confirm clustering result
sil <- silhouette(resC_consol$cluster,dist(Psi))
plot(sil,col=1:length(unique(cutX)), 
     main='Silhouette widths for consolidated clustering\n(\'Cuba\' as outlier in PCA analysis)')


# ############################################
# 5: Using catdes(), interpret cluster and represent them.
# ############################################

# 
cat_descript <- catdes(cbind(as.factor(resC_consol$cluster),X[-11,-9]),
       1,
       proba=0.05,
       row.w=NULL)

# ############################################
# 6: Assign Cuba to one of the defined clusters
# ############################################

# 
pcaX$ind.sup$coord[1:2]
