##  MIRI:     MVA
##  LAB #5:   Correspondence Analysis
##  Authors:  Cedric Bhihe <cedric.bhihe@gmail.com>
##            Santi Calvo <s.calvo93@gmail.com>  
##  Delivery: before 2018.05.04 - 23:55

##  Script name: lab5-script_mva.R

rm(list=ls(all=TRUE))

library(MASS)           # to use 'fractions'
library(FactoMineR)     # to use PCA method

require(graphics)       # enhanced graphics
library(ggplot2)        # to enhance graph plotting
library(ggrepel)        # to plot with well behaved labeling
#library("factoextra")   # ggplot2-based elegant visualization
 
options(scipen=6) # R switches to sci notation above 5 digits on plot axes
ccolors=c("red","green","blue","orange","cyan","tan1","darkred","honeydew2","violetred",
          "palegreen3","peachpuff4","lavenderblush3","lightgray","lightsalmon","wheat2")

set.seed(932178)
datestamp <- format(Sys.time(),"%Y%m%d-%H%M%S"); 

#setwd("C:/Users/calvo/Desktop/UPC/Courses/Third_semester/MVA/Exercises/Labs/")
setwd("~/Documents/Academic/UPC/MIRI/Subjects/MVA_multivariate-analysis/Labs/")

# ############################################
# 1: Read PCA_quetaltecaen data
#############################################

rm(list=ls(all=TRUE))

qttc = read.csv("Data/PCA_quetaltecaen.csv",
                      header = TRUE, 
                      quote = "\"", 
                      dec = ".",
                      sep="\t", 
                      encoding="UTF-8",
                      check.names=TRUE)

rownames(qttc) <- qttc[,1]
qttc <- qttc[,-1]
dim(qttc)

## Exploring the contigency table

# 1/ Pearson chi-square test for significant association (dependence) between row & col categs
# H0: joint cell counts distr. in 2D contingency table = (row marginals) x (col marginals)
# DFs = 49
chisq.test(qttc)
# Conclusion: We cannot reject the null hypothesis at the rist 5% of erring. 
# There is independence of row and column variables.
# Warning message in chisq.test(qttc) : Chi-squared approximation may be incorrect !!!

# 2/ Balloon plot
library ("gplots")
par(mfrow=c(1,1))
dt <- as.table(as.matrix(qttc))
balloonplot(t(dt), main ="\'¿Qué tal te caen?\' data set",
            xlab ="", ylab="",
            rowsrt=0, colsrt=45,
            label = FALSE, 
            show.margins = T)

# 3/ Mosaic plot
mosaicplot(dt, 
           shade=T,
           las=2,
           main = "\'¿Qué tal te caen?\' data set")
# Both representations show that there appear to be no difference between random data
# and the distribution of counts in the contigency table at hand !
par(mfrow = c(1,1))


#############################################
# 2: Perform CA on data. 
#    How many significant dimensions are there ?
#    Interpret projections on 1st factorial plane.
#############################################

# total nbr of counts
N <- sum(qttc)
# frequency matrix
fij <- qttc/N
# row weights in terms of perceived population similarity for each autonomous region
( fi <- rowSums(fij) )  
# col weights in terms of average similarity of a given population as perceived across autonomous regions
( fj <- colSums(fij) )  # see that catalanes have the lowest perceived similarity across all spanish regions
sum(fi) ; sum(fj) # check = 1

# row profile cloud with incorporated metric effect 
# -> deformed cloud of points
Di <- diag(fi)
Dj <- diag(fj)
# First, define the cloud of rows, with the matrix of conditional frequencies of rows:
( CFj_given_i <- solve(Di) %*% as.matrix(fij) )
# apply(CFj_given_i,1,sum)     # check that each row sums to 1

Fim <- solve(Di) %*% as.matrix(fij) %*% solve(sqrt(Dj))
# center of gravity of rows
( rowCentroid <- sqrt(fj) )
# ( apply(Fim,2,mean) )       # compare: close enough !
# 'Fim' is a deformed cloud where distances can now be calculated in the usual Euclidian
# sense owing to the metric correction introduced by 1/(f_.j) for rare occurence of "j" cat 
# vars along rows "i".

# We can now apply cloud centering and define a new centered data set, X_ctd
oneM <- matrix(1,nrow=dim(Fim), ncol=dim(Fim), byrow=T)

X_ctd <- Fim - oneM %*% sqrt(Dj)

rownames(X_ctd) <- rownames(qttc)
colnames(X_ctd) <- colnames(qttc)

# ... for which we can calculate eigenvalues and eigenvectors
Z <- sqrt(Di) %*% X_ctd
( eigX <- eigen(t(Z) %*% Z) )
evecs <- eigX$vectors        # eigenvectors
evals <- eigX$values         # eigenvalues
# As expected the last eval is 0. Along with the last column in evecs, it correspond to
# the centroid 'sqrt(fj)'

cumvarexp <- c(100*evals[1]/sum(evals))
for (ii in 2:(ncol(X_ctd))) {
  if(ii==2) cat("Eigenvalues","\t","Cumulative variance (%)\n")
  cumvarexp <- c(cumvarexp,cumvarexp[ii-1] + 100*evals[ii]/sum(evals))
  cat(format(evals[ii-1],scientific=T),"\t",round(cumvarexp[ii-1],6),"\n")
}

# Number of significant dimensions
# 1/ Point successively to every component of the list of evals sorted by decreasing value 
# 2/ When total explained inertia exceeds 80% chose eval index as nbr of significant dims
n_dim <- ncol(X_ctd)-1
for (nd in 1:n_dim) {  if( cumvarexp[nd] >= 80 ) break  }
cat("Number of significant dimensions:",nd)

# ############################
# ... or perform a "canned" CA, using FactoMineR
par(mfrow = c(1,2))
caX <- CA(qttc,ncp=ncol(qttc)-1,
          graph=T,
          axes = c(1,2),
          row.w = fi, excl=NULL)

# Above will plot **factor** maps in 1st factorial plane PC1-2, where row and col factors are printed 
# together with distinct colors for easier differentiation. 
# Factors, psi_alpha, are projections of either centered row "i" profiles, X_ctd, or of centered
# col "j" profiles on PC direction u_alpha (evec)
# psi <- t(X_ctd) %*% evecs   # for row profiles

caX$eig                   # no difference with previous result
evals <- caX$eig[,1]
# screeplot
plottitle = sprintf("PC explanatory power\n(\"¿Qué tal te caen?\" data set)")
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
     labels=paste0(as.character(round(caX$eig[,3],1)),"%"),
     cex=0.75,
     pos=1,
     col="red") # add labels
abline(h=mean(evals),col="gray")
text(x=4, y=mean(evals)-0.0001, 
     labels=paste0("Eigenvalue mean: ",as.character(format(mean(evals),scientific=T))),
     cex=0.75,
     pos=3,
     col="red")
grid()

# Conclusion: blue colored points, which are close together, have similar row profile 
# according to the Chi square metric.  Reiterate interpretation of Lab#3.
# Distances between same-color points are distances in the Chi square sense.
# A red point (col) is a barycenter of the blue points (row) weighted by the  col. 
# profile of that red points. And vice versa. 
# A red point may appear close to a blue point but that is only qualitative and no 
# conclusion can be drawn from that fact.
plot(caX,invisible="col", label="row")
plot(caX,invisible="row", label="col")

par(mfrow = c(1,1))

# ############################
# ... or perform a full PCA
par(mfrow = c(1,2))
pcaX <- PCA(X_ctd,
            scale.unit=FALSE,
            ncp = ncol(qttc)-1, 
            ind.sup = NULL, quanti.sup = NULL, quali.sup = NULL, 
            row.w = fi, col.w = NULL, 
            graph = TRUE, axes = c(1,2))
pcaX$eig

# Number of significant dimensions
# 1/ Point successively to every component of the list of evals sorted by decreasing value 
# 2/ When total explained inertia exceeds 80% chose eval index as nbr of significant dims
n_dim <- ncol(X_ctd)-1
for (nd in 1:n_dim) {  if(pcaX$eig[nd,3] >= 80) break  }
cat("Number of significant dimensions:",nd)

# Conclusion: Choose to retain nd=3 first evals = nbr of significant dims. 
# Together they account for about 83% of overall inertia/variance.
par(mfrow = c(1,1))



#############################################
# 3: Compute contribution of each cell to total inertia 
#    Compute contribution to inertia of diagonal cells.
#############################################

# tot intertia contrib ('tic') of each cell in contigency table
tic <- matrix(0,nrow=nrow(fij),ncol=ncol(fij),byrow=T)
colnames(tic)=colnames(fij); rownames(tic)=rownames(fij)

for (ii in 1:nrow(fij)) {
  for (jj in 1:ncol(fij)) {
    tic[ii,jj] <- ( fij[ii,jj] - fi[ii] * fj[jj] )^2 / (fi[ii]*fj[jj])
  }
}

View(format(tic,scientific=T))

# Compute total inertia
totI <- sum (tic)  # chi-square statistics / N  (slide 40 in MVA slides06)
# sum(apply(tic,1,sum))   # check
cat("\nTotal inertia:",totI,"\n\n")

# Discussion:
#   Earlier we saw that the Chi-square statistics led us to the conclusion that col
#   and rows were independent. This would be consistent with the cloud of row profiles 
#   being concentrated around the centroid and the inertia to be zero or very close, as:
#                fij[11,jj] = fi[ii] * fj[jj]
# Conclusions:
cat ("\nDiagonal's contribution to inertia:",
     round(100*sum(diag(tic))/totI,2),
     "%","\nConclusion: Overloaded diagonal and Gutmann effect.\n")

caX$row$contrib
# apply(caX$row$contrib,2,sum)   # checks to 100% col-wise
# Conclusions: 
#   Andalusia contributes almost 40% to inertia in the 5th direction and 14.6% in the first
#   Catalonia contributes almost 38% and 30% to inertia in the 1st and 2nd direction
#   ...
caX$row$cos2
# apply(caX$row$cos2,1,sum)    # checks to 1 row-wise
# Quality of representation of row profiles in the given directions (correlations with PCs)



#############################################
# 4: Nullify influence of overloaded diagonal on total inertia 
#    After doing so, compute new contributions to inertia of cells.
#############################################

eps <- 1e-4         # convergence_criterion on total inertia
qttc_new <- qttc    # new count table
N_new <- N          # temp N, used for compute
tic_new <- tic      # temp tic, idem
totI_old <- 1e4
err <- 1e4
cnt <- 0

while (err > eps) {
  cnt <- cnt +1
  fij_new <- qttc_new/N_new
  fi_new <- rowSums(fij_new)
  fj_new <- colSums(fij_new)
  diag(qttc_new) <- N_new * fi_new * fj_new
  N_new <- sum(qttc_new)
  # total inertia
  for (ii in 1:nrow(fij_new)) {
    for (jj in 1:ncol(fij_new)) {
      tic_new[ii,jj] <- ( fij_new[ii,jj] - fi_new[ii] * fj_new[jj] )^2 / (fi_new[ii]*fj_new[jj])
    }
  }
  totI_new <- sum(tic_new)
  err <- abs(totI_new - totI_old)
  totI_old <- totI_new
}

# New contribution to inertia of cells
View(format(tic_tmp,scientific=T))
cat("\nTotal inertia:",totI_tmp,"\n\n")
cat ("\nDiagonal's new contribution to inertia:",
     round(100*sum(diag(tic_tmp))/totI_tmp,2),
     "%","\nConclusion: We got rid of the diagonal overloading and (per Q-5) of the Gutmann effect.\n")



#############################################
# 5: Perform new CA on modfied table
#    Interpret results
#############################################

par(mfrow = c(1,2))
caX_new <- CA(qttc_new,ncp=ncol(qttc_new)-1,
              graph=T,
              axes = c(1,2),
              row.w = fi_new, excl=NULL)

# As before, the above will plot **factor** maps in 1st factorial plane PC1-2, where row 
# and col factors are printed together with distinct colors for easier differentiation. 
# Factors, psi_alpha, are projections of either centered row "i" profiles, X_ctd, or of centered
# col "j" profiles on PC direction u_alpha (evec)
# psi <- t(X_ctd) %*% evecs   # for row profiles

caX_new$eig                   # no difference with previous result
evals_new <- caX_new$eig[,1]

# screeplot
plottitle = sprintf("PC explanatory power\n\"¿Qué tal te caen?\" data set\n(corrected for Gutmann effect)")
plot(seq(1:length(evals_new)),evals_new,
     pch=15, 
     cex=1,
     col="blue",
     type="b",
     main=plottitle,
     sub="(Red labels show cumulative variance explanatory power)",
     xlab="Index (sorted)",
     ylab="Eigenvalues"
)
text(x=1:length(evals), y=evals_new, 
     labels=paste0(as.character(round(caX_new$eig[,3],1)),"%"),
     cex=0.75,
     pos=1,
     col="red") # add labels
abline(h=mean(evals_new),col="gray")
text(x=4, y=mean(evals_new), 
     labels=paste0("Eigenvalue mean: ",as.character(format(mean(evals),scientific=T))),
     cex=0.75,
     pos=3,
     col="red")
grid()

# Number of significant dimensions
# 1/ Point successively to every component of the list of evals sorted by decreasing value 
# 2/ When total explained inertia exceeds 80% chose eval index as nbr of significant dims
n_dim <- ncol(X_ctd)-1
for (nd_new in 1:n_dim) {  if( caX_new$eig[nd_new,3] >= 80 ) break  }
cat("Number of significant dimensions:",nd_new)

plot(caX_new,invisible="col", label="row", sub="(Corrected Gutmann effect)")
plot(caX_new,invisible="row", label="col", sub="(Corrected Gutmann effect)")
par(mfrow = c(1,1))




# ... and check with a new full PCA
Di_new <- diag(fi_new)
Dj_new <- diag(fj_new)
# First, define the cloud of rows, with the matrix of conditional frequencies of rows:
( CFj_given_i_new <- solve(Di_new) %*% as.matrix(fij_new) )
# apply(CFj_given_i,1,sum)     # check that each row sums to 1

Fim_new <- solve(Di_new) %*% as.matrix(fij_new) %*% solve(sqrt(Dj_new))
# center of gravity of rows
rowCentroid_new <- sqrt(fj_new)
# ( apply(Fim,2,mean) )       # compare: close enough !
# 'Fim' is a deformed cloud where distances can now be calculated in the usual Euclidian
# sense owing to the metric correction introduced by 1/(f_.j) for rare occurence of "j" cat 
# vars along rows "i".

# Apply cloud centering and define a new centered data set, X_ctd_new
X_ctd_new <- Fim_new - oneM %*% sqrt(Dj_new)
rownames(X_ctd_new) <- rownames(qttc_new)
colnames(X_ctd_new) <- colnames(qttc_new)

par(mfrow = c(1,2))
pcaX <- PCA(X_ctd_new,
            scale.unit=FALSE,
            ncp = ncol(qttc_new)-1, 
            ind.sup = NULL, quanti.sup = NULL, quali.sup = NULL, 
            row.w = fi_new, col.w = NULL, 
            graph = TRUE, axes = c(1,2))
pcaX$eig

# Number of significant dimensions
# 1/ Point successively to every component of the list of evals sorted by decreasing value 
# 2/ When total explained inertia exceeds 80% chose eval index as nbr of significant dims
n_dim <- ncol(X_ctd)-1
for (nd in 1:n_dim) {  if(pcaX$eig[nd,3] >= 80) break  }
cat("Number of significant dimensions:",nd)

par(mfrow = c(1,1))
# Conclusion: All checks. 
#   By way of conclusion explain in details certain points move so much.
