## ############################################
##  MIRI:     MVA
##  LAB #7:   Decision Trees
##  Authors:  Cedric Bhihe <cedric.bhihe@gmail.com>
##            Santi Calvo <s.calvo93@gmail.com>  
##  Delivery: before 2018.05.27 - 23:55
##  Script:   lab7-script_mva.R
## ############################################


rm(list=ls(all=TRUE))

# ############################################
## Environment and env. var.
# ############################################

options(scipen=6) # R switches to sci notation above 5 digits on plot axes
set.seed(932178)
setwd("~/Documents/Academic/UPC/MIRI/Subjects/MVA_multivariate-analysis/Labs/")
#setwd("C:/Users/calvo/Desktop/UPC/Courses/Third_semester/MVA/Exercises/Labs/")
ccolors=c("red","green","blue","orange","cyan","tan1","darkred","honeydew2","violetred",
          "palegreen3","peachpuff4","lavenderblush3","lightgray","lightsalmon","wheat2")
datestamp <- format(Sys.time(),"%Y%m%d-%H%M%S"); 



# ############################################
## Repos and libraries
# ############################################

setRepositories(ind = c(1:6,8))
#setRepositories()   # to specify repo on the fly:
#chooseCRANmirror()  # in case package(s) cannot be downloaded from default repo

#library(xlsx, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")

#install.packages(c("rpart", "rattle", "rpart.plot"))
library("rpart", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")        # tree building
library("rpart.plot", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")   # tree plotting
library("rattle", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")       # tree var importance

library("VIM", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
#library("mice", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")        # for imputations
#library(FactoMineR, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")    # to use PCA, CA and MCA method
require(graphics, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")       # enhanced graphics
library(ggplot2, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")        # to enhance graph plotting
library(ggrepel, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")        # to plot with well behaved labeling


# ############################################
## Functions
# ############################################
csvSaveF <- function(dataObj,targetfile) {
    write.table(dataObj,
                targetfile,
                append=F,
                sep=",",
                eol="\n",
                na ="NA",
                dec=".",
                row.names=F,
                col.names=T)
}    # save cvs to file on disk


# ############################################
## 1: Import file: audit.xlsx file
#     Convert to csv
# ############################################

# audit <- read.xlsx("Lab7/audit.xlsx",
#                  sheetIndex=1,
#                  sheetName=NULL, 
#                  header=T,
#                  as.data.frame=T, 
#                  encoding="UTF-8",
#                  keepFormulas=F)
# 
# rownames(audit) <- audit[,1]
# audit <- audit[,-1]
# targetfile <- "Data/audit.csv"
# csvSaveF(audit,targetfile)  # save to disk

# import csv file 
sourcefile <- "Data/audit.csv"
audit <- read.csv(sourcefile,
                  header=T,
                  quote = "\"", 
                  dec = ".",
                  sep=",", 
                  encoding="UTF-8",
                  check.names=TRUE,
                  stringsAsFactors = TRUE)

# inspect data for missings
dim(audit)
names(audit)
audit_aggr <- aggr(audit, 
                   numbers=TRUE,
                   bars=TRUE,
                   combined=FALSE,
                   prop=FALSE,
                   plot=TRUE,
                   axes=TRUE,
                   col=ccolors[6:7], 
                   labels=names(missing),
                   cex.axis=0.8,
                   ylab=c("Missing Data (count)","In-Variable Missing Data Distribution"))
# The majority of missings (100) concerns both the "Employment" and "Occupation" categorical variables
# Altogether 141 observations are affected by missings out of 2000 observation in total.
# We do not resort to imputation, because decision trees methods are known to be quite robust from the 
# point of view of missings and outliers. 
rm(audit_aggr)
# ############################################
## 2: Decide which predictors to use
#     Preprocess these variables if necessary
# ############################################

summary(audit)
audit_bak <- audit # backup

# We use observations' ID (1st column) as row names for the 2000 x 12 data set. 
# The data consists of:
# - 6 active categorical variables: c("Employment","Education", "Marital", "Occupation", "Gender", "Accounts")
names(audit_bak[which(sapply(audit_bak, is.factor))])
# - 4 active continuous variables:  c("Age","Income","Deductions","Hours")
names(audit_bak[which(! sapply(audit_bak[,1:10], is.factor))])
# - 1 supplementary continuous variable: c("Adjustment"), which is the result of the variable "Adjusted" being 1 
# - 1 supplementary categorical (yes/no or 0/1) variable: c("Adjusted"), which is our target/dependent variable

## discretize continuous variables
# Age
binsAge <- c(27,38,50)
modsAge <- c("Prime","Middle","Mature","Senior")
cat("Age:\n   below or equal to",binsAge[1],"years old - bin count:",length(which(audit$Age<=binsAge[1])),"\n")
cat("   between",binsAge[1]+1,"and",binsAge[2],"years old - bin count:", length(which(audit$Age>binsAge[1] & audit$Age<=binsAge[2])),"\n")
cat("   between",binsAge[2]+1,"and",binsAge[3],"years old - bin count:", length(which(audit$Age>binsAge[2] & audit$Age<=binsAge[3])),"\n")
cat("   above",binsAge[3],"years old - bin count:", length(which(audit$Age>binsAge[3])),"\n")

audit$Age[which(audit_bak$Age<=binsAge[1])] <- modsAge[1]
audit$Age[which(audit_bak$Age>binsAge[1] & audit_bak$Age<=binsAge[2])] <- modsAge[2] 
audit$Age[which(audit_bak$Age>binsAge[2] & audit_bak$Age<=binsAge[3])] <- modsAge[3]
audit$Age[which(audit_bak$Age>binsAge[3])] <- modsAge[4]

# Income
binsIncome <- c(34500,60000,115000)
modsIncome <- c("Low","Medium","High","Obscene")
cat("Income:\n   below or equal to",binsIncome[1],"USD/year - bin count:",length(which(audit$Income<=binsIncome[1])),"\n")
cat("   between",binsIncome[1]+1,"and",binsIncome[2],"USD/year - bin count:", length(which(audit$Income>binsIncome[1] & audit$Income<=binsIncome[2])),"\n")
cat("   between",binsIncome[2]+1,"and",binsIncome[3],"USD/year - bin count:", length(which(audit$Income>binsIncome[2] & audit$Income<=binsIncome[3])),"\n")
cat("   above",binsIncome[3]," USD/year - bin count:", length(which(audit$Income>binsIncome[3])),"\n")

audit$Income[which(audit_bak$Income<=binsIncome[1])] <- modsIncome[1]
audit$Income[which(audit_bak$Income>binsIncome[1] & audit_bak$Income<=binsIncome[2])] <- modsIncome[2] 
audit$Income[which(audit_bak$Income>binsIncome[2] & audit_bak$Income<=binsIncome[3])] <- modsIncome[3]
audit$Income[which(audit_bak$Income>binsIncome[3])] <- modsIncome[4]

# Deductions
cat("Deductions:\n   Audits with deductions:",length(which(audit$Deductions !=0)),paste0("(",100*length(which(audit$Deductions !=0))/ nrow(audit),"%)\n"))
max_y <- 300
histplot <- hist(audit$Deductions+1,
                 main="Deductions claimed", 
                 xlab="Amounts of Deduction (USD)",
                 ylab="Nbr of Deductions Claimants",
                 ylim=c(0,max_y),
                 border="blue", 
                 col="green",
                 las=1,
                 # breaks=10,
                 # proba=F,
                 freq=T)
text(x=(histplot$breaks[-1]+histplot$breaks[1:(length(histplot$breaks)-1)])/2,
     y=c(0.98*max_y,histplot$counts[-1]),
     labels=histplot$counts,
     col="blue",
     pos=3,
     cex=1)

histplot
rm(histplot)

# Hours
binsHours <- c(30,40)
modsHours <- c("Part-time","Reduced","Full")
cat("Hours:\n   below or equal to",binsHours[1],"- bin count:",length(which(audit$Hours<=binsHours[1])),"\n")
cat("   between",binsHours[1]+1,"and",binsHours[2],"- bin count:", length(which(audit$Hours>binsHours[1] & audit$Hours<=binsHours[2])),"\n")
cat("   above",binsHours[2],"hours/week - bin count:", length(which(audit$Hours>binsHours[2])),"\n")

audit$Hours[which(audit_bak$Hours<=binsHours[1])] <- modsHours[1]
audit$Hours[which(audit_bak$Hours>binsHours[1] & audit_bak$Hours<=binsHours[2])] <- modsHours[2]
audit$Hours[which(audit_bak$Hours>binsHours[2])] <- modsHours[3]

# ############################################
## 3: Select the last third of observations as testing data set
# ############################################

trainRows <- 1:floor(nrow(audit)*2/3)             # used to select optimal model
testRows <- ceiling(nrow(audit)*2/3):nrow(audit)  # used to realize model


# ############################################
## 4: Build the decision tree to predict var “Adjusted” using the training data. 
#     Determine cutoff value for decision making
# ############################################

#audit$Adjusted <- as.factor(audit$Adjusted)

## Keep  var "Deductions" to compare results and decision trees
auditTree1 <- rpart(Adjusted ~., 
                   data=audit[trainRows,-c(11)],  # take out illustrative var. "Adjustment"
                   weights=NULL,
                   method="class",
                   control=rpart.control(cp=0.001, xval=10))
# cp = complexity parameter = min value of alpha, when growing the initial tree
# xval = nbr of CV runs performed on training data

# plot raw tree
par(mfrow = c(1,1), xpd = NA)
plot(auditTree1)
text(auditTree1,use.n=T,cex=0.8,col="darkblue")
#par(mfrow = c(1,1))

plotcp(auditTree1)             # visualize cross validation results

# plot Relative tree impurity, R(T), vs Split number, aka Tree size
plot(auditTree1$cptable[,2],auditTree1$cptable[,3],
     type="l",
     xlab="size of tree",
     ylab="relative tree impurity",
     main="R(T)",
     col=ccolors[1])
lines(auditTree1$cptable[,2],auditTree1$cptable[,4],
      col=ccolors[2])
legend("topright",
       c("R(T) training","R(T) mean CV"),
       col=ccolors[1:2],
       lty=1)

# find the tree with the minimum CV error, minCVerror_idx
auditTree1$cptable # print complexity parameter table, where nsplit= nbr nodes -1

auditTree1$cptable <- as.data.frame(auditTree1$cptable)
minCVerror_idx <- which.min(auditTree1$cptable$xerror)
# print complexity parameter table
printcp(auditTree1,digits=3)
# find corresponding mean(CV_error) and sd(CV_error)
minCVerr_mean <- auditTree1$cptable$xerror[minCVerror_idx]
minCVerr_std <- auditTree1$cptable$xstd[minCVerror_idx]

# from root, find first tree whose mean(CV_error) <= minCVerr_mean + 1* minCVerr_std
i = 1
while(auditTree1$cptable$xerror[i] > (minCVerr_mean + minCVerr_std)) { i <- i+1 } 

# post-pruning of Tree
alpha = auditTree1$cptable$CP[i]
auditTree1_optimum <- prune(auditTree1,cp=alpha)

# plot and rules
rpart.plot(auditTree1_optimum, main="Optimal Tree")
asRules(auditTree1_optimum)

par(mfrow = c(1,1), xpd = NA)
plot(auditTree1_optimum)
text(auditTree1_optimum,use.n=T,cex=0.8,col="darkblue")

## Get rid of var "Deductions" to compare results and decision trees
auditTree2 <- rpart(Adjusted ~., 
                    data=audit[trainRows,-c(8,11)],  # take out illustrative var. "Adjustment", keep "Deductions"
                    weights=NULL,
                    method="class",
                    control=rpart.control(cp=0.001, xval=10))
# cp = complexity parameter = min value of alpha, when growing the initial tree
# xval = nbr of CV runs performed on training data

# plot raw tree
par(mfrow = c(1,1), xpd = NA)
plot(auditTree2)
text(auditTree2,use.n=T,cex=0.8,col="darkblue")
#par(mfrow = c(1,1))

plotcp(auditTree2)             # visualize cross validation results

# plot Relative tree impurity, R(T), vs Split number, aka Tree size
plot(auditTree2$cptable[,2],auditTree$cptable[,3],
     type="l",
     xlab="size of tree",
     ylab="relative tree impurity",
     main="R(T)",
     col=ccolors[1])
lines(auditTree2$cptable[,2],auditTree$cptable[,4],
      col=ccolors[2])
legend("topright",
       c("R(T) training","R(T) mean CV"),
       col=ccolors[1:2],
       lty=1)

# find the tree with the minimum CV error, minCVerror_idx
auditTree2$cptable <- as.data.frame(auditTree2$cptable)
minCVerror_idx <- which.min(auditTree2$cptable$xerror)
# print complexity parameter table
printcp(auditTree2,digits=3)
# find corresponding mean(CV_error) and sd(CV_error)
minCVerr_mean <- auditTree2$cptable$xerror[minCVerror_idx]
minCVerr_std <- auditTree2$cptable$xstd[minCVerror_idx]

# from root, find first tree whose mean(CV_error) <= minCVerr_mean + 1* minCVerr_std
i = 1
if(auditTree2$cptable$xerror[i] <= minCVerr_mean + minCVerr_std) {
    break
} else {
    i <- i+1
}

# post-pruning of Tree
alpha = auditTree2$cptable$CP[i]
auditTree2_optimum <- prune(auditTree2,cp=alfa)

# plot and rules
rpart.plot(auditTree2_optimum, main="Optimal Tree")
asRules(auditTree2_optimum)

# ############################################
## 5: Plot the importance of variables in the prediction.
# ############################################

auditTree1_optimum$variable.importance


# ############################################
## 6: Compute the accuracy, precision, recall and AUC on test data.
# ############################################




# ############################################
## 7:  Perform a Random Forest on test data
# ############################################

