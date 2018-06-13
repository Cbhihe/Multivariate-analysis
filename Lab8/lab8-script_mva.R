## ############################################
##  Lab #8:   Association Rules
##  Authors:  Cedric Bhihe, Santi Calvo 
##  Delivery: before 2018.06.__ - 23:55
##  Script:   script-ar.R
## ############################################


rm(list=ls(all=TRUE))

setwd("~/Documents/Academic/UPC/MIRI/Subjects/MVA_multivariate-analysis/Labs/")
set.seed(932178)

options(scipen=6) # R switches to sci notation above 5 digits on plot axes
ccolors=c("blue","red","green","orange","cyan","tan1","darkred","honeydew2","violetred",
          "palegreen3","peachpuff4","lavenderblush3","lightgray","lightsalmon","wheat2")

datestamp <- format(Sys.time(),"%Y%m%d-%H%M%S"); 


# ############################################
## Repos and libraries
# ############################################

setRepositories(ind = c(1:6,8))
#setRepositories()   # to specify repo on the fly:
#chooseCRANmirror()  # in case package(s) cannot be downloaded from default repo

library(arules, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")          # transaction rules
library(VIM, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")             # imputation and aggr() to view missings
library(FactoMineR, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")      # PCA, catgorical descriptions w/ catdes()
require(graphics, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")        # enhanced graphics
# library(ggplot2, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")         # to enhance graph plotting
# library(ggrepel, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")         # to plot with well behaved labeling


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
                row.names=T,
                col.names=T)
}    # save cvs to file on disk


# ############################################
## 1: Import file: Data/tic_tt.txt 
#     Check each variable's class
# ############################################

tictt1 <- read.csv("Data/tic_tt.txt",
                 header=T,
                 #encoding="UTF-8",
                 fileEncoding = "ISO-8859-1",
                 sep=";")

rownames(tictt1) <- tictt1[,1]
tictt1 <- tictt1[,-1]
classVar <- lapply(tictt1,class)   # class of each variable

# ############################################
## 2: Find the profiles of people who practice payments on the Internet
# ############################################

## convert every non-factor (logical) variable in a factor
tictt1[sapply(tictt1,is.logical)] <- lapply(tictt1[sapply(tictt1,is.logical)],as.factor)
str(tictt1)                       # verify conversion
length(unique(names(tictt1)))     # check =33

# inspect data for missings
aggr(tictt1, 
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

# lookup profile of people who buy in internet.
tictt1_catdes <- catdes(tictt1,
                        28,
                        proba=0.05,
                        row.w=NULL)


# ############################################
## 3: Transform `tictt` in a transaction file
# ############################################
tictt2 <- read.csv("Data/tic_tt.txt",
                   header=T,
                   #encoding="UTF-8",
                   fileEncoding = "ISO-8859-1",
                   sep=";")
rownames(tictt2) <- tictt2[,1]
tictt2 <- tictt2[,-1]
ticttr <- as(tictt2,"transactions")    # ensure 'arules' is correctly loaded before
rm(tictt2)


# ############################################
## 4: Define itemsets' parameters: min_support, min_confidence, max_size,
#     Run the apriori function.
# ############################################

min_support <- min(itemFrequency(ticttr))
sort(itemFrequency(ticttr), decreasing=F)[1:10]
length(itemFrequency(ticttr))
# itemFrequencyPlot(ticttr)
min_confidence <- min()
max_size <- max()
