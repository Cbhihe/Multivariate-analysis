## ############################################
##  Topic:    Association Rules
##  Authors:  Cedric Bhihe, Santi Calvo 
##  Date: 2018.06.15
##  Script:   script-ar.R
## ############################################


rm(list=ls(all=TRUE))

setwd("~/Documents/Academic/UPC/MIRI/Subjects/MVA_multivariate-analysis/Labs/")
# setwd("C:/Users/calvo/Desktop/UPC/Courses/Third_semester/MVA/mva-labs/")
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
csvSaveF <- function(dataObj,targetFile) {
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
tictt2 <- tictt1
tictt2[sapply(tictt1,is.logical)] <- lapply(tictt1[sapply(tictt1,is.logical)],as.factor)
str(tictt2)                       # verify conversion
length(unique(names(tictt2)))     # check =33

# inspect data for missings
aggr(tictt2, 
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
tictt2_catdes <- catdes(tictt2,
                        28,
                        proba=0.05,
                        row.w=NULL)


# ############################################
## 3: Transform `tictt` in a transaction file
# ############################################
ticttr <- as(tictt1,"transactions")    # ensure 'arules' is correctly loaded before

# Inspect the transctions variable
inspect(head(ticttr, n = 5))

# ############################################
## 4: Define itemsets' parameters: min_support, min_confidence, max_size,
#     Run the apriori function.
# ############################################

min_support <- min(itemFrequency(ticttr))
sort(itemFrequency(ticttr), decreasing=F)[1:10]
length(itemFrequency(ticttr))
# itemFrequencyPlot(ticttr)
ticTransRules = apriori(ticttr)  # execute with default parameter values:
                                 # support=.1,confidence=.8,maxlen=10,maxtime=5
summary(ticTransRules)           # check that the maximum item length is 9
attributes(ticTransRules)
class(ticTransRules)


# We choose:
min_support <- .01
min_confidence <- 0.75
max_size <- 5

ticTransRules = apriori(ticttr,parameter=list(support=min_support,
                                              confidence=min_confidence,
                                              maxlen=max_size,
                                              maxtime=10))
summary(ticTransRules)
inspect(head(ticTransRules, n = 10, by ="confidence",decreasing=F))
inspect(head(ticTransRules, n = 10, by ="lift",decreasing=T))
# ticTR <- subset(ticTransRules, subset = lhs != "")   # impossible to apply for object of type "rules"


# ############################################
## 5. List the 10 most frequent item sets (i.e. rules sorted by support)
# ############################################
inspect(head(ticTransRules, n = 10, by = "support",decreasing=T))


# ############################################
## 6. List the first 10 rules sorted by lift
# ############################################
inspect(head(ticTransRules, n = 10, by ="lift",decreasing=T))


# ############################################
## 7. List the first 10 rules sorted by lift with "Pagament.a.través.d.Internet." as Consequent
# ############################################
ticTR <- subset(ticTransRules, subset = rhs %in% c("Pagament.a.través.d.Internet."))
summary(ticTR)
inspect(head(ticTR,n = 10, by="lift", decreasing=T))
