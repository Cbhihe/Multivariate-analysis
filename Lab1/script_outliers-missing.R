# #####################################
## Topic:    Missing values and outliers
## Authors:  Cedric Bhihe, Santi Calvo
## Date:     2018.03.07 - 23:55
## Script:   script_outliers-missing.R
# #####################################


rm(list=ls(all=TRUE))
library("mice")
library("DMwR")
library("VIM")    # exec `install.packages("VIM")` in R shell first
library("dplyr")  # exec `install.packages("dplyr")` in R shell first
                  # split-apply-combine, specialized to df
library(chemometrics)  # exec `install.packages("chemometrics")` in R shell first
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
# Load data and inspect it
#############################################
list_files = as.vector(grep("^russet.*\\.csv$",
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
dim(russet_data)  # list nbr of observations followed by nbr of vars
class(russet_data)  # check that data set is of class 'dataframe'
summary(russet_data)

# Check missing ('NA') values in each component of the dataframe 'russet_data'
russet_aggr <- aggr(russet_data, 
                    numbers=TRUE,
                    bars=TRUE,
                    combined=FALSE,
                    prop=FALSE,
                    plot=TRUE,
                    axes=TRUE,
                    col=ccolors[6:7], 
                    labels=names(russet_data),
                    cex.axis=1,
                    ylab=c("Missing Data (count)","In-Variable Missing Data Distribution")
)

summary(russet_aggr)
names(russet_data)
#russet_aggr$combinations
russet_aggr$count
russet_aggr$percent
russet_aggr$missings
russet_aggr$tabcomb

#############################################
# Impute missing values in data + compare several imputations
#############################################

# Extract individuals concerned by imputations
cnt_missing <- apply(russet_data, 1, function(x) sum(is.na(x)))
russet_data <- cbind(russet_data,cnt_missing)
russet_data <- russet_data[1:11]

# Generate table of missing values in data set
country_missings <- cbind(russet_data[russet_data$cnt_missing != 0,][1],
                          russet_data[russet_data$cnt_missing != 0,][11])
first_col <- rownames(country_missings)

country_missings <- as.matrix.data.frame(country_missings)
country_missings <- cbind(first_col,as.character(country_missings[,1]),country_missings[,2],rep("",length(first_col)))
nbr_col=ncol(country_missings)

jj=0
for (index in country_missings[,1]) {
  jj <- jj+1
  missing_var = c()
  ii= as.numeric(index)
  missing_var <- names(russet_data[ii,])[is.na(russet_data[ii,])]
  country_missings[jj,nbr_col] <- paste(missing_var,collapse=", ") 
}

colnames(country_missings) <- c("df_row_index","Country","Nbr_missing","Var_missing")
rownames(country_missings) <- NULL

# Save table of missing values to disk
todisk_file <- sprintf("Lab1//Report//russet_missing-values.txt")
write.table(country_missings,file=todisk_file, append=FALSE, sep = "\t", eol="\n")


# Impute non assigned values using kNN based on Euclidian distance 
# of k nearest neighbors, for k in {1,3,5,7} 

# Dimension output matrix from nbr of kNN imputation trials
seq_lb=1;seq_ub=7;seq_inter=2

country_impute <- cbind(country_missings[,2],country_missings[,4])
column_names=c("Country","Var_missing")

# Imputation trials, kNN
for (kk in seq(seq_lb,seq_ub,seq_inter)) {
  column_names <- c(column_names,paste0("k=",kk))
  data_imp <- kNN(russet_data,
                  variable=colnames(russet_data[2:10]),
                  dist_var=colnames(russet_data[2:10]),
                  k=kk,
                  trace=TRUE,
                  useImputedDist=FALSE,
                  weightDist=FALSE)
  
  inter_df1 <-data_imp[data_imp$pais == country_missings[1,2],]
  inter_df2 <-data_imp[data_imp$pais == country_missings[2,2],]
  inter_df3 <-data_imp[data_imp$pais == country_missings[3,2],]
  inter_df4 <-data_imp[data_imp$pais == country_missings[4,2],]
  
  country_impute <- cbind(country_impute,c(as.numeric(inter_df1[names(inter_df1) == country_missings[1,4]]),
                                           as.numeric(inter_df2[names(inter_df2) == country_missings[2,4]]),
                                           as.numeric(inter_df3[names(inter_df3) == country_missings[3,4]]),
                                           as.numeric(inter_df4[names(inter_df4) == country_missings[4,4]]))
  )
}
colnames(country_impute) <- column_names

# Save table of imputed values to disk
todisk_file <- sprintf("Lab1//Report//russet_imputed-values_knn-trials.txt")
write.table(country_impute,file=todisk_file, append=FALSE, sep = "\t", eol="\n")

# >>> CHOOSE for instance k=5 for imputations
data_imp_knn <- kNN(russet_data,
                    variable=colnames(russet_data[2:10]),
                    dist_var=colnames(russet_data[2:10]),
                    k=5,
                    trace=TRUE,
                    useImputedDist=FALSE,
                    weightDist=FALSE)
data_imputed_knn <- as.matrix(data_imp_knn[,2:10])
rownames(data_imputed_knn) <- russet_data[,1]
class(data_imputed_knn) <- "numeric"

# Save table of imputed values for k=5 to disk
todisk_file <- sprintf("Lab1//Report//russet_imputed-values_knn-k5.txt")
write.table(data_imputed_knn,file=todisk_file, append=FALSE, sep = "\t", eol="\n")

# MICE imputation
data_imp_mice <- mice(russet_data[2:10],
                      where=is.na(russet_data[2:10]),
                      method="pmm",
                      m = 8,
                      maxint=10);
summary(data_imp_mice)
data_imputed_mice_df <- complete(data_imp_mice)
data_imputed_mice <- as.matrix(data_imputed_mice_df)
rownames(data_imputed_mice) <- russet_data[,1]
# Save table of imputed values for k=5 to disk
todisk_file <- sprintf("Lab1//Report//russet_imputed-values_mice-mpp.txt")
write.table(data_imputed_mice,
            file=todisk_file, 
            append=FALSE, 
            sep="\t", 
            eol="\n")


########################################
# Outliers analysis
########################################
fromdisk_file <- sprintf("Lab1//Report//russet_imputed-values_mice-mpp.txt")
data_imputed <- read.table(file=fromdisk_file,
                           header=TRUE,
                           sep="\t",
                           dec=".")

# Method 1: hierarchical clustering (HC)
HCoutlier <- outliers.ranking(data_imputed,
                              test.data=NULL,
                              method="sizeDiff",
                              clus=list(dist = "euclidean",alg = "hclust",meth = "average"),
                              power = 1, 
                              verb = F)
# Try for very different results (in the above) for the hierarchical clustering method: 
#    "ward.D" (instead of "average")

# top most outlier ranks ordered by decreasing score of outlyingness
HCoutlier$prob.outliers[HCoutlier$rank.outliers][1:3]  # 3 first potential candidates
HCoutlier$prob.outliers[HCoutlier$prob.outliers >=0.90]  # potential outliers w/ prob score >= 90%


HCoutlier <- outliers.ranking(data_imputed,
                              test.data=NULL,
                              method="sizeDiff",
                              clus=list(dist = "euclidean",alg = "hclust",meth = "ward.D"),
                              power = 1, 
                              verb = F)
# Try for very different results (in the above) for the hierarchical clustering method: 
#    "ward.D" (instead of "average")

# outlier ranking factors ordered by decreasing score of outlyingness
HCoutlier$prob.outliers[HCoutlier$rank.outliers][1:3]  # 3 first potential candidates
HCoutlier$prob.outliers[HCoutlier$prob.outliers >=0.90]  # potential outliers w/ prob score >= 90%


# Method 2: local outlier factor (LOF density based local detection)
#   for 3 and 5 neighbors
LOFoutlier_k3 <- lofactor(data_imputed, k=3)
LOFoutlier_k5 <- lofactor(data_imputed, k=5)
plot(density(LOFoutlier_k3))
# pick top outliers
outliers_k3 <- order(LOFoutlier_k3, decreasing=T)[1:5]
outliers_k5 <- order(LOFoutlier_k5, decreasing=T)[1:5]
if (identical(sort(outliers_k3),sort(outliers_k3))) {
  LOFtab <- rbind(round(LOFoutlier_k3[outliers_k3],3),
                  round(LOFoutlier_k5[outliers_k3],3)
                 )
  colnames(LOFtab) <- as.character(russet_data[outliers_k3,1])
# } else {
  # LOFtab <- rbind(LOFoutlier_k3[outliers_k3],
  #                 distribute LOFoutlier_k5 according to col-names below    
  #                )
  # colnames(LOFtab) <- as.character(russet_data[unique(union(outliers_k3,outliers_k5)),1])
}
rownames(LOFtab) <- c("k=3","k=5")
LOFtab

# Method 3: Mahalanobis distance outlier detection
MDoutliers <- Moutlier(data_imputed, quantile = 0.975, plot = TRUE)
# sorted outliers for MCD with rank .ge. cutoff
sort(round(MDoutliers$md[MDoutliers$md >= MDoutliers$cutoff],3),decreasing=T)
# 5 top most ranking outliers for MRD
MRD_index_ordered <- order(MDoutliers$rd, decreasing=T)
round(MDoutliers$rd[MRD_index_ordered][1:5],3)

##########################################
# Plots
##########################################

par(mfrow = c(1, 3))
index <- seq(1:length(russet_data[,1]))

# Hierarchical clustering
HC_index_ordered <- order(HCoutlier$prob.outliers, decreasing=T)
HC_df <- cbind.data.frame(index, round(HCoutlier$prob.outliers,3), russet_data$pais)
colnames(HC_df) <- c("Index", "HC_rank", "Country")
HC_plot <- plot(HCoutlier$prob.outliers, 
                 pch="o", 
                 cex=1, 
                 main="Potential HC outliers\n by hierarchical clustering (HC-Ward)",
                 ylab="HC Rank") 
HC_plot_cutoff <- 0.5*(HC_df[HC_index_ordered[4],]$HC_rank + HC_df[HC_index_ordered[5],]$HC_rank)
abline(h = HC_plot_cutoff, col="red")  # add cutoff line
text(x=1:length(index), 
     y=HC_df$HC_rank, 
     labels=ifelse(HC_df$HC_rank >= HC_df[HC_index_ordered[4],]$HC_rank,as.character(HC_df$Country),""),
     pos=4,
     col="red")  # add labels

# LOF
LOF_index_ordered <- order(LOFoutlier_k5, decreasing=T)
LOF_df <- cbind.data.frame(index, round(LOFoutlier_k5,3), russet_data$pais)
colnames(LOF_df) <- c("Index", "LOF_rank", "Country")
LOF_plot <- plot(LOFoutlier_k5, 
                pch="o", 
                cex=1, 
                main="Potential LOF outliers\n by local outliers factor analysis (LOF-k=5)",
                ylab="LOF Rank") 
LOF_plot_cutoff <- 0.5*(LOF_df[LOF_index_ordered[4],]$LOF_rank + LOF_df[LOF_index_ordered[5],]$LOF_rank)
abline(h = LOF_plot_cutoff, col="red")  # add cutoff line
text(x=1:length(index), 
     y=LOF_df$LOF_rank, 
     labels=ifelse(LOF_df$LOF_rank >= LOF_df[LOF_index_ordered[4],]$LOF_rank,as.character(LOF_df$Country),""),
     pos=4,
     col="red")  # add labels

# Mahalanobis distance
MRD_index_ordered <- order(MDoutliers$rd, decreasing=T)
MRD_df <- cbind.data.frame(index, round(MDoutliers$rd,3), russet_data$pais)
colnames(MRD_df) <- c("Index", "MRD_rank", "Country")
#MRD_plot <- ggplot(MRD_df, aes(Index,MRD), main="Malahanobis Robust Distance")
MRD_plot <- plot(MDoutliers$rd, 
                 pch="o", 
                 cex=1, 
                 main="Potential MRD outliers\n by Mahalanobis robust distance (MRD)",
                 ylab="MRD Rank") 
MRD_plot_cutoff <- 0.5*(MRD_df[MRD_index_ordered[4],]$MRD_rank + MRD_df[MRD_index_ordered[5],]$MRD_rank)
abline(h = MRD_plot_cutoff, col="red")  # add cutoff line
text(x=1:length(index), 
     y=MRD_df$MRD_rank, 
     labels=ifelse(MRD_df$MRD_rank >= MRD_df[MRD_index_ordered[4],]$MRD_rank,as.character(MRD_df$Country),""),
     pos=4,
     col="red")  # add labels
