##  Script: wparam-arbitrary.R
##  Author:   Cedric Bhihe <cedric.bhihe@gmail.com>

rm(list=ls(all=TRUE))
set.seed(932178)
#setwd("C:/Users/calvo/Desktop/UPC/Courses/Third_semester/MVA/Exercises/Lab1v2/")
setwd("~/Documents/Study/UPC/MIRI-subjects/MVA_multivariate-analysis/Labs/")

# Load data
list_files = as.vector(grep("^russet.*\\mpp.csv$",
                            list.files(path="Data",full.names=FALSE),
                            ignore.case = TRUE,
                            perl = TRUE,
                            fixed = FALSE,
                            inv=FALSE,
                            value=TRUE))
russet_imp_data = read.csv(paste0("Data/",list_files[1]),
                           header = TRUE, 
                           quote = "\"", 
                           dec = ".",
                           sep=",", 
                           check.names=TRUE)
dim(russet_imp_data)  # list nbr of observations followed by nbr of vars
russet_imp_data <- russet_imp_data[,-9]

# Toy example
wparam <- rep(1:5,10)
length(wparam) <- 47  # nrow(X)

# "Cuba" weight <- 0
cuba_row <- which(rownames(russet_imp_data)=="Cuba")
wparam <- c(rep(1,cuba_row-1),0,rep(1,nrow(russet_imp_data)-cuba_row))

# ##############################################
# Write data to disc
wparam <- as.matrix(wparam,nrow=47, byrow=T)
colnames(wparam) <- c("obs_weight")
datestamp=format(Sys.time(),"%Y%m%d-%H%M%S")
write.csv(wparam, file=sprintf("Data/%s_arbitrary-wparam.csv",datestamp))
