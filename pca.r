# author : Sandali Lokuge
# date : 12/11/2022
# get the PCA plots of data - before and after normalization
# set the working directory to main (project) file location

library(ggplot2)
library(ggfortify)

# set the working directory as the main project directory - setwd()

################################ functions ################################

plotPCA <- function(dataFile,typeCount,datasetNames){
  data <- read.csv(dataFile,sep='\t')
  data_woProbe <- data[,-1]
  data_tp <- as.data.frame(t(data_woProbe))
  
  # give the corresponding dataset names
  datasets<- c(rep(c(datasetNames[1]),times=c(typeCount[1])),
                  rep(c(datasetNames[2]),times=c(typeCount[3])),
                  rep(c(datasetNames[3]),times=c(typeCount[5])),
                  rep(c(datasetNames[4]),times=c(typeCount[2])),
                  rep(c(datasetNames[5]),times=c(typeCount[4])),
                  rep(c(datasetNames[6]),times=c(typeCount[6])))

  types <-c(rep(c("tumor"),times=c(typeCount[1]+typeCount[3]+typeCount[5])),
                rep(c("non-tumor"),times=c(typeCount[2]+typeCount[4]+typeCount[6])))
  
  df_type <- data.frame(types,datasets)
  
  autoplot(prcomp(data_tp, center = TRUE),df_type, colour="datasets", label = F,text.size = 15)
  
}

##########################################################################

# change these details according to the datasets
cancerType <- "Ovarian cancer"
dataset1_tumor <- 12
dataset1_nontumor <- 12
dataset2_tumor <- 4
dataset2_nontumor <- 4
dataset3_tumor <- 18
dataset3_nontumor <- 12

typeCount <- c(dataset1_tumor,dataset1_nontumor,dataset2_tumor,dataset2_nontumor,dataset3_tumor,dataset3_nontumor)
datasetNames <- c("GSE14407_tumor","GSE36668_tumor","GSE38666_tumor","GSE14407_non-tumor","GSE36668_non-tumor","GSE38666_non-tumor")


root_path <- getwd()
folder_path <- paste(root_path,"Expression_data",cancerType,sep="/")

# PCA for before normalization
plotPCA(paste(folder_path,"Ovarian cancer_beforeNorm.csv",sep="/"),typeCount,datasetNames)

# PCA for after normalization
plotPCA(paste(folder_path,"Ovarian cancer_afterNorm.csv",sep="/"),typeCount,datasetNames)


