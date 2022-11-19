# author : Sandali Lokuge
# date : 12/11/2022
# combine the three geo datasets 
# batch correction using Combat
# class specific quantile normalization
# set the working directory to main (project) file location

library(Biobase)
library(sva)
library("tibble")
library(preprocessCore)
library(dplyr)

# set the working directory as the main project directory - setwd()

################################ functions ################################

# combine three geo datasets
combineDatasets <- function(typeCount){
  
  df_list <- list()
  i <- 1
  for (geo_dataFile in geo_datasets){
    data <- read.csv(paste(folder_path,geo_dataFile,sep="/"),sep='\t')
    df_list[[i]] <- data
    i <- i+1
  }
  expdf_merge <- merge(df_list[[1]],df_list[[2]],by="Probe_ID")
  expdf_merge <- merge(expdf_merge,df_list[[3]],by="Probe_ID")
  
  # rearrange the tumor and non-tumor data
  tumor_df <- cbind(expdf_merge[,2:(2+typeCount[1]-1)],
                    expdf_merge[,(2+typeCount[1]+typeCount[2]):(2+typeCount[1]+typeCount[2]+typeCount[3]-1)],
                    expdf_merge[,(2+typeCount[1]+typeCount[2]+typeCount[3]+typeCount[4]):(2+typeCount[1]+typeCount[2]+typeCount[3]+typeCount[4]+typeCount[5]-1)])
  
  nontumor_df <- cbind(expdf_merge[,(2+typeCount[1]):(2+typeCount[1]+typeCount[2]-1)],
                       expdf_merge[,(2+typeCount[1]+typeCount[2]+typeCount[3]):(2+typeCount[1]+typeCount[2]+typeCount[3]+typeCount[4]-1)],
                       expdf_merge[,(2+typeCount[1]+typeCount[2]+typeCount[3]+typeCount[4]+typeCount[5]):(2+typeCount[1]+typeCount[2]+typeCount[3]+typeCount[4]+typeCount[5]+typeCount[6]-1)])
  
  res_df <- cbind(expdf_merge[1],tumor_df,nontumor_df)
  
  fileName <- paste(cancerType,"beforeNorm.csv",sep="_")
  write.table(res_df,paste(folder_path,fileName,sep="/"),sep="\t",row.names = FALSE)
  
  return(res_df)
}

# batch correction
batchEffectCorrection <- function(data,typeCount){
  df <- data[,2:ncol(data)]
  colnames_df <- colnames(df) 
  
  ex<-new("ExpressionSet",exprs=as.matrix(df))
  edata<-exprs(ex) 
  df<-as.data.frame(t(df))

  df["platform"]<- c(rep(c("GPL570_1"),times=c(typeCount[1])),
                     rep(c("GPL570_2"),times=c(typeCount[3])),
                     rep(c("GPL570_3"),times=c(typeCount[5])),
                     rep(c("GPL570_1"),times=c(typeCount[2])),
                     rep(c("GPL570_2"),times=c(typeCount[4])),
                     rep(c("GPL570_3"),times=c(typeCount[6])))
  df$platform=factor(df$platform)
  
  df["type"]<-c(rep(c("tumor"),times=c(typeCount[1]+typeCount[3]+typeCount[5])),
                rep(c("non-tumor"),times=c(typeCount[2]+typeCount[4]+typeCount[6])))
  df$type=factor(df$type)
  
  # remove Batch effect
  mod=model.matrix(~df$type)
  batch=df$platform
  batch=as.numeric(batch)
  cleandata<-ComBat(edata,batch,mod)
  
  new_cleandata <- as.data.frame(cleandata)
  new_cleandata <- add_column(new_cleandata, data['Probe_ID'], .before = 1)
  
  return(new_cleandata)
  
}

# class specific quantile normalization
quantileNormalization <- function(batchCorr_data,typeCount){
  
  tumor_df <- batchCorr_data[,2:(2+typeCount[1]+typeCount[3]+typeCount[5]-1)]
  nontumor_df <- batchCorr_data[,(2+typeCount[1]+typeCount[3]+typeCount[5]):(2+sum(typeCount)-1)]
  
  # class-specific quantile normalization
  tumor_df_norm <- as.data.frame(normalize.quantiles(as.matrix(tumor_df)))
  nontumor_df_norm <- as.data.frame(normalize.quantiles(as.matrix(nontumor_df)))
  
  res_df <- cbind(batchCorr_data[1],tumor_df_norm,nontumor_df_norm)
  colnames(res_df) <- c("Probe_ID",colnames(tumor_df),colnames(nontumor_df))
  
  return(res_df)
}

##########################################################################

#change these details according to the datasets
cancerType <- "Ovarian cancer"
dataset1_tumor <- 12
dataset1_nontumor <- 12
dataset2_tumor <- 4
dataset2_nontumor <- 4
dataset3_tumor <- 18
dataset3_nontumor <- 12

typeCount <- c(dataset1_tumor,dataset1_nontumor,dataset2_tumor,dataset2_nontumor,dataset3_tumor,dataset3_nontumor)

root_path <- getwd()
folder_path <- paste(root_path,"Expression_data",cancerType,sep="/")
geo_datasets <- list.files(path = folder_path)

# combine three datasets
data <- combineDatasets(typeCount)

# batch correction
batchCorr_data <- batchEffectCorrection(data,typeCount)

# class specific quantile normalization
normalized_data <- quantileNormalization(batchCorr_data,typeCount)

# save to csv file
fileName <- paste(cancerType,"afterNorm.csv",sep="_")
write.table(normalized_data,paste(folder_path,fileName,sep="/"),sep="\t",row.names = FALSE)


