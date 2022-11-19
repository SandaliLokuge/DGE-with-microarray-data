# author : Sandali Lokuge
# date : 12/11/2022
# read the CEL files and save the expression data in csv files
# set the working directory to main (project) file location
# the geo datsets should be separated into tumor and non-tumor conditions

library(limma)
library(affy)
library("annotate")
library("hgu133plus2.db") 

# set the working directory as the main project directory - setwd()

cancerType <- "Ovarian cancer"
root_path <- getwd()

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

geo_files <- list.dirs(path = paste(root_path,"microarray data",cancerType,sep="/"), full.names = TRUE, recursive = FALSE)

for (geo_file in geo_files){
  # get tumor and non-tumor CEL files separately 
  condition_files <- list.dirs(path =geo_file, full.names = TRUE, recursive = FALSE)
  
  expdf_list <- list()
  i <- 1
  
  for (condition_file in condition_files){
    # Set the path
    setwd(condition_file)
    
    # Reads the set of CEL files
    raw_data <- ReadAffy()
    
    # RMA normalization on our expression data across all samples
    rmaNorm_data<- affy::rma(raw_data)
    
    # Get the expression estimates for each array
    exp_data <- exprs(rmaNorm_data)
    
    exp_df <- as.data.frame(exp_data)
    exp_df <- cbind(Probe_ID = rownames(exp_df), exp_df)
    rownames(exp_df) <- 1:nrow(exp_df)
    
    expdf_list[[i]] <- exp_df
    i <- i+1
    
  }
  
  # merge tumor and non-tumor datasets
  expdf_merge <- merge(expdf_list[[2]],expdf_list[[1]],by="Probe_ID")
  # save the expression data
  filename <- unlist(strsplit(split_path(geo_file)[1], split= "_"))[1]
  filename <- paste(filename,"Exp.csv",sep="_")

  write.table(expdf_merge,paste(root_path,"Expression_data",cancerType,filename,sep="/"),sep="\t",row.names = FALSE)
  
}


