# author : Sandali Lokuge
# date : 12/11/2022
# find the gene symbols related to probes

library("annotate")
library("hgu133plus2.db") 
library("dplyr") 
library("stringr")
library("tidyverse")

# set the working directory as the main project directory - setwd()

################################ functions ################################
probe_to_gene <- function(readfile){
  
  annotationTable <- as.data.frame(AnnotationDbi::select(hgu133plus2.db, c(pull(readfile, Probe_ID)), c("SYMBOL","ENTREZID", "GENENAME")))
  
  # remove probes having two genes
  test <- annotationTable %>% group_by(PROBEID) %>% summarise(SYMBOL = toString(list(SYMBOL)))
  res_df <- add_column(readfile, test['SYMBOL'], .after = 1)
  res_df2 <- res_df %>% filter(!str_detect(SYMBOL, fixed("c(")))
  res_df2 <- res_df2[,2:ncol(res_df2)]
  
  return(res_df2)
  
}
cancerType <- "Ovarian cancer"

root_path <- getwd()
folder_path <- paste(root_path,"Limma_outputs",cancerType,sep="/")
outputFolder_path <- paste(root_path,"prob_annotation",cancerType,sep="/")

readfile <- read.csv(paste(folder_path,"upRegulated_genes.csv",sep="/"),sep='\t')
annotation_df <-probe_to_gene(readfile)
outputFile <- paste(outputFolder_path,"upRegulated_genes_Annotation.csv",sep="/")
write.table(annotation_df,outputFile, sep="\t",row.names = FALSE)

readfile <- read.csv(paste(folder_path,"downRegulated_genes.csv",sep="/"),sep='\t')
annotation_df <-probe_to_gene(readfile)
outputFile <- paste(outputFolder_path,"downRegulated_genes_Annotation.csv",sep="/")
write.table(annotation_df,outputFile, sep="\t",row.names = FALSE)


