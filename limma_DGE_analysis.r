# author : Sandali Lokuge
# date : 12/11/2022
# differential gene expression analysis with limma
# set the working directory to main (project) file location

library(limma)
library(pheatmap)


# set the working directory as the main project directory - setwd()

################################ functions ################################

limmaDGE <- function(data,typeCount){
  
  # read the normalized exp data
  df <- data[,2:ncol(data)]
  ex<-new("ExpressionSet",exprs=as.matrix(df))
  edata<-exprs(ex)
  
  # produce a design matrix (also known as model matrices) for a variety of linear models from limma package 
  conditons <- rep(c("tumor","nontumor"),times=c((typeCount[1]+typeCount[3]+typeCount[5]),(typeCount[2]+typeCount[4]+typeCount[6])))
  colnames(edata) <- conditons
  
  # design matrix of the microarray experiment
  design <- model.matrix(~0+factor(conditons))
  colnames(design) <- c("tumor","nontumor")
  
  y <- voom(edata, design, plot = T)
  
  # lmFit is a function from the limma package
  # It has two main arguments: the expression data and the design matrix. 
  # lmFit fits linear model for each gene given a series of arrays.
  
  fit <- lmFit(y, design)
  coef.fit <- fit$coefficients
  head(coef(fit))
  #gene expression in tumor compared to nontumor
  contr <- makeContrasts(tumor - nontumor, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  
  res_df <- as.data.frame(top.table)
  res_df <- add_column(res_df, data['Probe_ID'], .before = 1)
  
  return(res_df)
  
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

root_path <- getwd()
folder_path <- paste(root_path,"Expression_data",cancerType,sep="/")

data <- read.csv(paste(folder_path,"Ovarian cancer_afterNorm.csv",sep="/"),sep='\t')

limma_output <- limmaDGE(data,typeCount)

# save your results
outputFolder <- paste(root_path,"Limma_outputs",cancerType,sep="/")
write.table(limma_output,paste(outputFolder,"Ovarian cancer_limmaOutput.csv",sep="/"),sep="\t",row.names = FALSE)


