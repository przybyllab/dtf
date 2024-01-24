library(matrixStats)
library(dplyr)
library(clustvis)

#Creating a function to calculate the read count variation 
CalVar <- function(df) {
  
  #converting all read counts to numeric type as they are read as a character when read from a csv
  #5 is hard coded, this needs to be changed f the description data is changed
  valuesAsNumeric = mutate_all(df[5:dim(df)[1],], function(x) as.numeric(as.character(x)))
  variance = rowVars(as.matrix(valuesAsNumeric)) #calculating the variation 
  
  #adding the read counts variation
  #the read count variance vector will be added as a column to the dataframe data_of_interest,
  #therefore the vector has to have the same length as there are number of rows, so we have to add 
  #null values so that the vector
  #-1s will be added because variance can not be a -1 value
  rowOfData_of_interest =  dim(df)[1]
  NumberOfNullValuesToAdd = rowOfData_of_interest - length(variance)
  variance_column = c( rep(-1, each= NumberOfNullValuesToAdd), variance) # concatenating 26 -1s 
  df$variance = variance_column
  return(df)
}

#Creating a function to take the top 1% variation read counts 
top1Percent_Var <- function(df){
  
  #note that because the first rows of our dataframe hold description data, they will have a -1 read count variation value
  variance_column = df$variance #getting a vector of read count variation
  Is_neg1 = variance_column == -1 
  descriptionRowsCutOff = match(FALSE, Is_neg1 ) - 1 
  
  #separating the descriptor rows from RNAseq vsd rows 
  designData = df[1:descriptionRowsCutOff,]
  last_row = dim(df)[1]
  RNAseq = df[descriptionRowsCutOff+1 : last_row,]
  
  #sorting by variance
  sortedVarValues <- RNAseq [order(RNAseq$variance, decreasing = TRUE),]
  
  #getting top 1%
  index_of_top1percent = c(1:(dim(sortedVarValues)[1]*.01))
  top_1percent_var = sortedVarValues[index_of_top1percent,]
  
  #concatenating the design file and the top 1% variable read counts
  descript_top1percent = rbind(designData, top_1percent_var)
  
  last_column = dim(descript_top1percent)[2] #this is the beta Variance row
  descript_top1percent = descript_top1percent[,1:last_column -1] #removing the beta variance row
  
  return(descript_top1percent)
}

for (Patient_ID in c(1,5,8)){

  #reading in the file and processing it so that it is the proper format
  tdata <- read.csv(sprintf('[VSD value Pathway]\\Pt%d_vsd.csv', Patient_ID))
  
  rownames(tdata) <- make.names(tdata[,1], unique = TRUE)
  tdata <- tdata[,2: dim(tdata)[2]] #removing the first column
  designData <- read.csv(sprintf('[design data file pathway]\\design_pt%d.csv', Patient_ID), row.names = 1)
  print(dim(designData))
  designData <- t(designData)
  
  clustvis_format = rbind(designData, tdata)
  clustvis_format = clustvis_format[c(1:dim(clustvis_format)[1]),]
  vsd_with_variance <- CalVar(clustvis_format)
  top1PerVar_vsd <- top1Percent_Var(vsd_with_variance )
  
  #writing dataframe to csv
  
  write.csv(top1PerVar_vsd, 
            sprintf('[Path to where to save clustvis formatted RNA-seq data]\\RNA_seq_design_info_with_vsd_values_Top1Percent_Patient%d_4April2023.csv', Patient_ID))
}

pathways = c('PROGNOSTIC_MARKERS',)

for (Patient_ID in c(1,5,8)){
    design_file = sprintf("[design data file pathway]\\design_pt%d.csv", Patient_ID)
    file_path = sprintf("[VSD value Pathway]\\Pt%d_vsd_", Patient_ID)
    for(pathway in pathways){
    
      #reading in the file and processing it so that it is the proper format
      file_path = sprintf("[Path to RNA-seq vsd values]\\Pt%d_vsd_", Patient_ID)
      file_with_edata <- paste(file_path, pathway, ".csv", sep = "")
      edata <-read.csv(file_with_edata)
      rownames(edata) <- make.names(edata[,1], unique = TRUE) #issue when setting rowname using csv
      edata <- edata[,2: dim(edata)[2]] #removing the first column
      designData <- read.csv(design_file, row.names = 1)
      
      designData <- t(designData)
      
      clustvis_format = rbind(designData, edata)
      clustvis_format = clustvis_format[c(1:dim(clustvis_format)[1]),]
      
      
      clustvis_file_name = sprintf("design_info_Pt%d_vsd_%s_23March2023", Patient_ID, pathway)
      write.csv(clustvis_format, 
                sprintf('[Path to where to save clustvis formatted RNA-seq data]\\%s.csv', clustvis_file_name))
      
    
  }
}

