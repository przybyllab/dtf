library(matrixStats)
library(dplyr)
library(clustvis)



#Creating a function to calculate the beta value variation 
CalVar <- function(df) {
  
  #converting all beta values to numeric type as they are read as a character when read from a csv
  #3 is hard coded, this needs to be changed f the description data is changed
  valuesAsNumeric = mutate_all(df[3:dim(df)[1],], function(x) as.numeric(as.character(x)))
  variance = rowVars(as.matrix(valuesAsNumeric)) #calculating the variation 
  
  #adding the beta values variation
  #the beta_variance vector will be adding as a column to the dataframe data_of_interest,
  #therefore the vector has to have the same length as there are number of rows, so we have to add 
  #null values so that the vector has the correct length
  #-1s  will be added because variance can not be a -1 value
  rowOfData_of_interest =  dim(df)[1]
  NumberOfNullValuesToAdd = rowOfData_of_interest - length(variance)
  variance_column = c( rep(-1, each= NumberOfNullValuesToAdd), variance) # concatenating 26 -1s 
  df$variance = variance_column
  return(df)
}


#Creating a function to take the top 1% variation beta values 
top1Percent_Var <- function(df){
  
  #note that because the first rows of our dataframe hold description data, they will have a -1 beta_variation value
  variance_column = df$variance #getting a vector of beta variation
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
  
  #concatenating the design file and the top 1% variable beta values
  descript_top1percent = rbind(designData, top_1percent_var)
  
  last_column = dim(descript_top1percent)[2] #this is the beta Variance row
  descript_top1percent = descript_top1percent[,1:last_column -1] #removing the beta variance row
  
  return(descript_top1percent)
}


#_____________________________________________Creating Clustvis df____________________________________#

#reading in the file and processing it so that it is the proper format
tdata <- read.csv('C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\Pt1_vsd.csv')
rownames(tdata) <- make.names(tdata[,1], unique = TRUE) #issue when setting rowname using csv
tdata <- tdata[,2: dim(tdata)[2]] #removing the first column
designData <- read.csv('C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\design_pt1.csv', row.names = 1)

designData <- t(designData)

clustvis_format = rbind(designData, tdata)
clustvis_format = clustvis_format [c(1,c(3, c(4:dim(clustvis_format)[1]))),]


#_______________________________________________selecting data of interest___________________________________#


vsd_with_variance <- CalVar(clustvis_format)
top5PerVar_vsd <- top5Percent_Var(vsd_with_variance)
top1PerVar_vsd <- top1Percent_Var(vsd_with_variance )

#writing dataframe to csv
write.csv(top5PerVar_vsd, 
          'C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\RNA_seq_design_info_with_vsd_values_Top5Percent_Patient1_9Feb2023.csv')

write.csv(top1PerVar_vsd, 
          'C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\RNA_seq_design_info_with_vsd_values_Top1Percent_Patient1_9Feb2023.csv')





#_____________________________________________Creating Clustvis df____________________________________#

#reading in the file and processing it so that it is the proper format
tdata <- read.csv('C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\Pt5_vsd.csv')
rownames(tdata) <- make.names(tdata[,1], unique = TRUE) #issue when setting rowname using csv
tdata <- tdata[,2: dim(tdata)[2]] #removing the first column
designData <- read.csv('C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\design_pt5.csv', row.names = 1)

designData <- t(designData)

clustvis_format = rbind(designData, tdata)
clustvis_format = clustvis_format [c(1,c(3, c(4:dim(clustvis_format)[1]))),]


#_______________________________________________selecting data of interest___________________________________#


vsd_with_variance <- CalVar(clustvis_format)
top5PerVar_vsd <- top5Percent_Var(vsd_with_variance)
top1PerVar_vsd <- top1Percent_Var(vsd_with_variance )

#writing dataframe to csv
write.csv(top5PerVar_vsd, 
          'C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\RNA_seq_design_info_with_vsd_values_Top5Percent_Patient5_9Feb2023.csv')

write.csv(top1PerVar_vsd, 
          'C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\RNA_seq_design_info_with_vsd_values_Top1Percent_Patient5_9Feb2023.csv')




#_____________________________________________Creating Clustvis df____________________________________#

#reading in the file and processing it so that it is the proper format
tdata <- read.csv('C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\Pt8_vsd.csv')
rownames(tdata) <- make.names(tdata[,1], unique = TRUE) #issue when setting rowname using csv
tdata <- tdata[,2: dim(tdata)[2]] #removing the first column
designData <- read.csv('C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\design_pt8.csv', row.names = 1)

designData <- t(designData)

clustvis_format = rbind(designData, tdata)
clustvis_format = clustvis_format [c(1,c(3, c(4:dim(clustvis_format)[1]))),]


#_______________________________________________selecting data of interest___________________________________#


vsd_with_variance <- CalVar(clustvis_format)
top5PerVar_vsd <- top5Percent_Var(vsd_with_variance)
top1PerVar_vsd <- top1Percent_Var(vsd_with_variance )

#writing dataframe to csv
write.csv(top5PerVar_vsd, 
          'C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\RNA_seq_design_info_with_vsd_values_Top5Percent_Patient8_9Feb2023.csv')

write.csv(top1PerVar_vsd, 
          'C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\RNA_seq_design_info_with_vsd_values_Top1Percent_Patient8_9Feb2023.csv')

