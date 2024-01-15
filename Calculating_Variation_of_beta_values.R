library(matrixStats)
library(dplyr)

#Creating a function to calculate the beta value variation 
CalBetaVal <- function(df) {
  
  #converting all beta values to numeric type as they are read as a character when read from a csv
  #5 is hard coded, this needs to be changed if the description data is changed
  betavaluesAsNumeric = mutate_all(df[5:518437,], function(x) as.numeric(as.character(x)))
  beta_variance = rowVars(as.matrix(betavaluesAsNumeric)) #calculating the variation 
  
  #adding the beta values variation
  #the length of the beta_variance vector is 518 433 and will be adding it as a column to the dataframe data_of_interest,
  #therefore the vector has to have the same length are there are number of rows, so we have to add 
  #null values so that the vector has the correct length
  #I will be adding -1s because variance can not be a -1 value and we are selecting for the largest variance values
  rowOfData_of_interest =  dim(df)[1]
  NumberOfNullValuesToAdd = rowOfData_of_interest - length(beta_variance)
  beta_variance_column = c( rep(-1, each= NumberOfNullValuesToAdd), beta_variance) # concatenating 26 -1s 
  df$beta_variance = beta_variance_column
  return(df)
}


#Creating a function to take the top 1000 variation beta values 
top1000_BetaVar <- function(df){
  
  #note that because the first rows of our dataframe hold description data, they will have a -1 beta_variation value
  beta_variance_column = df$beta_variance #getting a vector of beta variation
  Is_neg1 = beta_variance_column == -1 
  descriptionRowsCutOff = match(FALSE, Is_neg1 ) - 1 
  
  #separating the descriptor rows from cg rows which hold the beta values
  designData = df[1:descriptionRowsCutOff,]
  last_row = dim(df)[1]
  CgRegions = df[descriptionRowsCutOff+1 : last_row,]
  
  #sorting by variance
  sortedBetaVarValues <- CgRegions[order(CgRegions$beta_variance, decreasing = TRUE),]
  
  #getting top 1000
  top_1000_var = sortedBetaVarValues[c(1:1000),]
  
  #concatenating the design file and the top 1000 var cg regions
  descript_top1000 = rbind(designData, top_1000_var)
  
  last_column = dim(descript_top1000)[2] #this is the beta Variance row
  descript_top1000 = descript_top1000[,1:last_column -1] #removing the beta variance row
  
  return(descript_top1000)
}

#Creating a function to take the top 1000 variation beta values 
topx_BetaVar <- function(df, x){
  
  #note that because the first rows of our dataframe hold description data, they will have a -1 beta_variation value
  beta_variance_column = df$beta_variance #getting a vector of beta variation
  Is_neg1 = beta_variance_column == -1 
  descriptionRowsCutOff = match(FALSE, Is_neg1 ) - 1 
  
  #separating the descriptor rows from cg rows which hold the beta values
  designData = df[1:descriptionRowsCutOff,]
  last_row = dim(df)[1]
  CgRegions = df[descriptionRowsCutOff+1 : last_row,]
  
  #sorting by variance
  sortedBetaVarValues <- CgRegions[order(CgRegions$beta_variance, decreasing = TRUE),]
  
  #getting top x
  top_x_var = sortedBetaVarValues[c(1:x),]
  
  #concatenating the design file and the top x var cg regions
  descript_topx = rbind(designData, top_x_var)
  
  last_column = dim(descript_topx)[2] #this is the beta Variance row
  descript_topx = descript_topx[,1:last_column -1] #removing the beta variance row
  
  return(descript_topx)
}

#Creating a function to take the top 1% variation beta values 
top1Percent_BetaVar <- function(df){
  
  #note that because the first rows of our dataframe hold description data, they will have a -1 beta_variation value
  beta_variance_column = df$beta_variance #getting a vector of beta variation
  Is_neg1 = beta_variance_column == -1 
  descriptionRowsCutOff = match(FALSE, Is_neg1 ) - 1 
  
  #separating the descriptor rows from cg rows which hold the beta values
  designData = df[1:descriptionRowsCutOff,]
  last_row = dim(df)[1]
  CgRegions = df[descriptionRowsCutOff+1 : last_row,]
  
  #sorting by variance
  sortedBetaVarValues <- CgRegions[order(CgRegions$beta_variance, decreasing = TRUE),]
  
  #getting top 1%
  index_of_top1percent = c(1:(dim(sortedBetaVarValues)[1]*.01))
  top_1percent_var = sortedBetaVarValues[index_of_top1percent,]
  
  #concatenating the design file and the top 1000 var cg regions
  descript_top1percent = rbind(designData, top_1percent_var)
  
  last_column = dim(descript_top1percent)[2] #this is the beta Variance row
  descript_top1percent = descript_top1percent[,1:last_column -1] #removing the beta variance row
  
  return(descript_top1percent)
}


####

#getting top 10%
#top_10percent_rows <- c(1:(dim(sortedData)[1]*.1))
#top_10percent_var = sortedData[top_10percent_rows,]

###


#Read in the data
edata <- read.csv('C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\20samples_detailed_and_design_info_DTFhet_with_rounded_beta_matrix18Jan2023.csv', row.names = 1)
dim(edata)
head(edata)

#selecting the data 
data_of_interest <- edata[c(6,7, 11,15,27:518459),] #selecting the patient number, type, mutation, year, and beta values

#separating the data by patient
data_of_interest_pt1 <- data_of_interest %>% select_if(function(col) col[1]== "pt1") # selecting patient 1
data_of_interest_pt5 <- data_of_interest %>% select_if(function(col) col[1]== "pt5") # selecting patient 5
data_of_interest_pt8 <- data_of_interest %>% select_if(function(col) col[1]== "pt8") # selecting patient 8

#removing the row that was used for selecting the data 
data_of_interest_pt1 <- data_of_interest_pt1[2: dim(data_of_interest_pt1)[1],]
data_of_interest_pt5 <- data_of_interest_pt5[2: dim(data_of_interest_pt5)[1],]
data_of_interest_pt8 <- data_of_interest_pt8[2: dim(data_of_interest_pt8)[1],]

#calculating the variation for patient 1
data_of_interest_pt1 <- CalBetaVal(data_of_interest_pt1)

write.csv(top1000_BetaVar(data_of_interest_pt1), "C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\top_1000_beta_value_var_20DTFsamples_Patient_1.csv") #writing to file

write.csv(top1Percent_BetaVar(data_of_interest_pt1), "C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\top_1percent_beta_value_var_20DTFsamples_Patient_1.csv") #writing to file

write.csv(topx_BetaVar(data_of_interest_pt1, 100),"C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\top_100_beta_value_var_20DTFsamples_Patient_1.csv")
write.csv(topx_BetaVar(data_of_interest_pt1, 200),"C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\top_200_beta_value_var_20DTFsamples_Patient_1.csv")

#calculating the variation for patient 5
data_of_interest_pt5 <- CalBetaVal(data_of_interest_pt5)

write.csv(top1000_BetaVar(data_of_interest_pt5), "C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\top_1000_beta_value_var_20DTFsamples_Patient_5.csv") #writing to file

write.csv(top1Percent_BetaVar(data_of_interest_pt5), "C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\top_1percent_beta_value_var_20DTFsamples_Patient_5.csv") #writing to file

#calculating the variation for patient 8
data_of_interest_pt8 <- CalBetaVal(data_of_interest_pt8)

write.csv(top1000_BetaVar(data_of_interest_pt8), "C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\top_1000_beta_value_var_20DTFsamples_Patient_8.csv") #writing to file

write.csv(top1Percent_BetaVar(data_of_interest_pt8), "C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Methylation\\top_1percent_beta_value_var_20DTFsamples_Patient_8.csv") #writing to file
