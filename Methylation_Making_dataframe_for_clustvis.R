library(matrixStats)


#reading in experimental data 
edata <- read.csv("[Pathway to beta values]\\20DTFsamples_rounded_beta_matrix_28Oct2022.csv', row.names = 1)

#reading design file 
ddata <- read.csv('[design file for methylation data]\\design_DTF_ALL_het_18Jan2023_design_DTF67_EE.csv')

#read in detailed methylation file 
dmdata <- read.csv('C[desing file for methylation data]\\20samples_detailed_info_DTFhet_18Jan2023.csv')


#turning STT into the index
rownames(ddata) <- ddata$STT #making STTs the index
rownames(dmdata) <- dmdata$STT #making STTs the index


#cleaning data for merge
ddata <- ddata[ddata$Methylation == 'YES',] #selecting those that were studied for methylation
ddata <- ddata[,c(2:17)] #dropping STT

design_methylation_detail_data <- cbind(ddata, dmdata)



#renaming sample name to have the same format as the methylation data
renameSample <- function(x) {
  # replace "-" with "." using gsub()
  gsub("-", ".", x) 
}
design_methylation_detail_data$Sample_name  <- lapply(design_methylation_detail_data$Sample_name, renameSample)

rownames(design_methylation_detail_data) <- design_methylation_detail_data$Sample_name #making sample_names the index

design_methylation_detail_data <- t(design_methylation_detail_data) #transpose of df to make the sample name the column names



#merging desing info with experimental
design_methylation__data <- rbind(design_methylation_detail_data, edata)

#unlist all columns as they were all lists
design_methylation__data_no_lists <-data.frame(lapply(design_methylation__data, function(x) unlist(x)))
rownames(design_methylation__data_no_lists) <- rownames(design_methylation__data)

#writing dataframe to csv
write.csv(design_methylation__data_no_lists, 
          '[merged experimental and desing info pathway]\\20samples_detailed_and_design_info_DTFhet_with_rounded_beta_matrix18Jan2023.csv')

