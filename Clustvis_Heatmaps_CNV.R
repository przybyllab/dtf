library(matrixStats)
library(dplyr)
library(clustvis)

#getting background information about the parameters the code was ran under
Package_version = packageVersion("clustvis")
Date = as.Date(Sys.Date(), "%Y_%m_%d")


#_______________________________________________creating pca and heatmap___________________________________#
#choosing the different type of clustering algorithms to use 
for (clusterType in c("complete","average","mcquitty", "ward.D")){ 
  
  #choosing the datasets that are used for this run
  for (merge_dvg in c( "_nomerge_dvgPartialNone")){
    
    #reading files for patient 5
    file = sprintf("[Insert path to cnv clustvis formatted file]\\clustvisPatient5_CNVs_filtered%s.csv", merge_dvg)
    
    #preparing data for clustvis 
    imp = importData(file, nbrRowAnnos = 0, nbrColAnnos = 5)
    proc = processData(imp, rowScaling = "none", rowCentering = FALSE)
    
    #heatmap generation
    hm = generateHeatmap(proc, clustDistRows = NA, clustMethodRows = "ward.D", clustDistCols = "euclidean", 
                         clustMethodCols = clusterType, showRownames= F, colorAnnoCol = c("Type", "Mutation","Year"), legendColorScheme =  "Set2", colorRangeMax = 0.5, colorRangeMin = -0.5)
    saveHeatmap(hm, file = sprintf("[Path to here to save cnv heatmaps]\\clustvisHeatmap_CNV_%s_20DTF_NOrowCenter_NOscalling_%s_Patient5_%s_clustvis_%s.pdf",clsuterType, merge_dvg,Date, Package_version))
    
    #reading files for patient 8 
    file = sprintf("[Insert pathway to cnv file]\\clustvisPatient8_CNVs_filtered%s.csv", merge_dvg)
    imp = importData(file, nbrRowAnnos = 0, nbrColAnnos = 5)
    proc = processData(imp, rowScaling = "none", rowCentering = FALSE)
    
    #heatmap generation
    hm = generateHeatmap(proc, clustDistRows = NA, clustMethodRows = "ward.D", clustDistCols = "euclidean", 
                         clustMethodCols = clusterType, showRownames= F, colorAnnoCol = c("Type", "Mutation","Year"),legendColorScheme =  "Set2",colorRangeMax = 0.5, colorRangeMin = -0.5)
    saveHeatmap(hm, file = sprintf("[Path to here to save cnv heatmaps]\\clustvisHeatmap_CNV_%s_20DTF_NOrowCenter_NOscalling_%s_Patient8_%s_clustvis_%s.pdf",clsuterType, merge_dvg,Date, Package_version))
  }
}


#patient 1 did not have enough data after dvg filtering for us to analyze using 
#clustvis therefore only the unfiltered data is being used for analysis 

#choosing the different type of grouping algorithms to use 
for (clsuterType in c("complete","average","mcquitty", "ward.D")){
  
  #choosing the datasets that are used for this run
  for (merge_dvg in c("_nomerge_dvgNotConsidered")){
  
    #reading files for patient 1
    file = sprintf("[Insert path to cnv clustvis formatted file]\\clustvisPatient1_CNVs_filtered%s.csv", merge_dvg)
    
    #preparing data for clustvis 
    imp = importData(file, nbrRowAnnos = 0, nbrColAnnos = 5)
    proc = processData(imp, rowScaling = "none", rowCentering = FALSE)
    
    #heatmap generation
    hm = generateHeatmap(proc, clustDistRows = NA, clustMethodRows = "ward.D", clustDistCols = "euclidean", 
                         clustMethodCols = clsuterType,showRownames= F, colorAnnoCol = c("Type", "Mutation","Year"), legendColorScheme =  "Set2", colorRangeMax = 0.8, colorRangeMin = -0.8)
    saveHeatmap(hm, file = sprintf("[Path to here to save cnv heatmaps]\\clustvisHeatmap_CNV_%s_20DTF_NOrowCenter_NOscalling_%s_Patient1_%s_clustvis_%s.pdf",clsuterType, merge_dvg,Date, Package_version))
  }
}





