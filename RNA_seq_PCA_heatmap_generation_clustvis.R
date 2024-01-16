library(matrixStats)
library(clustvis)
library(dplyr)


#getting background information about the parameters the code was ran under
Patient_ID = 5
Package_version = packageVersion("clustvis")
Date = as.Date(Sys.Date(), "%Y_%m_%d")

#creating a vector of the pathways under analysis 
pathways = c('PROGNOSTIC_MARKERS',)

#choosing the scalling that will be used for analysis
for (scalling_type in c("uv")){
  
  #choosing the different type of clustering algorithms to use 
  for (clusterType in c("complete","average","mcquitty", "ward.D")){
    
    #choosing the different patients that will be analyzed
    for (Patient_ID in c(1,5,8)){
      
      #choosing the pathways that will be analyzed
      for (pathway in pathways){
        
        #reading file 
        file=sprintf("C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\Clustvis_format\\design_info_Pt%d_vsd_%s_23March2023.csv", Patient_ID, pathway)
        
        #preparing data for clustvis 
        imp = importData(file)
        proc = processData(imp, rowScaling = scalling_type, rowCentering = TRUE)
        
        #pca generation
        pca = generatePCA(proc, pcx = 1, pcy = 2, switchDirX = FALSE,
                          switchDirY = FALSE, colorAnno = "Sample", colorScheme = "Set2",
                          showEllipses = FALSE, ellipseConf = 0.95, ellipseLineWidth = 1,
                          ellipseLineType = "solid", shapeAnno = NULL, shapeScheme = "various",
                          plotWidth = 20, plotRatio = 0.8, marginRatio = 0.05, pointSize = 8,
                          legendPosition = "right", fontSize = 20, axisLabelPrefix = "PC",
                          showVariance = TRUE, showSampleIds = FALSE, maxColorLevels = 20,
                          maxShapeLevels = 62)
        savePCA(pca, file = sprintf("C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Exploratory_Analysis\\Patient%d\\RNA_seq\\PCA\\clustvis_PCA_EuclidWard_20DTF_rowCenter_%s_scalling_%s_vsd_Patient%d_%s_clustvis_%s.pdf",Patient_ID,scalling_type,pathway,Patient_ID,Date, Package_version))
        
        
        #heatmap generation
        hm = generateHeatmap(proc, clustDistRows = NA, clustMethodRows = NA, clustDistCols = "euclidean", 
                             clustMethodCols = clusterType, legendColorScheme =  "Set2", fontSizeColnames = 20, 
                             fontSizeGeneral = 10, colorAnnoCol = c("Type", "Mutation","Year"), plotWidth = 15)
        saveHeatmap(hm, file = sprintf("C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Exploratory_Analysis\\Patient%d\\RNA_seq\\Heatmap\\clustvisHeatmap_Euclidean_%s_20DTF_rowcenter_none_scalling_%s_vsd_Patient%d_%s_clustvis_%s.pdf", Patient_ID,  clusterType, pathway,Patient_ID,Date, Package_version))
      }
    }
  }
  
  #choosing the different type of clustering algorithms to use 
  for (clusterType in c("complete","average","mcquitty", "ward.D")){
    
    #choosing the different patients that will be analyzed
    for (Patient_ID in c(1,5,8)){
      
      #choosing the datasets that are used for this run
      for (topxvar in  c("Top1Percent",)){
        
        #reading file
        file = sprintf("C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\RNA_seq\\Clustvis_format\\RNA_seq_design_info_with_vsd_values_%s_Patient%d_4April2023.csv", topxvar, Patient_ID)
        
        #preparing data for clustvis
        imp = importData(file)
        proc = processData(imp, rowScaling = scalling_type, rowCentering = TRUE)
        
        #pca generation
        pca = generatePCA(proc, pcx = 1, pcy = 2, switchDirX = FALSE,
                          switchDirY = FALSE, colorAnno = "Sample", colorScheme = "Set2",
                          showEllipses = FALSE, ellipseConf = 0.95, ellipseLineWidth = 1,
                          ellipseLineType = "solid", shapeAnno = "Type", shapeScheme = "various",
                          plotWidth = 20, plotRatio = 0.8, marginRatio = 0.05, pointSize = 8,
                          legendPosition = "right", fontSize = 26, axisLabelPrefix = "PC",
                          showVariance = TRUE, showSampleIds = FALSE, maxColorLevels = 15,
                          maxShapeLevels = 62)
        savePCA(pca, file = sprintf("C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Exploratory_Analysis\\Patient%d\\RNA_seq\\PCA\\clustvis_PCA_Euclidean_Ward_20DTFsamples_rowCenter_%s_scalling_%s_vsd_Patient%d_%s_clustvisVersion_%s.pdf",Patient_ID,scalling_type, topxvar,Patient_ID,Date, Package_version))
        
        #heatmap generation
        hm = generateHeatmap(proc, clustDistRows = NA, clustMethodRows = NA, clustDistCols = "euclidean",
                             colorAnnoCol = c("Type", "Mutation","Year"),
                             clustMethodCols = clusterType, legendColorScheme =  "Set2", fontSizeColnames = 20,  fontSizeGeneral = 15)
        saveHeatmap(hm, file = sprintf("C:\\Users\\chels\\OneDrive\\Documents\\Pryzbyl_Lab\\Exploratory_Analysis\\Patient%d\\RNA_seq\\Heatmap\\clustvisHeatmap_Euclidean_%s_20DTFsamples_rowcenter_none_scalling_%s_vsd_Patient%d_%s_clustvisVersion_%s.pdf", Patient_ID, clusterType, topxvar,Patient_ID,Date, Package_version))
      }
    }
  }
}


