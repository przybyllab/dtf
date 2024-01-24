library(matrixStats)
library(clustvis)
library(dplyr)


#getting background information about the parameters the code was ran under
Package_version = packageVersion("clustvis")
Date = as.Date(Sys.Date(), "%Y_%m_%d")

#choosing the scalling that will be used for analysis
for (scalling_type in c("none")){
  
  #choosing the different patients that will be analyzed
  for (Patient_ID in c(1,5,8)){ 
    
    #choosing the different type of grouping algorithms to use 
    for (clsuterType in c("complete","average","mcquitty", "ward.D")){
      
      #choosing the datasets that are used for this run
      for (topxvar in  c("top_1percent")){
        
        #reading files 
        file = sprintf("[pathwaywat to clustvis file]\\%s_beta_value_var_20DTFsamples_Patient_%d.csv", topxvar, Patient_ID)
        
        #preparing data for clustvis 
        imp = importData(file)
        proc = processData(imp, rowScaling = scalling_type, rowCentering = F)
        
        #pca generation
        pca = generatePCA(proc, pcx = 1, pcy = 2, switchDirX = FALSE,
                          switchDirY = FALSE, colorAnno = "Sample", colorScheme = "Dark2",
                          showEllipses = FALSE, ellipseConf = 0.95, ellipseLineWidth = 1,
                          ellipseLineType = "solid", shapeAnno = NULL, shapeScheme = "various",
                          plotWidth = 20, plotRatio = 0.8, marginRatio = 0.05, pointSize = 5,
                          legendPosition = "right", fontSize = 20, axisLabelPrefix = "PC",
                          showVariance = TRUE, showSampleIds = FALSE, maxColorLevels = 30,
                          maxShapeLevels = 62)
        savePCA(pca, file = sprintf("[Pathway to save PCA]\\clustvis_PCA_Euclidean_Ward_20DTFsamples_NOrowCentering_%s_scalling_%s_beta_matrix_Patient%d_%s_clustvisVersion_%s.pdf",Patient_ID,scalling_type, topxvar,Patient_ID,Date, Package_version))
        
        #heatmap generation
        hm = generateHeatmap(proc, clustDistRows = NA, clustMethodRows = NA, clustDistCols = "euclidean", 
                             clustMethodCols = clsuterType,fontSizeColnames = 20, fontSizeGeneral = 15, legendColorScheme =  "Dark2", colorAnnoCol = c("Type", "Mutation","Year"))
        saveHeatmap(hm, file = sprintf("[Pathway to heatmaps]\\clustvisHeatmap_Euclidean_%s_20DTFsamples_NOrowcentering_%s_scalling_%s_beta_matrix_Patient%d_%s_clustvisVersion_%s.pdf", Patient_ID, clsuterType,scalling_type, topxvar,Patient_ID,Date, Package_version))
        
      }
    }
  }
}




