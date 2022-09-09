rawCounts_toTMM <- function(path_to_directory) {
  library(edgeR)
  files <- list.files(path=path_to_directory, pattern="*.csv", full.names=TRUE,
                      recursive=FALSE)
  
  for (i in files){
    #Import count data. Assumes the first column contains the gene identifier
    countData <- read.csv(i, sep=",", row.names=1)
    
    #Removes categorical data from count dataframe
    noCat <- sapply(countData, is.numeric) #comment this line out if removing categorical variables is messing up your dataframe
    countData <- countData[noCat] #comment this line out if removing categorical variables is messing up your dataframe
    
    y <- DGEList(counts=countDataCleaned)
    
    #Filters out genes with low counts. This is typically done before TMM normalization
    keep <- filterByExpr(y) #comment out this line if you do not want to remove genes with low counts
    y <- y[keep,,keep.lib.sizes=FALSE] #comment out this line if you do not want to remove genes with low counts
    dgeList <- calcNormFactors(y, method="TMM")
    tmm <- cpm(dgeList)
    
    newName <- sub(".csv", '_TMM.csv',i)
    
        write.csv(as.data.frame(tmm), 
                  file=newName)  
    
  }

}

#Example usage
rawCounts_toTMM("direct/path/to/your/count/directory")



