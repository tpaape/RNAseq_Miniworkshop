rawCounts_toTMM <- function(counts, metadata, treatment_Column) {
  library(edgeR)
  library(stringr)
  #Import count data. Assumes the first column contains the gene identifier
  countData <- read.csv(counts, row.names=1, sep='\t')
  
  #Import metadata. Assumes the first column contains the sample identifier
  metaData <- read.csv(metadata, sep='\t')
  
  metaData$Name <- paste('X', metaData$Name, sep='')
  
  for_groups <- data.frame(colnames(countData))
  metaData$Name <- str_replace(metaData$Name, '-', '.')
  for_groups <- merge(for_groups, metaData, by.x='colnames.countData.', by.y='Name', all.x=TRUE)
  
  col <- grep(treatment_Column, colnames(for_groups))

  
  for_groups <- for_groups[match(colnames(countData), for_groups$colnames.countData.),]

  groups <- factor(for_groups[[col]])

  #Removes categorical data from count dataframe
  noCat <- sapply(countData, is.numeric) #comment this line out if removing categorical variables is messing up your dataframe
  countData <- countData[noCat] #comment this line out if removing categorical variables is messing up your dataframe

  y <- DGEList(counts=countData, group=groups)

  #Filters out genes with low counts. This is typically done before TMM normalization
  keep <- filterByExpr(y) #comment out this line if you do not want to remove genes with low counts
  y <- y[keep,,keep.lib.sizes=FALSE] #comment out this line if you do not want to remove genes with low counts
  dgeList <- calcNormFactors(y, method="TMM")
  tmm <- cpm(dgeList)
  
  newName <- sub(".csv", '_TMM.csv',counts)
  
  write.csv(as.data.frame(tmm), 
            file=newName)
}


#Example usage  
rawCounts_toTMM(counts='13320077-counts.txt', 
                metadata='test_metadata.txt', 
                treatment_Column="Treatment")

rawCounts_toTMM(counts="C:/Users/mclea/OneDrive/Brookhaven/RNAseq/Time_Series/Poplar/Raw Output/20400Tim-raw_genes_counts_INT_both_outliers_removed.csv", 
                metadata="C:/Users/mclea/OneDrive/Brookhaven/RNAseq/Time_Series/All_metadata_pop_noOutliers.csv", 
                treatment_Column="TissueXTrtXTime")