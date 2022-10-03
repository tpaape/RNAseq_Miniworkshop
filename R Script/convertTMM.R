rawCounts_toTMM <- function(counts, metadata, treatment_Column) {
  library(edgeR)
  library(stringr)
  #Import count data. Assumes the first column contains the gene identifier
  if (grepl('.csv', counts, fixed=TRUE) == TRUE) {
  countData <- read.csv(counts, row.names=1)
  } else if (grepl('.tsv', counts, fixed=TRUE) == TRUE) {
    countData <- read.csv(counts, row.names=1, sep='\t')
  } else if (grepl('.txt.', counts, fixed=TRUE) == TRUE) {
    print('.txt file provided as counts input. Assumming it is tab delimited.')
    countData <- read.csv(counts, row.names=1, sep='\t')
  }
  
  #Import metadata. Assumes the first column contains the sample identifier
  if (grepl('.csv', metadata, fixed=TRUE) == TRUE) {
    metaData <- read.csv(metadata)
  } else if (grepl('.tsv', metadata, fixed=TRUE) == TRUE) {
    metaData <- read.csv(metadata, sep='\t')
  } else if (grepl('.txt.', metadata, fixed=TRUE) == TRUE) {
    print('.txt file provided as metadat input. Assumming it is tab delimited.')
    metaData <- read.csv(metadata, sep='\t')
  }
  
  #metaData$Name <- paste('X', metaData$Name, sep='') #an X may be added to sample names by default if they start with a number. Uncomment this line if you need an X added to metaData sample names to match count file. 
  
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

