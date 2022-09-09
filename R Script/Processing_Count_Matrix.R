library(DESeq2)
library(pheatmap)

#Processing Count Data RNAseq Miniworkshop
#setwd("~/full_data/deseq2")

#Load the count matrix data download from the instance
counts <- read.delim('13320077-counts.txt', header=T, sep='\t', row.names=1)

#Load metadata file
metaData <- read.delim('test_metadata.txt', header=T, sep='\t')

# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(metaData) <- metaData$Samples

#Create a DESeq2 object from the matrix file
dds <- DESeqDataSetFromMatrix(countData = counts,
                                         colData = metaData,
                                         design = ~ Treatment)

# Number of genes before filtering:
nrow(dds)

# Filter
dds <- dds[rowSums(counts(dds)) > 10, ]


# Number of genes left after low-count filtering:
nrow(dds)

#Fit model to the data
dds_model <- DESeq(dds)

#Calculates log2 normalized counts
norm_counts_log2 <- log2(counts(dds_model, normalized = TRUE)+1)
write.table(norm_counts_log2, "normalized_counts_log2.txt", quote=F, col.names=T, row.names=T, sep="\t")

#Calculates TPMs
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
################################################################
#Visualize the counts with between-sample distance matrix heatmap
################################################################

# Try with the vst transformation
vst_norm <- vst(dds_model)

# calculate between-sample distance matrix
sampleDistMatrix <- as.matrix(dist(t(assay(vst_norm))))

# create figure in PNG format
png("sample_distance_heatmap_star.png")
pheatmap(sampleDistMatrix)
# close PNG file after writing figure in it
dev.off() 

################################################################
#Visualize the counts with PCA ordination
################################################################

png("PCA_star.png")
plotPCA(object = vst_norm,
        intgroup = "Treatment")
dev.off()

################################################################
#Calculate Differential expression statistics between treatments
################################################################

# check results names: depends on what was modeled. Here it was the "Time"
resultsNames(dds_model)

# extract results for t25 vs t0
# contrast: the column from the metadata that is used for the grouping of the samples (Time), then the baseline (t0) and the group compared to the baseline (t25) -> results will be as "t25 vs t0"
de <- results(object = dds_model, 
              name="Treatment_mercury_vs_control")

write.table(de, "Mercury_vs_control_DE.txt", quote=F, col.names=T, row.names=T, sep="\t")

################################################################
#Calculate Differential expression statistics with shrinkage of
#LFC estimates toward 0
################################################################

# processing the same results as above but including the log2FoldChange shrinkage
# useful for visualization and gene ranking
de_shrink <- lfcShrink(dds = dds_model,
                       coef="Treatment_mercury_vs_control",
                       type="apeglm")

write.table(de_shrink, "Mercury_vs_control_DE_apeGLM.txt", quote=F, col.names=T, row.names=T, sep="\t")
