library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(tidyr)

mat1 <- read.csv('Zscore_of_goEnrichment_Pop_trimmed.csv', row=1)

metadata <- read.csv('All_metadata.csv')

#Change "TissueXTrtXTime" string with treatment column from your metadata file
#that corresponds the names on your data matrix file
#################################################################################
metadata <- metadata %>% distinct(TissueXTrtXTime, .keep_all = TRUE) #Change this
#mat %>% select(matches())

mat1<-na.omit(mat1)
#Comment out the line below if you are working with Gene expression data
mat1 <- t(mat1) #comment out this line if working with gene expression data

mat1 <- data.matrix(mat1)

column_tree = hclust(dist(t(mat1)))
column_order <- column_tree$labels

#Change "TissueXTrtXTime" to the column that corresponds to your count labels
#######################################################################################
meta2 <- arrange(metadata, factor(metadata$TissueXTrtXTime, levels=column_order)) #change this

####################################################################################
####################################################################################
#Create vectors for heatmap annotation, 1 vector per annotation
trt <- meta2$Treatment
tiss <- meta2$Tissue
time <- meta2$Time

#Change colors and categories based on the treatment (column)
ha <- HeatmapAnnotation(
  tissue = tiss,
  treatment = trt,
  timepoint = time,
  col=list(tissue=structure(names=c("Leaf","Root"), c("green","brown")),
           treatment=structure(names=c("FeEx","Control","FeLim", "ZnEx", "ZnLim"), c("red2","royalblue","limegreen", "gold", "purple")),
           timepoint=structure(names=c("0h","1h", "2d", "4d", "7d",
                                       "14d", "21d"), brewer.pal(7,"OrRd"))),
  border=TRUE,
  show_legend=c(TRUE,TRUE,TRUE),
  show_annotation_name=FALSE,
  annotation_legend_param = list(
    tissue= list(title="Tissue Type"),
    treatment= list(title="Metal Treatment"),
    timepoint= list(title="Time Point"))
)

#Change row_km and column_km for the number of clusters you want to generate
ht_list <- Heatmap(mat1, name = "Zscore", row_km=10, column_km=7, column_gap=unit(2, "mm"),
                   bottom_annotation=ha,
                   column_names_gp=gpar(fontsize=8),
                   row_names_gp=gpar(fontsize=8))

###################################################################################
##################################################################################

#Change filename
pdf("Zscore_goEnich_full_annotated_7k_trimmed.pdf",width=12,height=12) #Change this

draw(ht_list, annotation_legend_side = "left", heatmap_legend_side = "left",
     padding = unit(c(2, 5, 2, 20), "mm"))

dev.off()

############################################################################################