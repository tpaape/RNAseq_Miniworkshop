library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(tidyr)

mat1 <- read.csv('Zscore_of_goEnrichment_Pop_trimmed.csv', row=1)
mat2 <- read.csv('Zscore_of_goEnrichment_Pop_FeOnly.csv', row=1)
mat3 <- read.csv('Zscore_of_goEnrichment_Pop_ZnOnly.csv', row=1)
mat4 <- read.csv('Zscore_of_goEnrichment_Pop_FeOnly_leaf.csv', row=1)
mat5 <- read.csv('Zscore_of_goEnrichment_Pop_FeOnly_root.csv', row=1)
mat6 <- read.csv('Zscore_of_goEnrichment_Pop_ZnOnly_leaf.csv', row=1)
mat7 <- read.csv('Zscore_of_goEnrichment_Pop_ZnOnly_root.csv', row=1)

metadata <- read.csv('All_metadata.csv')

metadata <- metadata %>% distinct(TissueXTrtXTime, .keep_all = TRUE)
#mat %>% select(matches())



###############################################################################################
mat1<-na.omit(mat1)
mat1 <- t(mat1)

mat1 <- data.matrix(mat1)

column_tree = hclust(dist(t(mat1)))
column_order <- column_tree$labels

meta2 <- arrange(metadata, factor(metadata$TissueXTrtXTime, levels=column_order))

trt <- meta2$Treatment
tiss <- meta2$Tissue
time <- meta2$Time

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


ht_list <- Heatmap(mat1, name = "Zscore", row_km=10, column_km=7, column_gap=unit(2, "mm"),
                   bottom_annotation=ha,
                   column_names_gp=gpar(fontsize=8),
                   row_names_gp=gpar(fontsize=8))

pdf("Zscore_goEnich_full_annotated_7k_trimmed.pdf",width=12,height=12)

draw(ht_list, annotation_legend_side = "left", heatmap_legend_side = "left",
     padding = unit(c(2, 5, 2, 20), "mm"))

dev.off()

############################################################################################