# RNAseq Mini-Course #
- This repository is meant to be a resource for working with AWS EC2 instances and processing RNAseq data

## Schedule ##

### Thursday 9/7/22 ###
1.  Introduction to AWS and working with EC2 clusters
2.  Introduction to RNAseq processing (Brief Overview)
3.  Introduction to linux and bash scripting (First half of powerpoint)
4.  Extract assigned reads for processing
5.  Start processing RNAseq data (Step 1 of )
	1.  FastQC of Raw Reads
	2.  Run trimmomatic to remove adapters and trim for quality
6.  Introduction to linux and bash scripting (Second half of power point)
7.  Continue processing RNAseq data (Step 2 of )
	1.  FastQC trimmed and cleaned reads
	2.  Run STAR alignment on cleaned reads
8. Finish RNAseq processing powerpoint
9. (Time permitting) Review FastQC output:
	1. Run MultiQC on FastQC reports
	2. Look through FastQC outputs
	3. How to make STAR genome index directories
	4. Concatenating genomes for dual-transcriptomics
10. (Time permitting) Discuss PhyloFlash as a tool for troubleshooting low mapping percentage 
	1. Review example PhyloFlash output file

### Friday 9/8/22 ###
1.   Introduction to RNAseq count processing
2.   Continue processing RNAseq data with bash and R (Step 3 of )
	1.   Aggregate count summary data
	2.   Aggregate count files 
	3.   Count filtering and Normalization in R
		1.   Different methods of count normalization
	4.   Visualize sample clustering with distance matrix
	5.   Visualize sample clustering by PCA
	6.   Differential expression analysis with DESeq
	7.   GO Enrichment of Differentially Expressed Genes
3.   More advanced processing (Time permitting)
	1.   Creating heatmaps in R with ComplexHeatmap
	2.   Working with python programs within a conda environment
	3.   More advanced differential expression analysis with model building
	4.   Working with pre-exisiting R-Shiny apps
	5.   Cloning github repositories