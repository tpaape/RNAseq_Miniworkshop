# RNAseq Mini-Course #
- This repository is meant to be a resource for working with AWS EC2 instances and processing RNAseq data
- How to clone this repository (and other private repositories you have access to):

      git clone https://username@github.com/mclear73/RNAseq_Miniworkshop.git

- Make sure to change "username" to your username

## Moving data in and out of AWS S3 buckets ##
- It is expensive to stop instances if the only purpose of maintaining the instance is to store data. The most efficient AWS option for storing data is using an S3 bucket. 
- I have created a bucket called **qpsi-team-bucket** that can be used for temporarily storing data
- To transfer **from EC2 to S3** use the following template while on your ec2 instance:

      aws s3 cp testFile.txt s3://qpsi-team-bucket

- To transfer **from S3 to EC2** use the following template while on your ec2 instance:

      aws s3 cp s3://qpsi-team-bucket/nph18078-sup-0002-tabless1-s11.xlsx ~/Temp/

- This will create a new folder, "Temp" and download the nph...xlsx file there. Alternatively, you can specify an existing directory and download the file into there
- **To transfer entire directories,** there are 2 methods:
- Sync:

      aws s3 sync native_Directory s3://qpsi-team-bucket

- Recursive:

      aws s3 cp my-local-folder s3://qpsi-team-bucket/ --recursive

- To **transfer files based on wild cards:**

      aws s3 cp s3://qpsi-team-bucket/ . --recursive --exclude "*" --include "*.fastq.gz"

- This will copy all files ending in .fastq.gz from the qpsi-team-bucket to your current directory on your cluster

      aws s3 cp /Temp/sub_temp/ s3://qpsi-team-bucket/ --recursive --exclude "*" --include "*fq"

- This will copy all files from the sub_temp directory on your EC2 instance that end in "fq" to the qpsi-team-bucket bucket

## Schedule ##

### Thursday 9/7/22 ###
1.  Introduction to AWS and working with EC2 clusters
2.  Introduction to RNAseq processing (Brief Overview)
3.  Introduction to linux and bash scripting (First half of powerpoint)
4.  Extract assigned reads for processing
5.  Start processing RNAseq data (Step 1 of  3)
	1.  FastQC of Raw Reads
	2.  Run trimmomatic to remove adapters and trim for quality
	3.  Create STAR genome indices
6.  Introduction to linux and bash scripting (Second half of power point)
7.  Continue processing RNAseq data (Step 2 of 3)
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
2.   Continue processing RNAseq data with bash and R (Step 3 of 3)
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