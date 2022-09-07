**File Download via ftp**



1. Log onto maple on the cluster and use the wget function to collect the data with the following command:

    	wget -r -np -R “index.html” --user <user> --password <pass> <http://example.com/>



	- This script is meant to grab all of the files from a particular directory if you want to download individual files then you need to delete the "-r -np -R "index.html" portion of the code



1.  Once all of the file have been downloaded use fastQC to **QC all of the data**
1.  Make a directory called FastQC_Output with the following command:

    	mkdir FastQC_Output

1.  In maple write the following as a .pbs script:
	#Initial QC run     
		#!/bin/bash

    	#PBS -N QC_RAWdata

    	#PBS -l ncpus=4

    	#PBS -l mem=8gb

    	#PBS -m ae

    	#PBS -M mclear.go.olemiss.edu
    
    	umask 007
    	cd $PBS_O_WORKDIR
    	for f in $(ls *fastq.gz)
		do
    	zcat ${f} | ~/ptmp/FastQC/fastqc -t 4 -o FastQC_Output ${f}
		done 
    	
	1. This will loop through all of the files in the directory and put them through fastQC
2. Take a look at the fastQC files using multiQC. It is beneficial to get an initial look at your samples because you can check depth, quality, etc. But most importantly, it can help give you an idea of the appropriate adapter file that you need to include in the adapter trimming step. 
	1. To easily view all of the multiQC output at the same time. First navigate to the directory where you deposited all of the fastQC output files. 
	2. Once in the appropriate directory, make a .pbs file according to the following code:
	#Code for MultiQC.pbs
    	#!/bin/bash
		#PBS -N MultiQC
		#PBS -l ncpus=4
		#PBS -l mem=8gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		umask 007
		cd $PBS_O_WORKDIR

		module load python

		multiqc .
	1. This will compile all of the fastQC files sitting in the current directory
	2. view the files and check the adapter sequences and compare that to the sequences sitting in the trimmomatic folder to select the correct adapter file
3. **Trim the adapters** using Trimmomatic and the following code. **Ensure that you change the adapter file depending on the correct adapters you identified using fastQC/multiQC and/or google**
 	#Code for trimmomatic.pbs
    	#!/bin/bash
    	#PBS -N clean_adapters
    	#PBS -l ncpus=8
    	#PBS -l mem=16gb
    	#PBS -m ae
    	#PBS -M mclear.go.olemiss.edu
    
    	umask 007
    	cd $PBS_O_WORKDIR
    
    	module load java
    
    	for f in $(ls *1_001.fastq.gz | sed 's/1_001.fastq.gz//' | sort -u)
    	do
    	java -jar ~/ptmp/StickingTransJune19/AspergillusLaneA/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ${f}1_001.fastq.gz ${f}2_001.fastq.gz  ${f}1_paired.fq.gz ${f}1_unpaired.fq.gz ${f}2_paired.fq.gz ${f}2_unpaired.fq.gz ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10
    	done
	1. Note that this code is pointing to a directory labeled "adapters" you need to go into the Trimmomatic source code and copy the adpaters directory into the directory where your .pbs file and reads are found.
		1. Also ensure that you are using the correct adapters from the adapter folder
		2. The trailing three numbers are for the following:
			1. 2 -> seed mismatches: maximum mismatch count whichwill still allow full match to be performed
			2. 30 -> palindrome clip threshhold: how accurate the match between the two 'adapterligated' reads must be for PE palindrome read alignment
			3. 10 -> simple clip threshhold: how accurate the match between ay adapter sequecne must be against a read
			4. It appears that trimmomatic is commonly used for adapter trimming and QC of raw reads for C. reinharditt see Zones et al., 2015. Default settings are also commonly used for N. crassa and A. nidulans
5. Run fastQC and multiQC on these files again to make sure that the adapters have been trimmed. If not change the adpater files in the trimmomatic .pbs script
	6. **Trim reads for quality using Sickle**
		1. While still in the same directory with all fo the trimmed reads run the following .pbs scrip:

		#Code for sickle.pbs
			#!/bin/bash
			#PBS -N sickle
			#PBS -l ncpus=4
			#PBS -l mem=8gb
			#PBS -m ae
			#PBS -M mclear.go.olemiss.edu

			umask 007
			cd $PBS_O_WORKDIR

			for f in $(ls *R1_paired.fq.gz | sed 's/R1_paired.fq.gz//' | sort -u)
			do
			sickle pe -f ${f}R1_paired.fq.gz -r ${f}R2_paired.fq.gz \
			-t sanger -o ${f}R1_sickled-noadapt.fq -p ${f}R2_sickled-noadapt.fq \
			-s ${f}sickle-noadapt_singletons.fq
			done


		1. This will loop through all of the trimmomatic output files and clean the reads based on the default sickle settings. Note that sickle doesn't allow manipulation of its trim settings. This technique was also employed in Ban et al., 2019 with C. reinhardtti. 
	2. **Final Read QC** with FastQC. 
		1. Use the following code to ensure that the reads are now an acceptible quality for mapping:
	#Final fastQC.pbs
		#!/bin/bash
		#PBS -N final_fastQC
		#PBS -l ncpus=4
		#PBS -l mem=8gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		umask 007
		cd $PBS_O_WORKDIR

		for f in $(ls *noadapt.fq | sed 's/noadapt.fq//')
		do
		~/ptmp/FastQC/fastqc -t 4 -o FastQC_Output ${f}noadapt.fq
		done
1. **MultiQC** is a helpful tool for visualizing fastQC data to run MultiQC, perform the following:
	1. Navigate to your FastQC_Output directory by:

			cd FastQC_Output

	1. Run the following .pbs script:
	#MultiQC.pbs
		#!/bin/bash
		#PBS -N MultiQC
		#PBS -l ncpus=4
		#PBS -l mem=8gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		umask 007
		cd $PBS_O_WORKDIR

		module load python

		multiqc .
	1. This will output a .html that can be viewed in your browser

2. **Creating Indexed Genomes for STAR**
 

	1. Prior to mapping the trimmed reads to the fungal/algal genomes, we need to create indexed genomes for each pairing. For this experiment I used the following annotations/assemblies from JGI:
		
		1. C. reinhardtii: v5.6 -> Creinhardtii_281_v5.0.fa **&** Creinhardtii_281_v5.6.gene_exons.gff3
		2. A. nidulans: Aspnid1 -> Aspnid1_GeneCatalog_CDS_20110130.fasta **&** Aspnid1.AspGD_genes.gff3
		3. N. crassa: (Neurospora crassa OR74A v2.0) Neucr2 -> Neucr2_GeneCatalog_CDS_20130412.fasta **&** Neucr2.BroadModels.gff3
 
	1. First, we need to concatenate our two assembly files and two annotation files for the corresponding fungus with algae. To do this, we first need to remove the head from the fungus annotations
	#Remove head from fungus files
		#To remove the top two lines from a gff3 file, use the code below. 
		#to remove a different number of header lines, change the number to 1+ the number of lines you want to remove		
		tail -n +3 file.gff3 > file_NOHEAD.gff3
	
	1. Now, we can concatenate our files:
	#Concatenate assemblies and annotations
		cat file1.gff3 file2_NOHEAD.gff3 >> combined.gff3
		cat file1.fasta file2.fasta >> combined.fasta
	1. Now, run the following .pbs scripts to generate combined indexed genomes as well as individual indexed genomes for each organism
	#Combined Indexed genome
		mkdir Asp_COMB_Genome
		
	#STAR Index.pbs
		#!/bin/bash
		#PBS -N STAR2.7_COMB_Index
		#PBS -l ncpus=4
		#PBS -l mem=8gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		umask 007
		cd $PBS_O_WORKDIR

		~/STAR-2.7.1a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
		--genomeDir ~/ptmp/StickingTransJune19/Assemblies_Annotations/Aspergillus/Asp_COMB_Genome \
		--genomeFastaFiles ~/ptmp/StickingTransJune19/Assemblies_Annotations/Aspergillus/ASP_CHLAM_ASSEMBLY.fa \
		--sjdbGTFfile ~/ptmp/StickingTransJune19/Assemblies_Annotations/Aspergillus/ASP_CHLAM_ANNOTATION.gff3 \
		--sjdbOverhang 100 \
		--genomeSAindexNbases 13 \
		--runThreadN 8
	1. Modify the code above to generate indexes for N. crassa, A. nidulans, & C. reinhardtii individually as well as the concatenated N. crassa & C. reinhardtii index
	2. PLEASE change the genomeSAindexNbases setting as you work with different genomes. "13" is used for the combined genomes. the equation for determining this number is log2(genomeLength)/2-1
		1. N. crassa -> 12, genome size: 40Mb
		2. A. nidulans -> 12, genome size: 29.8Mb
		3. C. reinhardtii -> 13, genome sizeL: 110 Mb
1. Now we are ready to **Map sequences to Genomes** with STAR. To do so use the following .pbs script:
	#STAR.pbs
		#!/bin/bash
		#PBS -N STAR2.7_COMB
		#PBS -l ncpus=8
		#PBS -l mem=32gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		umask 007
		cd $PBS_O_WORKDIR

		for i in $(ls *R1_sickled-noadapt.fq | sed 's/R1_sickled-noadapt.fq//')
		do
		~/STAR-2.7.1a/bin/Linux_x86_64/STAR --runMode alignReads \
		--outSAMtype BAM Unsorted \
		--genomeDir ~/ptmp/StickingTransJune19/Assemblies_Annotations/Neurospora/Asp_COMB_Genome \
		--readFilesIn ${i}R1_sickled-noadapt.fq ${i}R2_sickled-noadapt.fq \
		--runThreadN 8 \
		--outFileNamePrefix ${i}COMB
		--alignIntronMin 20
		--alignIntronMax 3000
		done
	1. Notice the location of the STAR file. This will need to change based on your system
	2. the last two intron related settings were taken from Zones et al. 2015. These settings don't noticeably change the mapping results compared to default
	3. This mapping process needs to be repeated for the different genome indexes (i.e. fungus only & chlamy only)
3. Now we first need to sort our output files from STAR so that we can input them into HTSeq
	#SAMsort.pbs
		#!/bin/bash
		#PBS -N SAM_SORT
		#PBS -l ncpus=4
		#PBS -l mem=4gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		umask 007
		cd $PBS_O_WORKDIR

		for i in $(ls *COMBAligned.out.bam | sed 's/.bam//')
		do
		/usr/local/apps/sam/samtools-1.6/bin/samtools sort ${i}.bam  -n -o \
		${i}.sorted.bam
		done


1. Finally, we can get gene counts using **HTSeq**
	#HTSeq.pbs
		#!/bin/bash
		#PBS -N HTSeq
		#PBS -l ncpus=4
		#PBS -l mem=4gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		umask 007
		cd $PBS_O_WORKDIR

		module load python/2.7
		source activate htseq

		for i in $(ls *.sorted.bam | sed 's/.sorted.bam//')
		do
		htseq-count -s no -t exon -f bam ${i}.sorted.bam ~/ptmp/StickingTransJune19/Assemblies_Annotations/Aspergillus/ASP_CHLAM_ANNOTATION.gtf > ${i}.counts
		done
	1. Note that we are using the combined *.gtf* file for this step this is because the id labels do not agree in the .gff3 files so we had to also generate concatenated .gtf files 
	2. -s no -> indicates that this data is **not** from a strand specific library
	2. -t exon -> indicates that the key word describing exons is exon
	3. -f bam -> indicates that we are inputing .bam files
	4. This will export .count files for every sample. These files can be opened and modified in excel
1. Now, we want to run the following python file that will compile all of our count data automatically as well as generate an initial metadata file
	#Step1_Count compilation.py
		import glob
		import pandas as pd

		#Creates a list of filenames from files in the working directory that end in .csv
		filenames = glob.glob('*.counts')

		#creates a list of dataframes from the working directoryo
		dfs = [pd.read_csv(filename, names = ['Gene', filename], sep='\t') for filename in filenames]

		#Combines all dataframes into one
		combinedDF = pd.concat(dfs, axis=1)
		#Drops duplicate gene columns that are remnants from the concat function
		#This works by transposing the dataframe, dropping duplicates, then re-transposing
		combinedDF = combinedDF.T.drop_duplicates().T
		#Sets the index to the gene column
		combinedDF = combinedDF.set_index('Gene')

		#Starts a dataFrame with the filenames as the index
		metaDF = pd.DataFrame(filenames).rename(columns={0:'Sample'})
		metaDF = metaDF.set_index('Sample')

		#Saves the dataframes in a single excel file in the current directory
		writer = pd.ExcelWriter('Count Export.xlsx')
		combinedDF.to_excel(writer,'Counts')
		metaDF.to_excel(writer,'MetaData')
		writer.save()
	1. Notice that this will compile all of the .counts files in a single directory. For this reason, it is important to have a new directory for every experiment. 
	2. It is expected that the user will populate the metadata sheet with the appropriate information manually in excel. To do so, use the following template.


