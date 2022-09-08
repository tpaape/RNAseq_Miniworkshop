# Processing Raw RNAseq Data #

## Downloading the data from JGI ##
### Identifying the necessary portals on JGI ###
- Navigate to genome.jgi.doe.gov
- Log in to your JGI account
- Press on the Genome Portal button
- 

### Downloading files to EC2 instance###
- Use the following template to login to your JGI account
- Make sure you replace username and password with your username and password:
		curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=USER_NAME' --data-urlencode 'password=USER_PASSWORD' -c cookies > /dev/null
- To download multiple files simultaneously from multiple portals using a specific file format as a pattern, use the following code as a template:
		curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/bulk/request' -b cookies --data-urlencode 'portals=portalName1,portalName2,..' --data-urlencode 'fileTypes=Sequence'
- In order to download all of the files associated with our Medicago metal RNAseq project, we will use the following script:
		curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/bulk/request' -b cookies --data-urlencode 'portals=MedtruEProfiling_2_FD,MedtruEProfiling_7_FD,MedtruEProfiling_3_FD,MedtruEProfiling_6_FD,MedtruEProfiling_4_FD,MedtruEProfiling_8_FD,MedtruEProfiling_10_FD,MedtruEProfiling_5_FD,MedtruEProfiling_9_FD' --data-urlencode 'fileTypes=Sequence'
- When you run this command, it will give you a url to check the status of your download. To check the status, use the following code as a template. Make sure you replace the example url with the url that was provided as an output
		curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/bulk/NNNNN-MM/status' -b cookies
- When the download is finished, you will get an output like this after running the code found above:
		Download request completed.
		Data URL: https://genome.jgi.doe.gov/portal/ext-api/downloads/bulk/NNNNN-MM/zip
- You can then use that link to download your files
	- Note: this will take a long time if it is more than 1-2 files. It is best to use a job handler to download these files so that the download continues if your local machine's connection gets disconnected for any reason. 
- (Option 1:) To download the files without a job handler, run the following command:
		curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/bulk/NNNNN-MM/zip' -b cookies > your_file.zip
- (Option 2:) To download the files with PBS, create a new text file by using your favorite text editor. This example will use nano
	- Create a new file called download_JGI.pbs
			nano download_JGI.pbs
	- This will open a blank text file. Copy the following into the text file:
			#Download JGI files 
			#!/bin/bash
		    	#PBS -N download_JGI
		    	#PBS -l ncpus=4
		    	#PBS -l mem=8gb
		    	#PBS -m ae
		    	#PBS -M mclear.go.olemiss.edu
		    
		    	cd $PBS_O_WORKDIR
	
			curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/bulk/NNNNN-MM/zip' -b cookies > your_file.zip

	- Note: make sure you change the following lines:
		- `#PBS -M mclear.go.olemiss.edu` so that your email address is included
		- Change the url after `wget` to your url
		- Change the filename `your_file.zip` to a name that is relevant to you but keep it as a .zip extension
	-  press cmd + x
	-  press `y`
	-  press enter
	- Run your .pbs script
			qsub download_JGI.pbs
- Once the download has finished, transfer the zip file into a new project directory
		mkdir Medicago
		mv your_file.zip Medicago/
- Unzip the files
		unzip Medicago.zip
		
## Running FastQC on the raw files ##
- Navigate into the project folder that you created
		cd Medicago
- Make a directory to deposit FastQC files in 
		mkdir FastQC_Raw_Output
- Create a new fastQC_Raw.pbs file
		nano fastQC_Raw.pbs
- This will open up an empty text file
- Insert the following code in to the fastQC_Raw.pbs text file:
		#Initial QC run     
		#!/bin/bash
    	#PBS -N QC_RAWdata
    	#PBS -l ncpus=4
    	#PBS -l mem=8gb
    	#PBS -m ae
    	#PBS -M mclear.go.olemiss.edu
    
    	cd $PBS_O_WORKDIR
    	for f in $(ls *fastq.gz)
		do
    	zcat ${f} | ~/FastQC/fastqc -t 4 -o FastQC_Raw_Output ${f}
		done 

-  press cmd + x
-  press `y`
-  press enter
-  Note: make sure you change the email portion `#PBS -M mclear.go.olemiss.edu`
- This file will loop through all of the files in the directory and put them through fastQC
- To run fastQC_Raw.pbs
		qsub fastQC_Raw.pbs
- To check the progress on the job:
		qstat -a

### Running MultiQC on the FastQC output ###
- To easily view all of the multiQC output at the same time. First navigate to the directory where you deposited all of the fastQC output files. 
		cd fastQC_Raw_Output

- Create a new multiQC.pbs file
		nano runMultiQC.pbs

- Once in the appropriate directory, make a .pbs file according to the following template:
		#Code for MultiQC.pbs
    	#!/bin/bash
		#PBS -N MultiQC
		#PBS -l ncpus=4
		#PBS -l mem=8gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		cd $PBS_O_WORKDIR

		module load py3.7

		multiqc .
- Make sure you change the email portion of the script to your email
- press cmd + x
- press `y`
- press enter
- Make sure you are in the directory where the FastQC output files are located
- To run multiQC
		qsub runMultiQC.pbs
- This will create multiple output files, but we would like to download the `multiqc_report.html` file
- Download the file to your local machine and open it in your browser to view the output

## Running Trimmomatic on the raw files ##
- Create a new file for running Trimmomatic
		nano trimmomatic.pbs
- This will open a blank text file
- Copy and paste the following text into the blank template. Make sure you change the email!:
 	   #Code for trimmomatic.pbs
    	#!/bin/bash
    	#PBS -N clean_adapters
    	#PBS -l ncpus=8
    	#PBS -l mem=16gb
    	#PBS -m ae
    	#PBS -M mclear.go.olemiss.edu
    
    	cd $PBS_O_WORKDIR
       
    	for f in $(ls *1_001.fastq.gz | sed 's/1_001.fastq.gz//' | sort -u)
    	do
    	java -jar ~/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ${f}1_001.fastq.gz ${f}2_001.fastq.gz  ${f}1_paired.fq.gz ${f}1_unpaired.fq.gz ${f}2_paired.fq.gz ${f}2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    	done

### Running FastQC and MultiQC on the trimmed and cleaned files ###
- Create a directory to deposit the new FastQC output files
- Within the `Medicago` directory:
		mkdir  FastQC_Trimmed_Output

- Make a new file for running FastQC
	nano  fastQC_trimmed.pbs

- This will open a blank text document. Copy and paste the following code into the document. Make sure you change the email address!
		#Post Trimmomatic fastQC.pbs
		#!/bin/bash
		#PBS -N trimmed_fastQC
		#PBS -l ncpus=4
		#PBS -l mem=8gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		cd $PBS_O_WORKDIR

		for f in $(ls *noadapt.fq | sed 's/noadapt.fq//')
		do
		~/FastQC/fastqc -t 4 -o FastQC_Trimmed_Output ${f}noadapt.fq
		done
- Press cmd + x
- Press `y`
- Press Enter
- To run the new FastQC script:
		qsub fastQC_trimmed.pbs
- To view the output navigate into the output directory, first copy the runMultiQC.pbs file into your new output folder
		cd FastQC_Raw_Output
		cp runMultiQC.pbs ../FastQC_Trimmed_Output

- Navigate into your new output directory
		cd ../FastQC_Trimmed_Output
- Run MultiQC
		qsub runMultiQC.pbs
- When completed, download the `multiqc_report.html` file to your local machine and view on a browser
- Make sure that sequence quality looks good and ensure that adapters have been removed

## Mapping the cleaned reads to reference genome(s) ##
### Creating a STAR indexed reference genome ###
- Create a new directory for the indexed genome
		mkdir M_trunV5

- Create a new pbs file for generating the M. truncatula indexed genome for STAR mapping
		nano medicago_STAR_index.pbs

- Copy and paste the following into the empty text file, make sure you change the email!		
		#STAR medicagoIndex.pbs
		#!/bin/bash
		#PBS -N medicago_Index
		#PBS -l ncpus=8
		#PBS -l mem=64gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		cd $PBS_O_WORKDIR

		~/STAR-2.7.10a/bin/Linux_x86_64_static/STAR --runMode genomeGenerate \
		--genomeDir ~/Medicago/M_truncV5 \
		--genomeFastaFiles ~/Genomes/MtrunA17r5.0-20161119-ANR.genome.fasta \
		--sjdbGTFfile ~/Genomes/MtrunA17r5.0-ANR-EGN-r1.8.gff3 \
		--sjdbOverhang 100 \
		--genomeSAindexNbases 13 \
		--runThreadN 8

- Close and save the file
- To generate the index file run
		qsub medicago_STAR_index.pbs

### concatenating genomes (optional) ###
### Mapping reads to indexed reference genome ###
- Create a new pbs file for mapping the Medicago reads to the indexed genome file
		nano STAR_map_medicago.pbs

- Copy and paste the following into the empty text file, make sure you change the email!
		#STAR.pbs
		#!/bin/bash
		#PBS -N STAR_medicago
		#PBS -l ncpus=8
		#PBS -l mem=32gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		cd $PBS_O_WORKDIR

		for i in $(ls *R1_sickled-noadapt.fq.gz | sed 's/R1_sickled-noadapt.fq.gz//')
		do
		~/STAR-2.7.10a/bin/Linux_x86_64_static/STAR --runMode alignReads \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--genomeDir ~/Medicago/M_truncV5 \
		--readFilesIn ${i}R1_sickled-noadapt.fq ${i}R2_sickled-noadapt.fq.gz \
		--outReadsUnmapped Fastx \
		--runThreadN 8 \
		--quantMode GeneCounts \
		--outFileNamePrefix ${i}Mtrunc \
		--alignIntronMax 3000
		done

- To run the STAR mapping script
		qsub STAR_map_medicago.pbs

## Creating count files from STAR output (*.bam) files ##
### Sorting *.bam files with samtools ###
- Open a new text file for a samsort .pbs script
		nano samSort.pbs

- Copy and paste the following into the empty text file, makes sure you change the email!
		#SAMsort.pbs
		#!/bin/bash
		#PBS -N SAM_SORT
		#PBS -l ncpus=4
		#PBS -l mem=8gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		cd $PBS_O_WORKDIR

		for i in $(ls *MtruncAligned.out.bam | sed 's/.bam//')
		do
		/usr/local/apps/sam/samtools-1.6/bin/samtools sort ${i}.bam  -n -o \
		${i}.sorted.bam
		done

- Run the samsort pbs file
		qsub samSort.pbs

### Creating count files with HTseq ###
		#HTSeq.pbs
		#!/bin/bash
		#PBS -N HTSeq
		#PBS -l ncpus=4
		#PBS -l mem=8gb
		#PBS -m ae
		#PBS -M mclear.go.olemiss.edu

		cd $PBS_O_WORKDIR

		module load py3.7

		for i in $(ls *.sorted.bam | sed 's/.sorted.bam//')
		do
		htseq-count -s no -t exon -f bam ${i}.sorted.bam ~/Genomes/ASP_CHLAM_ANNOTATION.gtf > ${i}.counts
		done
### Aggregate all of the count files with bash ###
- Create a new text file to run a bash command to aggregate the count files


### Aggregating a count summary file with bash ###
- Create a new text file to run a bash command to aggregate the STAR log files
		nano agg_STAR_log.sh

- Copy and paste the following into the empty text file
		for i in $(grep -l 'Uniquely mapped reads number' *Log.final.out); do 
			export percent=$(sed -n '/Uniquely mapped reads %/p' $i |cut -f2)
			export num=$(sed -n '/Uniquely mapped reads number/p' $i| cut -f2)
			export file=$(paste <(echo "$i") <(echo "$num") <(echo "$percent"))
			echo "$file" >> compiled_STAR_stats.txt
		done

- Save and close the file
- Make the bash script executable
		chmod u+x agg_STAR_log.sh

-While in the directory where all of the STAR output files are located run the bash script
		`./agg_STAR_log.sh`

### Aggregate the count files with Python script ###
- Open a new text file to aggregate counts
		nano aggregate_counts.py
- In the empty text document, copy and paste the following code
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

- To run the python script, make sure that you are operating in the correct conda environment
		conda activate py3.7

- Make sure that the python script is in the directory where all of the HTseq output files are located
- To move the file
		cp aggregate_counts.py /path/to/directory/htseq_counts

- Navigate to the directory where the HTseq counts are located using `cd`
- Run the python script
		python aggregate_counts.py
- This will create a file called 'Count Export.xlsx' that can be downloaded to your local machine

## Exporting the count and summary files ##
- Download the aggregated count file
		scp -i "AWS_key.pem" centos@ec2-34-227-60-128.compute-1.amazonaws.com:/home/centos/Count Export.xlsx ~/Downloads

- This will place the file into your Downloads directory
- Download the summary file
		scp -i "AWS_key.pem" centos@ec2-34-227-60-128.compute-1.amazonaws.com:/home/centos/compiled_STAR_stats.txt ~/Downloads

- This will place the file into your Downloads directory

## What to do with the remaining data ##
### Checking for contamination with phyloFlash ###