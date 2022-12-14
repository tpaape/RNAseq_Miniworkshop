# Processing Raw RNAseq Data #

## Downloading the data from JGI ##
### Identifying the necessary portals on JGI ###
- Navigate to genome.jgi.doe.gov
- Log in to your JGI account
- Press on the Genome Portal button
- Search for the projects you are interested in downloading
- Press the search button
- Check the boxes located to the left of the projects that you want to download
- Once you have selected all of the projects that you want to download, press the **Save Selected Results** button
- Select the Set name you want to save the results to or **Create New**
- Press the **My Favorites** button in the top-left corner
- Select the name of search where you saved the results and press the **Reports** button
- Press the **Download Project Overview Report** button
- This will download a .csv file of all the project information and metadata
- Open the .csv file and look at the **Portal ID** column. These are the IDs you will want to use when downloading in bulk from JGI

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
		
	-  press ctrl + x
	-  press `y`
	-  press enter
	- Run your .pbs script
	
          qsub download_JGI.pbs

- Once the download has finished, transfer the zip file into a new project directory

      mkdir Medicago
      mv your_file.zip Medicago/

- Unzip the files

      unzip Medicago.zip

## (Optional, Necessary For JGI Downloads ## 
- If files have been downloaded from JGI, they are interleaved files with both PE reads in the same file. These files need to be split before moving forward with processing. If your files are already split, they will be labed something like `sample_R1.fastq.gz` & `sample_R2.fastq.gz`. If your files are already split in to R1 and R2 or if they are SE reads, you can skip this step.
- Create a new split_Reads.pbs script

      nano split_Reads.pbs

- Copy and paste the following into the empty file:

      #Initial QC run     
      #!/bin/bash
      #PBS -N split
      #PBS -l ncpus=16
      #PBS -l mem=32gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov

      source activate py3.7
    
      cd $PBS_O_WORKDIR

      for f in $(ls *fastq.gz | sed 's/fastq.gz//')
      do
      zcat ${f}fastq.gz | reformat.sh in=${f}fastq.gz out1=${f}_R1.fastq out2=${f}_R2.fastq 
      done 

      gzip *.fastq

- Press ctrl + x
- Press `y`
- Press enter
- Navigate to the directory with all of your interleaved fastq.gz files and run the following command **(make sure to change the path to your file path)**:

      qsub PATH/to/your/split_Reads.pbs

- To check the progress on the job:

      qstat -a

## Running FastQC on the raw files ##
- Navigate into the project folder that you created where your split reads are located

      cd Medicago

- Create a new fastQC_Raw.pbs file

      nano fastQC_Raw.pbs

- This will open up an empty text file
- Insert the following code in to the fastQC_Raw.pbs text file:

      #Initial QC run     
      #!/bin/bash
      #PBS -N QC_RAW
      #PBS -l ncpus=4
      #PBS -l mem=8gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov
    
      cd $PBS_O_WORKDIR

      mkdir FastQC_Raw_Output

      for f in $(ls *.fastq.gz)
      do
      zcat ${f} | ~/FastQC/fastqc -t 4 -o FastQC_Raw_Output ${f}
      done 

-  press ctrl + x
-  press `y`
-  press enter
-  Note: make sure you change the email portion `#PBS -M mclear@bnl.gov`
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
      #PBS -M mclear@bnl.gov

      cd $PBS_O_WORKDIR

      source activate py3.7

      multiqc .

- Make sure you change the email portion of the script to your email
- press ctrl + x
- press `y`
- press enter
- Make sure you are in the directory where the FastQC output files are located
- To run multiQC

      qsub runMultiQC.pbs

- This will create multiple output files, but we would like to download the `multiqc_report.html` file
- Download the file to your local machine and open it in your browser to view the output

## Running Trimmomatic on the raw files ##
- Create a new file for running Trimmomatic

      nano run_trimm.pbs

- This will open a blank text file
- Copy and paste the following text into the blank template. Make sure you change the email!:

      #Code for trimmomatic.pbs
      #!/bin/bash
      #PBS -N clean_adapters
      #PBS -l ncpus=8
      #PBS -l mem=16gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov
    
      cd $PBS_O_WORKDIR
      
      for f in $(ls *1.fastq.gz | sed 's/1.fastq.gz//' | sort -u)
      do
      java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
      ${f}1.fastq.gz ${f}2.fastq.gz \
      ${f}1_paired.fq.gz ${f}1_unpaired.fq.gz \
      ${f}2_paired.fq.gz ${f}2_unpaired.fq.gz \
      ILLUMINACLIP:/home/centos/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
      done

### Running FastQC and MultiQC on the trimmed and cleaned files ###
- Make a new file for running FastQC

      nano  fastQC_trimmed.pbs

- This will open a blank text document. Copy and paste the following code into the document. Make sure you change the email address!

      #Initial QC run     
      #!/bin/bash
      #PBS -N QC_Trim
      #PBS -l ncpus=4
      #PBS -l mem=8gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov
    
      cd $PBS_O_WORKDIR

      mkdir FastQC_Trim_Output

      for f in $(ls *.fq.gz)
      do
      zcat ${f} | ~/FastQC/fastqc -t 4 -o FastQC_Trim_Output ${f}
      done 

- Press ctrl + x
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
### Creating a STAR indexed reference genome (*Medicago truncatula* example) ###
- Create a new pbs file for generating the M. truncatula indexed genome for STAR mapping

      nano medicago_STAR_index.pbs

- Copy and paste the following into the empty text file, make sure you change the email!		

      #STAR medicagoIndex.pbs
      #!/bin/bash
      #PBS -N medicago_Index
      #PBS -l ncpus=8
      #PBS -l mem=58gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov

      cd $PBS_O_WORKDIR

      STAR --runMode genomeGenerate \
      --genomeDir ~/Medicago/M_truncV5 \
      --genomeFastaFiles ~/Genomes/Mtrunc_v5/mtruna17r5.0-20161119-anr.genome.fasta \
      --sjdbGTFfile ~/Genomes/Mtrunc_v5/MtrunA17r5.0-ANR-EGN-r1.8.gff3 \
      --sjdbGTFtagExonParentGene ID \
      --sjdbOverhang 100 \
      --genomeSAindexNbases 13 \
      --runThreadN 8

- Close and save the file
- To generate the index file run

      qsub medicago_STAR_index.pbs

### Creating a STAR indexed reference genome (Rhizobia example) ###
- Create a new pbs file for generating the M. truncatula indexed genome for STAR mapping

      nano rhizobia_STAR_index.pbs

- Copy and paste the following into the empty text file, make sure you change the email!

      #STAR rhizobiaIndex.pbs
      #!/bin/bash
      #PBS -N SMp0_Index
      #PBS -l ncpus=8
      #PBS -l mem=58gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov

      cd $PBS_O_WORKDIR

      STAR --runMode genomeGenerate \
      --genomeDir ~/Genomes_Indexed/SMp01_Indexed \
      --genomeFastaFiles ~/Genomes/SMp01/SMp01_genome_trd.fasta \
      --sjdbGTFfile ~/Genomes/SMp01/SMp01_491863.gff \
      --sjdbOverhang 100 \
      --genomeSAindexNbases 10 \
      --sjdbGTFfeatureExon CDS \
      --sjdbGTFtagExonParentGene ID \
      --runThreadN 8

- Close and save the file
- To generate the index file run

      qsub rhizobia_STAR_index.pbs

### Concatenating genomes (optional) ###
### Mapping reads to indexed reference genome ###
- Create a new pbs file for mapping the Medicago reads to the indexed genome file

      nano STAR_map_medicago.pbs

- Copy and paste the following into the empty text file, make sure you change the email!

      #STAR.pbs
      #!/bin/bash
      #PBS -N STAR_align
      #PBS -l ncpus=8
      #PBS -l mem=32gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov

      cd $PBS_O_WORKDIR

      for i in $(ls *R1_paired.fq.gz | sed 's/R1_paired.fq.gz//')
      do
      STAR --runMode alignReads \
      --readFilesCommand zcat \
      --outSAMtype BAM Unsorted \
      --genomeDir ~/Genomes_Indexed/STf07_Indexed \
      --readFilesIn ${i}R1_paired.fq.gz  ${i}R2_paired.fq.gz \
      --outReadsUnmapped Fastx \
      --runThreadN 8 \
      --quantMode GeneCounts \
      -???sjdbGTFfile ~/Genomes/STf07/Stf07_491978.gff \
      --sjdbGTFfeatureExon CDS \
      --sjdbGTFtagExonParentGene ID \
      --outFileNamePrefix ${i}STf07 \
      --alignIntronMax 3000
      done

- Close and save the file
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
      #PBS -M mclear@bnl.gov

      cd $PBS_O_WORKDIR

      for i in $(ls *Aligned.out.bam | sed 's/.bam//')
      do
      samtools sort ${i}.bam  -n -o \
      ${i}.sorted.bam
      done

- Close and save the file
- To run the samsort pbs file:

      qsub samSort.pbs

### Creating count files with HTseq ###
- Open a new file called htseq_counts.pbs

      sudo nano htseq_counts.pbs

- Copy and paste the following into the blank file

      #HTSeq.pbs
      #!/bin/bash
      #PBS -N HTSeq
      #PBS -l ncpus=4
      #PBS -l mem=8gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov

      cd $PBS_O_WORKDIR

      source activate py3.7

      for i in $(ls *.sorted.bam | sed 's/.sorted.bam//')
      do
      htseq-count -s no -t CDS -f bam -i ID ${i}.sorted.bam /home/centos/Genomes/VTc07/Vtc07.436.gff > ${i}.counts
      done

- Close and save the file
- To run the htseq_couts file:

      qsub htseq_counts.pbs

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

      ./agg_STAR_log.sh

### Aggregate the count files with Python script ###
- Open a new text file to aggregate counts

       nano aggregate_counts.py

- In the empty text document, copy and paste the following code

      import glob
      import pandas as pd
  
      #Creates a list of filenames from files in the working directory that end in .csv
      filenames = glob.glob('*ReadsPerGene.out.tab')
  
      #creates a list of dataframes from the working directoryo
      dfs = [pd.read_csv(filename, names = ['Gene', filename], sep='\t', skiprows=4, usecols=[0,1]) for filename in filenames]
  
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
      combinedDF.to_csv('CompiledCounts.csv')
      metaDF.to_csv('MetaData.csv')

- To run the python script, make sure that you are operating in the correct conda environment

      conda activate py3.7

- Make sure that the python script is in the directory where all of the HTseq output files are located
- To move the file

      cp aggregate_counts.py /path/to/directory/htseq_counts

- Navigate to the directory where the HTseq counts are located using `cd`
- Run the python script

      python aggregate_counts.py

- This will create a files called '**CompiledCounts.csv**' & '**MetaData.csv**' that can be downloaded to your local machine

## Exporting the count and summary files ##
- Download the aggregated count file

      scp -i "AWS_key.pem" centos@ec2-34-227-60-128.compute-1.amazonaws.com:/home/centos/CompiledCounts.csv ~/Downloads

- This will place the file into your Downloads directory
- Download the summary file

      scp -i "AWS_key.pem" centos@ec2-34-227-60-128.compute-1.amazonaws.com:/home/centos/compiled_STAR_stats.txt ~/Downloads

- This will place the file into your Downloads directory

## What to do with the remaining data ##
### Checking for contamination with phyloFlash ###
- To run phyloFlash create a pbs script

      nano run_phyloFlash.pbs

- Copy and paste the following into the blank text file **(make sure to change the input file names to your files)**:

      #run_phyloFlash.pbs
      #!/bin/bash
      #PBS -N phyloFlash
      #PBS -l ncpus=8
      #PBS -l mem=32gb
      #PBS -m ae
      #PBS -M mclear@bnl.gov

      source activate py3.7

      cd $PBS_O_WORKDIR

      phyloFlash.pl -lib run01 -read1 reads_F.fq.gz \
      -read2 -reads_R.fq.gz \
      -almosteverything