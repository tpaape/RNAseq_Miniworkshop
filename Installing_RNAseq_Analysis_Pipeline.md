# Installing RNAseq Analysis Pipeline #

- This has been written assuming CentOS7 is the operating system and CrowdStrike has already been installed
- If you need to install CrowdStrike, refer to the end of the [Creating_EC2_instance](https://github.com/mclear73/RNAseq_mini_workshop/blob/main/Creating_EC2_instance.md) instructions

## Install dependencies ##
- Install wget, nano, git, and java
		sudo yum install -y wget nano git  java-1.8.0-openjdk

- Install openPBS dependencies (part 1)
	  sudo yum install -y gcc make rpm-build libtool hwloc-devel \
  	libX11-devel libXt-devel libedit-devel libical-devel \
  	ncurses-devel perl postgresql-devel postgresql-contrib python3-devel tcl-devel \
  	tk-devel swig expat-devel openssl-devel libXext libXft \
  	autoconf automake gcc-c++

- Install openPBS dependencies (part 2)
	  sudo yum install -y expat libedit postgresql-server postgresql-contrib python3 \
      sendmail sudo tcl tk libical

- Install openPBS dependencies (part 3)
		sudo yum install -y postgresql-server postgresql-contrib

## Install openPBS ##
- Clone openPBS repository from github
		git clone https://github.com/openpbs/openpbs.git

- Move into the openpbs directory
		cd openpbs

- Run the following commands:
		#First run this
		./autogen.sh

		#Then this
		./configure --prefix=/opt/pbs

		#Then this
		make

		#Finally, this
		sudo make install

- Now run this post-install command:
		sudo /opt/pbs/libexec/pbs_postinstall

- Change config file
		sudo nano /etc/pbs.conf

- Edit the config file line: `PBS_START_MOM=0` to `PBS_START_MOM=1`
- - Save and close the pbs.conf file
- Change file permissions:

		sudo chmod 4755 /opt/pbs/sbin/pbs_iff /opt/pbs/sbin/pbs_rcp

- Start PBS services:

		sudo /etc/init.d/pbs start

## Install Anaconda (python) ##
- Install anaconda
		curl -O https://repo.anaconda.com/archive/Anaconda3-5.3.1-Linux-x86_64.sh

		bash Anaconda3-5.3.1-Linux-x86_64.sh

- Follow prompts for installation

## Install FastQC ##
- Download the install file
		wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
- Unzip the file
		unzip fastqc_v0.11.9.zip
- Make fastqc executable
		cd FastQC

		chmod 755 fastqc

- Permanently add fastqc to your path by open the following file:

		sudo nano ~/.bashrc
- At the bottom of the file add the following line:
		export PATH="/home/centos/FastQC:$PATH"

-If you close your current terminal and reconnect, you should be able to run FastQC by running the command `fastqc`


## Install MultiQC ##
- Create a conda environment for multiQC
		conda create --name py3.7 python=3.7
		conda activate py3.7
- Install MultiQC
		conda install -c bioconda/label/cf201901 multiqc bbmap

## Install Trimmomatic ##
- Download binary
		wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
- Unzip the binary
		unzip Trimmomatic-0.39.zip
- You can now run Trimmomatic by running `java -jar Trimmomatic-0.39/trimmomatic-0.39.jar`

## Install STAR ##
- Download the most recent version
		wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz

- Unpack the file
		tar -xzf 2.7.10a.tar.gz

- Navigate to source directory
		cd STAR-2.7.10a/source

- Install
		make STAR

- Permanently add fastqc to your path by open the following file:
		sudo nano ~/.bashrc

- At the bottom of the file add the following line:
		export PATH="/home/centos/STAR-2.7.10a/bin/Linux_x86_64_static:$PATH"

-If you close your current terminal and reconnect, you should be able to run STAR by running the command `STAR`

## Install samtools ##
- Install samtools dependencies
		sudo yum group install -y "Development Tools"
		sudo yum install -y ncurses-devel bzip2-devel xz-devel

- Move to your home directory
		cd
- Download samtools
		wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2

- Unpack the file
		tar -xvjf samtools-1.9.tar.bz2

- Navigate into the samtools directory
		cd samtools-1.9

- Configure for install
		./configure --prefix=/usr/local

- Compile and install
		make
		sudo make install

- Navigate to your home directory
		cd

## Install HTseq and other count processing tools ##
- If not already activated, activate the py3.7 environment
		conda activate py3.7
- Install HTseq
		conda install -c bioconda htseq
- Install pandas
		conda install -c anaconda pandas

## Install PhyloFlash (optional) ##
- Note: this install can be time consuming and difficult. Only attempt if you NEED to use phyloFlash

- Activate your conda environment
		conda activate py3.7

- Install dependencies from bioconda:		
		conda install -c bioconda bowtie vsearch spades mafft bedtools

- Install dependencies from Anaconda:
		conda install -c anaconda scipy cython

- Install dependencies from conda-forge:
		conda install -c conda-forge biopython

- Download phyloFlash
		wget https://github.com/HRGV/phyloFlash/archive/pf3.4.tar.gz
		tar -xzf pf3.4.tar.gz

- Permanently add fastqc to your path by open the following file:
		sudo nano ~/.bashrc

- At the bottom of the file add the following line:
		export PATH="/home/centos/phyloFlash-pf3.4:$PATH"

-If you close your current terminal and reconnect, you should be able to run phyloFlash by activating your conda environment `conda activate py3.7` and  running the command `phyloFlash.pl -check_env`

- Create the database locally
		phyloFlash_makedb.pl --remote

- This step will take a while. Make sure you leave time for this installation to take place. You will occasionally be prompted to accept downloads etc.

- Confirm that your phyloFlash environment is operational by running the command:
		phyloFlash.pl -check_env

- We did not install EMIRGE as part of this process. EMIRGE has proven difficult to install successfully so make sure you do not attempt to use EMIRGE during your phyloFlash analyses

		