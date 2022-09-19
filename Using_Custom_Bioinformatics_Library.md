# Using Custom Bioinformatics Library #
- The majority of my data manipulation and visualization code is written in python. 
- The easiest way to use python and manage packages is using Anaconda 
- Prior to using my bioinformatics library functions you will need to:
	1. Install Anaconda
	2. Create a conda environment
	3. Download the dependencies into your conda environment
	4. Download the BioinformaticsLibrary.py script

## Install python with Anaconda ##

### Install Anaconda Mac ###
- Navigate to [this URL](https://www.anaconda.com/downloads#macos) to download the Anaconda distribution for your Mac OS
- Double click the downloaded file and follow the prompts to install
- It is not necessary to install DataSpell after you install  Anaconda 

### Install Anaconda PC ###
- Navigate to [this URL](https://www.anaconda.com/) to download the Anaconda distribution for your OS
- Double click the downloaded file and follow the prompts to install
- It is not necessary to install DataSpell after you install Anaconda 

## Create a conda environment  and install dependencies ##
- **For Mac:** Open up a terminal window and follow the following steps where it says **"Mac/PC instructions Same"**
- **For PC:** Press the windows icon on your computer
	- Type in anaconda
	- The Anaconda Prompt should pop up
	- Select the Anaconda Prompt program
	- Follow the following steps
- **Mac/PC instructions Same:**
- Create a conda environment

      conda create -n BioLib python=3.9

- Activate your conda environment

      conda activate BioLib

- Install bioconda dependencies

      conda install -c bioconda gprofiler-official
      
- Install anaconda dependencies

      conda install -c anaconda pandas numpy scipy spyder seaborn requests

- Install conda-forge dependencies

      conda install -c conda-forge matplotlib

## Download BioinformaticsLibrary ##
- Download the [BioinformaticsLibrary.py](https://github.com/mclear73/RNAseq_Miniworkshop/blob/main/Python%20Scripts/BioinformaticsLibrary.py) script from github
- Make sure that this python is located in the working directory where plan on running your analyses

## Typical python use case ##
- You will always activate python from either your **terminal (Mac)** or **Anaconda prompt (PC)**
- Navigate to the appropriate prompt based on your OS (terminal or anaconda prompt)
- Activate your conda environment

      conda activate BioLib

- Run the following command:

      spyder

- This will open Spyder within your conda environment
- Spyder is setup nearly identically to RStudio and should feel very similar
- Change your working directory to where you downloaded the BioinformaticsLibrary.py script
- Every python script you write should have the following line at the top:

      import BioinformaticsLibrary as BL

- This will import all of the functions I have written. Here is an example for how you would call one of the functions:

      BL.compile_Counts(directory="<insert/directory/of/files/here>")