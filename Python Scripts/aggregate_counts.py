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