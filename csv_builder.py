# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:48:25 2020

@author: holtj

This executable script creates a csv file from the Salmon outputs.
It relies on the reference folder and quant folder being accessible and having 
the right internal folder structures. See write up and graph. 
"""

#We're using glob to grab filenames. I know this isn't ideal. 
import glob
#Pandas is used to write to csv. It was easy this way. 
import pandas as pd

"""
Directories used by the program. 
This should be an input tbh.
Also should probably not do scratch. 
"""
#ref folder name, ref genomes go here
ref_folder = 'References'
#csv output folder name
csv = ''

"""
Silly function. Uses glob to rgab file names and then edits the string a lil. 
I should rewrite so that this is done inline. 
"""
#grabs stuff and returns a list of them-used throughout
def grab_ref(ref_folder1):
    files=glob.glob(ref_folder1+'/*')
    return(files)

"""
Build CSV for each index
This is a copy of Alex's code, but using os and glob. 
It creates a csv file for each index with samples for rows and genes for columns. 
"""
#Grab the index names 
index_names =[] 
for file_name in glob.glob('INDEX_'+ref_folder+'/*'):
    index_name = file_name.replace('INDEX_references/','').replace('.ffn.gz','')
    index_names.append(index_name)
    
#loop through the list of index names and build dataframe for each one
#directory structure is really important. See powerpoints/journal. 
df_list = []
for index in index_names: 
    glob_list = []
    replace = 'quants/'+index+'/'
    for z in grab_ref('quants/'+index):
        for file in glob.glob(z+'/*.sf'):
            glob_list.append(file)
    expression_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
        rename(file.replace(replace, '').replace('_quant/quant.sf',''))
    for file in glob_list)
    print(len(expression_df))
    df_list.append(expression_df)
    #save dataframe as a csv file
    #name of csv contains reference index 
    expression_df.to_csv('aligned_to_'+index.replace("/","_"), sep='\t')

