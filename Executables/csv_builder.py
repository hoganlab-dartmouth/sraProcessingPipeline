# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:48:25 2020

@author: holtj

This executable script creates a csv file from the Salmon outputs.
"""

import os
import glob
import pandas as pd
import numpy as np

"""
Directories
this should be an input tbh
"""
#ref folder name, ref genomes go here
ref_folder = '/dartfs-hpc/scratch/Jake/references'
#data folder name, samples and runs go here
data = '/dartfs-hpc/scratch/Jake/data'
#csv output folder name
csv = '/dartfs-hpc/scratch/Jake/Ex'
"""
Silly function
should get rid of this
"""
#grabs stuff and returns a list of them-used throughout
def grab_ref(ref_folder1):
    files=glob.glob(ref_folder1+'/*')
    #print(files)
    return(files)

"""
Build CSV for each index
This is a copy of Alex's code, but using os and glob. 
It creates a csv file for each index with samples for rows and genes for columns. 
"""
#Grab the index names 
index_names =[] 
for y in glob.glob('INDEX_'+ref_folder+'/*'):
    index_name = y.replace('INDEX_references/','').replace('.ffn.gz','')
    print('Index: '+index_name)
    index_names.append(index_name)
    
#loop through the list of index names and build dataframe for each one
df_list = []
for i in index_names: 
    glob_list = []
    replace = 'quants/'+i+'/'
    #print(i)
    for z in grab_ref('quants/'+i):
        #print(z)
        for file in glob.glob(z+'/*.sf'):
            glob_list.append(file)
            #print(file)
    expression_df = pd.DataFrame(
    pd.read_csv(file, sep="\t", index_col=0)["TPM"].
        rename(file.replace(replace, '').replace('_quant/quant.sf',''))
    for file in glob_list)
    df_list.append(expression_df)
    #save dataframe as a csv file
    #name of csv contains reference index 
    expression_df.to_csv(csv+'/aligned_to_'+i, sep='\t')
print(df_list)

