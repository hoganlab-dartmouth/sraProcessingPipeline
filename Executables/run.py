# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:44:29 2020

@author: holtj

This executable is the pipeline management. It calls the others and submits jobs. 
"""

import os
import glob
import pandas as pd
import numpy as np
"""
Make sure Salmon and SRAToolkit are installed and configured properly
"""


"""
Setup directories
"""
#ref folder name, ref genomes go here
ref_folder = 'references'
os.system('mkdir '+ref_folder)
#data folder name, samples and runs go here
data = 'data'
os.system('mkdir '+data)
#csv output folder name
csv = 'Ex'
os.system('mkdir '+csv)

"""
Create dictionary
"""
#read in SRR#s from metadata table provided by Georgia
sra_run_table = pd.read_csv(r'SraRunTable.csv', sep=",", header=0, index_col=0)
#print(sra_run_table['Run'][3:7])
#create a dict key of biosamples : srr numbers(could do SRX:SRR instead)
bio_sample_dic = sra_run_table.groupby('BioSample').agg({'Run':lambda x:x.tolist()})['Run'].to_dict()

"""
Call the indexer executable
"""

"""
Start submitting jobs using exec #2
"""
for i in bio_sample_dic.keys():
    #build string input for quantifier.py
    #write pbs file
    #qsub pbs ???
    print(i)
    
"""
Call the csv_builder executable
"""
