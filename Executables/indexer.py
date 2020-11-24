# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 07:37:53 2020

@author: holtj

This executable script builds the indexes for Salmon. 
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
ref_folder = '/dartfs-hpc/scratf004ky5ch/Jake/references'
#data folder name, samples and runs go here
data = '/dartfs-hpc/scratch/Jake/data'
#csv output folder name
csv = '/dartfs-hpc/scratch/Jake/Ex'
"""
Silly function
should get rid of this
"""
def grab_ref(ref_folder1):
    files=glob.glob(ref_folder1+'/*')
    #print(files)
    return(files)
"""
Build the indexes
"""
#builds Salmon index from given file names, puts indexes into folder, INDEX_references 
def salmon_index(genome):
    os.system('salmon index -t '+str(genome)+' -i INDEX_'+str(genome))
    #print('salmon '+ genome)
    return 

#call grab_ref() and salmon_index() to build index for each reference genome
for i in grab_ref(ref_folder):
    salmon_index(i)
    
