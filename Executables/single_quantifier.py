# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 08:34:36 2020

@author: holtj

This single_quantifier script handles the single SRR number per SRX/BioSample# case. 
Use the other one if you want to combine samples. 

This script downloads data, feeds it into Salmon using the correct 
parameters and deletes downloaded data after Salmon is finished.
To do this, it calls SRAToolkit fastq-dump and Salmon quant. os.system() is used
to call this and other command line functions. 

Command line input when calling this script should be: 
python quantifier.py -l "SRR,SRR" -> this was a rushed decision. Originally I was doing SRX, SRR
but because of time and some issues that popped up with that, I just switched to doing SRR,
Please note the variable names are now messed up. SRX is actually an SRR and SRR is an SRR number. 

"""

#os is used for command line interface
import os
#glob is used to find files by name
import glob
#numpy is used to iterate through lists. Could probably rewrite to not use numpy
import numpy as np
#argparse is used for comman line input 
import argparse

"""
Using arparse to take in command line arguments

Command line input when calling this script should be: 
python quantifier.py -l "SRX,SRR,SRR,SRR"
"""
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--list', help='delimited list input', type=str)
args = parser.parse_args()
input_lis = [str(item) for item in args.list.split(',')]

"""
Directories

Tells the script which folder your references, data, and outputs are stored in.

This should probably be an input? It's used in every script. 
"""
#ref folder name, ref genomes go here
ref_folder = 'References'
#data folder name, samples and runs go here
data = '/dartfs-hpc/scratch/f004ky5/data'
if os.path.isdir(data)==False:
  os.system('mkdir '+data)
#csv output folder name
csv = 'Ex'

"""
Silly function. Should rewrite to not use this. 
"""
#grabs stuff and returns a list of them. This is used throughout scripts.
def grab_ref(ref_folder1):
    files=glob.glob(ref_folder1+'/*')
    return(files)

"""
Downloads data using fastq-dump, quantifies data using Salmon, and deletes 
downloaded data. 
"""
srx = input_lis[0]
srr = input_lis[1]
for index in glob.glob('INDEX_'+ref_folder+'/*'):
    index_name = index.replace('INDEX_references/','').replace('.ffn.gz','')
    print('Index: '+index_name)
    dump_directory = data+'/'+srr
    os.system('mkdir '+dump_directory)
    print('Dump Directory: '+dump_directory)
    os.system('fastq-dump --outdir '+dump_directory+' --gzip --skip-technical --readids --split-e --clip '+srr)
    run_files = grab_ref(data+'/'+srr)
    print('SRR Files: '+run_files[0])
    output_name='quants/'+index_name+'/'+srr+'_quant'
    print('Quant output:'+output_name)
    if len(run_files)<1:
        print('Oops, no runs for '+srr)
    elif len(run_files)==1:
      print('One run for '+srr)
      input_name='-r'
      input_name +=' '+run_files[0]
      print('Salmon input: '+input_name)
      os.system('salmon quant -i '+index+' -l A '+input_name+' -p 1 --validateMappings -o '+output_name)
    elif len(run_files)>1:
      print('Multiple runs for '+srr)
      input_name='-r'
      for z in np.arange(0, len(run_files)):
        input_name+= ' '+run_files[z]
      print('Salmon Input: '+input_name)
      os.system('salmon quant -i '+index+' -l A '+input_name+' -p 1 --validateMappings -o '+output_name)
os.system('rm -r '+data+'/'+srr)