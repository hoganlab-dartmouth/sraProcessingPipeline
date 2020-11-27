# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 08:34:36 2020

@author: holtj

The directories are probbaly a little messed up for running on Discovery. 

This single_quantifier script handles the single SRR number per SRX/BioSample# case. 
Use the other one if you want to combine samples. 

This script downloads data, feeds it into Salmon using the correct 
parameters and deletes downloaded data after Salmon is finished.
To do this, it calls SRAToolkit fastq-dump and Salmon quant. os.system() is used
to call this and other command line functions. 

Command line input when calling this script should be: 
python quantifier.py -l "SRX,SRR,SRR,SRR"

I would like to freeze this script as an executable to get around having to add 
the libraries to Discovery. 

The paths on discovery might mess this up. On Discovery we want to use 
/global/scratch since we are cleaning things up after. 
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
input_lis = [int(item) for item in args.list.split(',')]

"""
Directories

Tells the script which folder your references, data, and outputs are stored in.

This should probably be an input? It's used in every script. 
"""
#ref folder name, ref genomes go here
ref_folder = '/dartfs-hpc/scratch/Jake/references'
#data folder name, samples and runs go here
data = '/dartfs-hpc/scratch/Jake/data'
#csv output folder name
csv = '/dartfs-hpc/scratch/Jake/Ex'

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
for y in glob.glob('INDEX_'+ref_folder+'/*'):
    index_name = y.replace('INDEX_references/','').replace('.ffn.gz','')
    print('Index: '+index_name)
    dump_directory = 'data/'+srx+'/'+srr
    os.system('mkdir '+dump_directory)
    print('Dump Directory: '+dump_directory)
    os.system('fastq-dump --outdir '+dump_directory+' -v --gzip --skip-technical --readids --split-3 --clip '+srr)
    run_folders = grab_ref('data/'+srx)
    output_name='quants/'+index_name+'/'+srx+'_quant'
    print('Quant output:'+output_name)
    if len(run_folders)<1:
        print('Oops, no runs for '+srx)
    if len(run_folders)==1:
        if len(grab_ref(run_folders[0]))<1:
            print('Oops, no runs for '+srx)
        elif len(grab_ref(run_folders[0]))==1:
            print('One run for '+srx)
            input_name='-r'
            input_name +=' '+grab_ref(run_folders[0])[0]
            print('Salmon input: '+input_name)
            os.system('salmon quant -i '+y+' -l A '+input_name+' -p 8 --validateMappings -o '+output_name)
        elif len(grab_ref(run_folders[0]))>1:
            runs = grab_ref(run_folders[0])
            print('Multiple runs for ')
            input_name='-r'
            for z in np.arange(0, len(runs)):
                input_name+= ' '+runs[z]
            print('Salmon Input: '+input_name)
        os.system('salmon quant -i '+y+' -l A '+input_name+' -p 8 --validateMappings -o '+output_name)
        os.system('rm -r data/'+srx)