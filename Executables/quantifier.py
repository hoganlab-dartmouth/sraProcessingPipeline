# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 09:12:51 2020

@author: holtj

This script downloads data, feeds it into Salmon using the correct 
parameters and deletes downloaded data after Salmon is finished.
To do this, it calls SRAToolkit fastq-dump and Salmon quant. os.system() is used
to call this and other command line functions. 

Command line input when calling this script should be: 
python quantifier.py -l "SRX,SRR,SRR,SRR"

I would like to freeze this script as an executable to get around having to add 
the libraries to Discovery. 

The paths on discovery might mess this up. On Discovery we want to use 
/global/scratch.

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
ref_folder = 'references'
#data folder name, samples and runs go here
data = 'data'
#csv output folder name
csv = 'Ex'

"""
Silly function. Should rewrite to not have this. 
"""
#grabs stuff and returns a list of them. This is used throughout scripts.
def grab_ref(ref_folder1):
    files=glob.glob(ref_folder1+'/*')
    return(files)

"""
Download data using fastq-dump, quantifies data using Salmon, and deletes 
downloaded data. 

This is the heart of the project. So ugly. 
"""

#dummy variables until command line input is setup
w = input_lis[0]

for y in glob.glob('INDEX_'+ref_folder+'/*'):
    index_name = y.replace('INDEX_references/','').replace('.ffn.gz','')
    print('Index: '+index_name)
    for q in input_lis[1:]:
        dump_directory = 'data/'+w+'/'+q
        os.system('mkdir data/'+w+'/'+q)
        print('Dump Directory: '+dump_directory)
        os.system('fastq-dump --outdir '+dump_directory+' -v --gzip --skip-technical --readids --split-3 --clip '+q)
    run_folders = grab_ref('data/'+w)
    output_name='quants/'+index_name+'/'+w+'_quant'
    print('Quant output:'+output_name)
    if len(run_folders)<1:
        print('Oops, no runs for '+w)
    if len(run_folders)==1:
        if len(grab_ref(run_folders[0]))<1:
            print('Oops, no runs for '+w)
        elif len(grab_ref(run_folders[0]))==1:
            print('One run for ')
            print('Salmon input: '+grab_ref(run_folders[0])[0])
            os.system('salmon quant -i '+y+' -l A '+grab_ref(run_folders[0])[0]+' -p 8 --validateMappings -o '+output_name)
        elif len(grab_ref(run_folders[0]))>1:
            runs = grab_ref(run_folders[0])
            print('Multiple runs for ')
            if "_1" in runs[0] or "_2" in runs[0]:
                print('Multiple paired runs for '+w)
                input_name='-1'
                for z in np.arange(0, len(runs)):
                    if "_1" in runs[z]:
                        input_name+= ' '+runs[z]
                input_name+=' -2'
                for z in np.arange(0, len(runs)):
                    if "_2" in runs[z]:
                        input_name+= ' '+runs[z]
            else:
                print('Multiple runs for '+z)
                input_name='-r'
                for z in np.arange(0, len(runs)):
                    input_name+= ' '+runs[z]
            print('Salmon Input: '+input_name)
    if len(run_folders)>1:
        print('Multiple Folder Case')
        first_read_str = ''
        second_read_str = ''
        orphan_read = '-r'
        for run_folder in run_folders:
            if len(grab_ref(run_folder))<1:
                print('Opps, no runs for '+run_folder)
            elif len(grab_ref(run_folder))==1:
                print('One run for '+run_folder)
                print('Salmon input: '+grab_ref(run_folder)[0])
                orphan_read += ' '+grab_ref(run_folder)[0]
                os.system('salmon quant -i '+y+' -l A '+grab_ref(run_folders[0])[0]+' -p 8 --validateMappings -o '+output_name)
            elif len(grab_ref(run_folder))>1:
                print('Multiple runs for '+run_folder)
                runs=grab_ref(run_folder)
                if "_1" in runs[0] or "_2" in runs[0]:
                    print('Multiple paired runs for '+w)
                    for z in np.arange(0, len(runs)):
                        if "_1" in runs[z]:
                            first_read_str += ' '+runs[z]
                    for z in np.arange(0, len(runs)):
                        if "_2" in runs[z]:
                            second_read_str+= ' '+runs[z]
                else:
                    print('Multiple runs for '+run_folder)
                    for z in np.arange(0, len(runs)):
                        orphan_read+=runs[z]
        if len(orphan_read) > 2: 
            if len(first_read_str)>1 or len(second_read_str)>1:
                input_name = orphan_read+first_read_str+second_read_str
            else: 
                input_name = orphan_read
        elif len(orphan_read) == 2:
            input_name = '-1'+first_read_str+' -2'+second_read_str
        print('Salmon Input: '+input_name)
        os.system('salmon quant -i '+y+' -l A '+input_name+' -p 8 --validateMappings -o '+output_name)
        os.system('rm -r data/'+w)