# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:44:29 2020

@author: holtj

This executable is the pipeline management. It submits the jobs for each script
that we want to call. There are three python scripts: indexer, quantifier, and csv_builder

I am not sure how to make sure the jobs happen in the right order. 

As of 11/24, the best way to run this is to do the indexing locally and setup 
the directories with the right files on Discovery. Then you can start writing an submitting jobs using this run file.
Once your jobs are done, you then move the Ex folder to your local machine and run the csv_builder.  Eventually I'd like
for everything to be run on discovery, but for testing this should be fine. 

Home directory: 
DartFS: /dartfs-hpc/rc/home/e/d12345e (50GB)

/scratch is local to nodes. Job data can be stored here but should be deleted. 
"""

import pandas as pd
from subprocess import Popen, PIPE
import time

"""
Make sure Salmon and SRAToolkit are installed and configured properly
Copied Alex's code, with some modifications, for this.
Not really sure if it works yet. 

Will probably end up just doing this manually before calling scripts. 
"""
# Download latest version of compiled binaries of NCBI SRA toolkit 
# if not os.path.exists("sratoolkit.current-centos_linux64.tar.gz"):
#     wget.download("ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz")


# Extract tar.gz file 
# if os.path.exists("sratoolkit.current-centos_linux64.tar.gz"):
#     os.system('tar -xzf sratoolkit.current-centos_linux64.tar.gz')
#     os.system('export PATH=$PATH:sratoolkit.2.10.7-centos_linux64/bin')

# Now SRA binaries added to path and ready to use

"""
Setup directories
These really should be input options
"""
#ref folder name, ref genomes go here
ref_folder = '/dartfs-hpc/scratch/Jake/references'
#os.system('mkdir '+ref_folder)
#data folder name, samples and runs go here
data = '/dartfs-hpc/scratch/Jake/data'
#os.system('mkdir '+data)
#csv output folder name
csv = '/dartfs-hpc/scratch/Jake/Ex'
#os.system('mkdir '+csv)

"""
Create dictionary
"""
#read in SRR#s from metadata table provided by Georgia
sra_run_table = pd.read_csv(r'SraRunTable.csv', sep=",", header=0, index_col=0)
#create a dict key of biosamples : srr numbers(could do SRX:SRR instead)
bio_sample_dic = sra_run_table.groupby('BioSample').agg({'Run':lambda x:x.tolist()})['Run'].to_dict()

"""
Submit indexer job
This call creates the Salmon index files. 
"""
# Open a pipe to the qsub command.
proc = Popen('qsub', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

# Customize your options here
job_name = "indexer"
walltime = "1:00:00"
processors = "nodes=1:ppn=1"
command = "python indexer.py"

job_string = """#!/bin/bash
#PBS -q largeq
#PBS -A NCCCsub 
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -l feature='cellf'
#PBS -m bea
#PBS -M jacob.d.holt.gr@dartmouth.edu
#PBS -o ./output/%s.out
#PBS -e ./error/%s.err

cd $PBS_O_WORKDIR

module load Salmon
%s""" % (job_name, walltime, processors, job_name, job_name, command)

# Send job_string to qsub
#if (sys.version_info > (3, 0)):
proc.stdin.write(job_string.encode('utf-8'))

#proc.stdin.write(job_string)
out, err = proc.communicate()

# Print your job and the system response to the screen as it's submitted
print(job_string)
print(out)

time.sleep(0.1)
    
"""
Start submitting jobs that call quantifier script once the indexer job is done. 
"""
for i in list(bio_sample_dic.keys()):
    # Open a pipe to the qsub command.
    proc = Popen('qsub -hold_jid "indexer*"', shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
    
    # Customize your options here
    job_name = "single_quant_%s" % bio_sample_dic[i]
    walltime = "01:00:00"
    processors = "nodes=1:ppn=1"
    command = "python single_quantifier.py -l "+i+","+bio_sample_dic[i]

    job_string = """#!/bin/bash
    #PBS -q largeq
    #PBS -N %s 
    #PBS -l walltime=%s
    #PBS -l %s
    #PBS -l mem=2gb
    
    cd "/dartfs-hpc/rc/home/5/f004ky5/sraProcessingPipeline"

    module load python/anaconda3

    module load salmon

    module load sratoolkit
    
    %s
    
    """ % (job_name, walltime, processors, command)
    
    # Send job_string to qsub
    #if (sys.version_info > (3, 0)):
    proc.stdin.write(job_string.encode('utf-8'))
    
    #proc.stdin.write(job_string)
    out, err = proc.communicate()
    
    # Print your job and the system response to the screen as it's submitted
    print(job_string)
    print(out)
    
    time.sleep(1)
    
"""
CSV Builder Job submit
"""
# Open a pipe to the qsub command.
proc = Popen('qsub -hold_jid "quant*"', shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

# Customize your options here
job_name = "csv_builder"
walltime = "00:20:00"
processors = "nodes=1:ppn=1"
command = "csv_builder.py "

job_string = """#!/bin/bash
#PBS -q largeq
#PBS -N %s 
#PBS -l walltime=%s
#PBS -l %s
#PBS -l mem=2gb

cd "/dartfs-hpc/rc/home/5/f004ky5/sraProcessingPipeline"

module load python/anaconda3

module load salmon

module load sratoolkit

%s

""" % (job_name, walltime, processors, command)

# Send job_string to qsub
#if (sys.version_info > (3, 0)):
proc.stdin.write(job_string.encode('utf-8'))

#proc.stdin.write(job_string)
out, err = proc.communicate()

# Print your job and the system response to the screen as it's submitted
print(job_string)
print(out)

time.sleep(0.1)