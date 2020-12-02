# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 12:54:57 2020

@author: holtj
"""

import pandas as pd
from subprocess import Popen, PIPE
import time

"""
Create dictionary
"""
#read in SRR#s from metadata table provided by Georgia
sra_run_table = pd.read_csv(r'SraRunTable.csv', sep=",", header=0, index_col=0)
#create a dict key of biosamples : srr numbers(could do SRX:SRR instead)
bio_sample_dic = sra_run_table.groupby('BioSample').agg({'Run':lambda x:x.tolist()})['Run'].to_dict()

"""
Start submitting jobs that call quantifier script 
Needs inexer job to be finished. 
"""
for i in bio_sample_dic.keys()[600:601]:
    # Open a pipe to the qsub command.
    proc = Popen('qsub -hold_jid my_job_indexer', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
    
    # Customize your options here
    job_name = "my_job_%d" % i
    walltime = "1:00:00"
    processors = "nodes=1:ppn=1"
    command = "python quantifier.py -l "+i+","+bio_sample_dic[i]

    job_string = """#!/bin/bash
    #PBS -q tyestq
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
    
    module load sratoolkit
    
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

    