# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 21:05:03 2020

@author: holtj

To test out python submitting jobs.
"""

from subprocess import Popen, PIPE
import time

"""
Submit indexer job
This call creates the Salmon index files. 
"""
# Open a pipe to the qsub command.
proc = Popen('qsub', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

# Customize your options here
job_name = "indexer_test"
walltime = "1:00:00"
processors = "nodes=1:ppn=1"
command = "python indexer.py"

job_string ="""#!/bin/bash -l
#PBS -q default
#PBS -N %s 
#PBS -l walltime=%s
#PBS -l %s
#PBS -l mem=2gb
#PBS -m bea
#PBS -M jacob.d.holt.gr@dartmouth.edu
#PBS -l feature='cellk'

cd "/dartfs-hpc/rc/home/5/f004ky5/sraProcessingPipeline"

module load python/anaconda3

module load salmon

%s

exit 0

""" % (job_name, walltime, processors, command)

# Send job_string to qsub

proc.stdin.write(job_string.encode('utf-8'))

print(proc.communicate())