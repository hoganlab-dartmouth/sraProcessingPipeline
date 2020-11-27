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
job_name = "my_job_indexer"
walltime = "1:00:00"
processors = "nodes=1:ppn=1"
command = "python ~/sraProcessingPipeline/indexer/indexer.py"

job_string = """#!/bin/bash
#PBS -q testq
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