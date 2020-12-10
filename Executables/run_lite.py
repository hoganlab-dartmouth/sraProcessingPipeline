# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 12:54:57 2020

@author: holtj

I was playing around with using afterok here to get around the limits on how 
many jobs I could submit at one time. I think because Georgia is attached to a group
she will not have this same issue. If you doing more than a couple hundred samples, it's probably
a good idea to contactt he discovery people though to let them know what you're doing and 
make sure they don't lock your account.
"""

import csv
from subprocess import Popen, PIPE
import time

"""
Create dictionary
"""
with open('GeorgiaTesterRemainder.csv', mode='r') as infile:
    reader = csv.reader(infile)
    bio_sample_dic = {rows[0]:rows[1] for rows in reader}
"""
Start submitting jobs that call quantifier script 
Needs inexer job to be finished first. 
Runs csv_builder.py after

I did this during debugging to make things easier. The indexer and csv_builder
 were pretty easy to get up and 
"""
previous_job = ""
for i in list(bio_sample_dic.keys())[0:1]:
    # Open a pipe to the qsub command.
    proc = Popen('qsub', shell=False, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
    
    # Customize your options here
    job_name = "single_quant_%s" % bio_sample_dic[i]
    walltime = "00:05:00"
    processors = "nodes=1:ppn=1"
    command = "python single_quantifier.py -l "+i+","+bio_sample_dic[i]

    job_string = """#!/bin/bash
    #PBS -q testq
    #PBS -N %s 
    #PBS -l walltime=%s
    #PBS -l %s
    #PBS -l mem=2gb
    #PBS W depend=afterok:%s
    
    cd "/dartfs-hpc/rc/home/5/f004ky5/sraProcessingPipeline"

    module load python/anaconda3

    module load salmon

    module load sratoolkit
    
    %s
    
    """ % (job_name, walltime, processors, previous_job, command)
    
    # Send job_string to qsub
    #if (sys.version_info > (3, 0)):
    proc.stdin.write(job_string.encode('utf-8'))
    
    #proc.stdin.write(job_string)
    out, err = proc.communicate()
    
    # Print your job and the system response to the screen as it's submitted
    print(job_string)
    print(out)
    previous_job = out.decode('utf-8')[0:7]
    print("job id:"+str(out[0:7])[2:8])
    
    time.sleep(1)
    