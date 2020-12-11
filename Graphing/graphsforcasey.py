# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 07:48:38 2020

@author: holtj
"""

"""
modules to be used
"""
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import scipy.stats as stats
import math
import seaborn as sns
from sklearn.linear_model import LinearRegression
from bioinfokit import analys, visuz
import random
import gzip
from Bio import SeqIO
from collections import OrderedDict

"""
Functions to be used. First one is copied from Alex. 
"""

def dict_gene_num_to_ids(fasta_file):
    """
    Returns mapping between gene number created from SRA (i.e. PGD######)
    and gene locus (PA######). The PA###### is the position of the gene on
    the Pseudomonas genome and is more interpretable compared to the
    SRA-generated PGD######
    Arguments
    ----------
    fasta_file: str
        Reference transcriptome file 
    """

    seq_id_to_gene_id = {}

    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id
            seq_record_desc = record.description
            seq_record_desc_lst = seq_record_desc.split(";")
            tagged_item = [
                item.split("=")[1]
                for item in seq_record_desc_lst
                if "locus_tag" in item
            ][0]

            seq_id_to_gene_id[seq_id] = tagged_item

    return seq_id_to_gene_id

def freedman_diaconis(data, returnas="bins"):
    """
    Use Freedman Diaconis rule to compute optimal histogram bin width. 
    ``returnas`` can be one of "width" or "bins", indicating whether
    the bin width or number of bins should be returned respectively. 


    Parameters
    ----------
    data: np.ndarray
        One-dimensional array.

    returnas: {"width", "bins"}
        If "width", return the estimated width for each histogram bin. 
        If "bins", return the number of bins suggested by rule.
    """
    data = np.asarray(data, dtype=np.float_)
    IQR  = stats.iqr(data, rng=(25, 75), scale="raw", nan_policy="omit")
    N    = data.size
    bw   = (2 * IQR) / np.power(N, 1/3)

    if returnas=="width":
        result = bw
    else:
        datmin, datmax = data.min(), data.max()
        datrng = datmax - datmin
        result = int((datrng / bw) + 1)
    return(result)

sns.color_palette()

#%%%
"""
Load in data from csv file of interest (see files in folder) and some processing. 
Note: tpp.csv is the data from CLC
All other csv files are from Salmon
Files are using TPM, but we are going to change to log10(TPM) to compare across samples. 
"""
#make data frame from csv
file_name = "PA14_K31_Paired_Bias_Acc"
file_name_2 = "PA14_K31_R2_R1_Reversed"
salmon = pd.read_csv (file_name+'.csv')

#rename genes as PA14 and make it so row names are column names
seq_id_to_gene_id_pa14 = dict_gene_num_to_ids("PA14.gz")
salmon.rename(mapper=seq_id_to_gene_id_pa14, axis="columns", inplace=True)
salmon = salmon.transpose()
new_header = salmon.iloc[0] 
salmon = salmon[1:] 
salmon.columns = new_header 

#add 1 to all points in dataset before taking log10
salmon = salmon+1
for i in ['a','b', 'wt', 'wt_b']:
    salmon[i] = pd.to_numeric(salmon[i], downcast="float")
salmon = np.log10(salmon)

#make data frame from csv
salmon2 = pd.read_csv (file_name_2+'.csv')

#rename genes as PA14 and make it so row names are column names
salmon2.rename(mapper=seq_id_to_gene_id_pa14, axis="columns", inplace=True)
salmon2 = salmon2.transpose()
new_header = salmon2.iloc[0] 
salmon2 = salmon2[1:] 
salmon2.columns = new_header 

#add 1 to all points in dataset before taking log10
salmon2 = salmon2+1
for i in ['a','b', 'wt', 'wt_b']:
    salmon2[i] = pd.to_numeric(salmon2[i], downcast="float")
salmon2 = np.log10(salmon2)

#this is really ugly, but we're gonna rename columns and make three new dataframes
clc = pd.read_csv(r'clc.csv')
clc = clc.set_index('genes')
df2=pd.DataFrame()
df2[['clc_a','clc_b','clc_wt','clc_wt_b']]=clc[['a','b','wt','wt_b']]
df1=pd.DataFrame()
df1[['s_a_1','s_b_1','s_wt_1', 's_wt_b_1']]=salmon[['a','b','wt','wt_b']]
df3=pd.DataFrame()
df3[['s_a_2','s_b_2','s_wt_2', 's_wt_b_2']]=salmon2[['a','b','wt','wt_b']]

#df containing salmon+clc
df = pd.concat([df2, df1, df3], axis=1)
df = df.dropna()

#%%
"""
Quick scatter plot colored by salmon or salmon2
"""
file_name3= "Paired\n GC Bias Corr"
file_name_4 = "Unpaired\n No GC Bias Corr.\n Reverse Complement"
#df = df[df > 0]
plt.scatter(df['clc_a'],df['s_a_1'], alpha=0.5, label=file_name3, color='k')
plt.scatter(df['clc_b'],df['s_b_1'], alpha=0.5, label=file_name3, color='k')
plt.scatter(df['clc_wt'],df['s_wt_1'], alpha=0.5, label=file_name3, color='k')
plt.scatter(df['clc_wt_b'],df['s_wt_b_1'], alpha=0.5, label=file_name3, color='k')
plt.scatter(df['clc_a'],df['s_a_2'], alpha=0.5, label=file_name_4, color='r')
plt.scatter(df['clc_b'],df['s_b_2'], alpha=0.5, label=file_name_4, color='r')
plt.scatter(df['clc_wt'],df['s_wt_2'], alpha=0.5, label=file_name_4, color='r')
plt.scatter(df['clc_wt_b'],df['s_wt_b_2'], alpha=0.5, label=file_name_4, color='r')
x = np.linspace(0,5)
y = x
plt.plot(x,y, color='k', linestyle='dashed')



plt.xlabel('CLC Log10(TPM)')
plt.ylabel('Salmon, Log10(TPM)')
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(),bbox_to_anchor=(1.05, 1), loc='upper left')
#plt.title(file_name+', No Bias Correction')
plt.show()
#%%
"""
Quick scatter plot colored by salmon or salmon2
"""
file_name3= "Uncorrected"
file_name_4 = "Polycistronic Corrected"
#df = df[df > 0]
# plt.scatter(df['s_a_1'], df['clc_a'], alpha=0.5, label=file_name3, color='k')
# plt.scatter(df['s_b_1'], df['clc_b'], alpha=0.5, label=file_name3, color='k')
# plt.scatter(df['s_wt_1'], df['clc_wt'], alpha=0.5, label=file_name3, color='k')
# plt.scatter(df['s_wt_b_1'], df['clc_wt_b'], alpha=0.5, label=file_name3, color='k')
plt.scatter(df['s_a_2'], df['clc_a'], alpha=0.5, label=file_name_4, color='r')
plt.scatter(df['s_b_2'], df['clc_b'], alpha=0.5, label=file_name_4, color='r')
plt.scatter(df['s_wt_2'], df['clc_wt'], alpha=0.5, label=file_name_4, color='r')
plt.scatter(df['s_wt_b_2'], df['clc_wt_b'], alpha=0.5, label=file_name_4, color='r')
x = np.linspace(0,5)
y = x
plt.plot(x,y, color='k', linestyle='dashed')

plt.ylabel('CLC, '+r'$\mathtt{log}_{10}(\mathtt{TPM})$') 
plt.xlabel('Salmon, '+r'$\mathtt{log}_{10}(\mathtt{TPM})$')
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(),loc='upper left')
#plt.title(file_name+', No Bias Correction')
plt.show()