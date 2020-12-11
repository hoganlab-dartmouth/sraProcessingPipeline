# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:28:34 2020

@author: holtj

code for creating graphs
"""

"""
modules to be used
"""
import pandas as pd 
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
file_name = "PA14_K31_R2"
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

#this is really ugly, but we're gonna rename columns and make three new dataframes
clc = pd.read_csv(r'clc.csv')
clc = clc.set_index('genes')
df2=pd.DataFrame()
df2[['clc_a','clc_b','clc_wt','clc_wt_b']]=clc[['a','b','wt','wt_b']]
df1=pd.DataFrame()
df1[['s_a','s_b','s_wt', 's_wt_b']]=salmon[['a','b','wt','wt_b']]

#df containing salmon+clc
df = pd.concat([df2, df1], axis=1)
df = df.dropna()

#%%%
"""
Quick scatter plot of Salmon vs CLC indexed by PA14 gene id
"""
plt.scatter(df['clc_a'],df['s_a'], alpha=0.5, label='mut_rep1')
plt.scatter(df['clc_b'],df['s_b'], alpha=0.5, label='mut_rep2')
plt.scatter(df['clc_wt'],df['s_wt'], alpha=0.5, label='wt_rep1')
plt.scatter(df['clc_wt_b'],df['s_wt_b'], alpha=0.5, label='wt_rep2')
x = np.linspace(0,5)
y = x
plt.plot(x,y, color='red', linestyle='dashed')



plt.xlabel('CLC Log10(TPM)')
plt.ylabel('Salmon Log10(TPM)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title(file_name+', No Bias Correction')
plt.show()

#%%%
"""
Quick scatter plot comparing forward and reverse reads
"""

#load in reverse reads, reverse reads are R1 in the tester_data
file_name2 = "PA14_K31_Unpaired"
salmon_r = pd.read_csv (file_name+'.csv')

#rename genes as PA14 and make it so row names are column names
seq_id_to_gene_id_pa14 = dict_gene_num_to_ids("PA14.gz")
salmon_r.rename(mapper=seq_id_to_gene_id_pa14, axis="columns", inplace=True)
salmon_r = salmon_r.transpose()
new_header = salmon_r.iloc[0]
salmon_r = salmon_r[1:] 
salmon_r.columns_r = new_header
salmon_r.columns=['a','wt','b','wt_b']



#add 1 to all points in dataset before taking log10
salmon_r = salmon_r+1
for i in list(salmon_r.columns.values):
    salmon_r[i] = pd.to_numeric(salmon_r[i], downcast="float")
salmon_r = np.log10(salmon_r)

#this is really ugly but we're gonna rename columns by making new dataframes
salmon_r2=pd.DataFrame()
salmon_r2[['mut_rep1_r','mut_rep2_r','wt_rep1_r','wt_rep2_r']]=salmon_r[['a','b','wt','wt_b']]
salmon_f=pd.DataFrame()
salmon_f[['mut_rep1_f','mut_rep2_f','wt_rep1_f','wt_rep2_f']]=salmon[['a','b','wt','wt_b']]

f_and_r = pd.concat([salmon_r2, salmon_f], axis=1)

#make the scatters
for i in np.arange(0,len(list(salmon_r2.columns.values))):
    plt.scatter(f_and_r[list(salmon_r2.columns.values)[i]],f_and_r[list(salmon_f.columns.values)[i]], alpha=0.5, label=list(salmon_r2.columns.values)[i][:-2])
x = np.linspace(0,5)
y = x
plt.plot(x,y, color='red', linestyle='dashed')



plt.xlabel(file_name2+', Log10(TPM)')
plt.ylabel(file_name+', Log10(TPM)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('Salmon, No Bias Correction')
plt.show()
#%%
"""
Residual Squared plot with RMSD
Not sure if RMSD is approprtiate or that I am doing it correctly? 
"""
res_sum = 0
plt.scatter(df['clc_a'],(df['s_a']-df['clc_a'])**2, alpha=0.5, label='mut_rep1')
res_sum += ((df['clc_a']-df['s_a'])**2).sum()
plt.scatter(df['clc_b'],(df['s_b']-df['clc_b'])**2, alpha=0.5, label='mut_rep2')
res_sum += ((df['clc_b']-df['s_b'])**2).sum()
plt.scatter(df['clc_wt'],(df['s_wt']-df['clc_wt'])**2, alpha=0.5, label='mut_rep1')
res_sum += ((df['clc_wt']-df['s_wt'])**2).sum()
plt.scatter(df['clc_wt_b'],(df['s_wt_b']-df['clc_wt_b'])**2, alpha=0.5, label='mut_rep2')
res_sum += ((df['clc_wt_b']-df['s_wt_b'])**2).sum()
plt.xlabel('CLC')
plt.ylabel('Residual Squared')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title(file_name+', No Bias Correction: RMSD='+str(round(((res_sum/(len(df-2)))**(1/2))/4,2))+', N='+str(len(df)))
plt.show()
#%%
"""
Here we are pulling up the genes that are highly expressed in Salmon, but not so in CLC. 

Once we have those genes we rename them (if we can)
"""

#create a dictionary that lets us rename genes
gene_key = pd.read_csv(r'PA14_Keys.csv')
new_header = gene_key.iloc[1] #grab the first row for the header
gene_key = gene_key[2:] #take the data less the header row
gene_key.columns = new_header #set the header row as the df header
gene_key = gene_key.reset_index()
gene_key['LocusTag'] = gene_key['LocusTag'].str.replace(r'"', '')
gene_keys = dict(zip(gene_key.LocusTag,gene_key.Name))
del_list = []
for i in gene_keys.keys():
    if type(gene_keys[i]) != str:
        del_list.append(i)
for z in del_list:
    del gene_keys[z]
    

discrep_output_name = file_name+"_discrep_output"
discrep_frame = pd.DataFrame()
header_lis = []
header_lis2 = []
heads = list(salmon.columns.values)

for header in heads:
    tpm_lis = []
    gene_lis = []
    for index in df.index.tolist():
        if round(clc[header][index], 2) < 0.5 and round(salmon[header][index], 2) > 1.0 :
            tpm_lis.append(10**salmon[header][index])
            if index not in gene_keys.keys():
                gene_lis.append(index)
            else:
                gene_lis.append(gene_keys[index])
    header_lis.append(tpm_lis)
    header_lis2.append(gene_lis)
longest_list = max(len(header) for header in header_lis)

for i in np.arange(0,len(header_lis)):
    if len(header_lis[i])==longest_list:
        discrep_frame[heads[i]]=pd.Series(header_lis[i])
        discrep_frame[heads[i]+'_genes']=pd.Series(header_lis2[i])

for i in np.arange(0,len(header_lis)):
    if len(header_lis[i])==longest_list:
        continue        
    else:
        discrep_frame[heads[i]]=pd.Series(header_lis[i]).fillna(0)
        discrep_frame[heads[i]+'_genes']=pd.Series(header_lis2[i]).fillna(0)

discrep_frame.to_csv(discrep_output_name+'.csv', header=True)

#%%
"""
Graph the highest 25 expressed genes

I should redo these so that the shape reflects the alignment method. 
"""
metric ='clc_wt'
df = df.sort_values(by=metric, ascending=False).dropna()
plt.scatter(df.index[0:25], df['s_wt'][0:25], alpha=0.5, label='Salmon_wt_rep1')
plt.scatter(df.index[0:25], df['clc_wt'][0:25], alpha=0.5, label='CLC_wt_rep1')
plt.scatter(df.index[0:25], df['s_wt_b'][0:25], alpha=0.5, label='Salmon_wt_rep2')
plt.scatter(df.index[0:25], df['clc_wt_b'][0:25], alpha=0.5, label='CLC_wt_rep2')
plt.scatter(df.index[0:25], df['s_a'][0:25], alpha=0.5, label='Salmon_mut_rep1')
plt.scatter(df.index[0:25], df['clc_a'][0:25], alpha=0.5, label='CLC_mut_rep1')
plt.scatter(df.index[0:25], df['s_b'][0:25], alpha=0.5, label='Salmon_mut_rep2')
plt.scatter(df.index[0:25], df['clc_b'][0:25], alpha=0.5, label='CLC_mut_rep2')
plt.xticks(rotation=90)
plt.xlabel('Genes')
plt.ylabel('Log10(TPM)')
plt.title(file_name+', No Bias Correction, Sorted by: '+metric)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

#%%
"""
Graph the lowest 25 expressed genes
"""
metric ='clc_wt'
df = df.sort_values(by=metric, ascending=True).dropna()
plt.scatter(df.index[0:25], df['s_wt'][0:25], alpha=0.5, label='Salmon_wt_rep1')
plt.scatter(df.index[0:25], df['clc_wt'][0:25], alpha=0.5, label='CLC_wt_rep1')
plt.scatter(df.index[0:25], df['s_wt_b'][0:25], alpha=0.5, label='Salmon_wt_rep2')
plt.scatter(df.index[0:25], df['clc_wt_b'][0:25], alpha=0.5, label='CLC_wt_rep2')
plt.scatter(df.index[0:25], df['s_a'][0:25], alpha=0.5, label='Salmon_mut_rep1')
plt.scatter(df.index[0:25], df['clc_a'][0:25], alpha=0.5, label='CLC_mut_rep1')
plt.scatter(df.index[0:25], df['s_b'][0:25], alpha=0.5, label='Salmon_mut_rep2')
plt.scatter(df.index[0:25], df['clc_b'][0:25], alpha=0.5, label='CLC_mut_rep2')
plt.xticks(rotation=90)
plt.xlabel('Genes')
plt.ylabel('Log10(TPM)')
plt.title(file_name+', No Bias Correction, Sorted by: '+metric)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

#%%
"""
Graph highest expressed genes of interest. 
"""
genes = (pd.read_csv(r'H_pho_reg_PA14_2.csv', names=['gene']))
genes=list(genes['gene'])
s_df = df[df.index.isin(genes)]

metric ='clc_wt'

df = df.sort_values(by=metric, ascending=False).dropna()
plt.scatter(s_df.index[0:25], s_df['s_wt'][0:25], alpha=0.5, label='Salmon_wt_rep1')
plt.scatter(s_df.index[0:25], s_df['clc_wt'][0:25], alpha=0.5, label='CLC_wt_rep1')
plt.scatter(s_df.index[0:25], s_df['s_wt_b'][0:25], alpha=0.5, label='Salmon_wt_rep2')
plt.scatter(s_df.index[0:25], s_df['clc_wt_b'][0:25], alpha=0.5, label='CLC_wt_rep2')
plt.scatter(s_df.index[0:25], s_df['s_a'][0:25], alpha=0.5, label='Salmon_mut_rep1')
plt.scatter(s_df.index[0:25], s_df['clc_a'][0:25], alpha=0.5, label='CLC_mut_rep1')
plt.scatter(s_df.index[0:25], s_df['s_b'][0:25], alpha=0.5, label='Salmon_mut_rep2')
plt.scatter(s_df.index[0:25], s_df['clc_b'][0:25], alpha=0.5, label='CLC_mut_rep2')
plt.xticks(rotation=90)
plt.xlabel('Genes')
plt.ylabel('Log10(TPM)')
plt.title(file_name+', No Bias Correction, Pho Genes sorted by: '+metric)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()
    
#%%
"""
Make a list of the genes that we are interested in and use it to make a new dataframe interest_df
"""
genes = (pd.read_csv(r'H_pho_reg_PA14_2.csv', names=['gene']))
genes=list(genes['gene'])
interest_df = df[df.index.isin(genes)]
#interest_df.index.replace(gene_keys, inplace=True)

#%%
"""
Heatmap of genes of interest (see above)
"""
fg = sns.clustermap(interest_df.dropna())
plt.show()

#%%
for i in list(df.columns.values):
    plt.hist(df[i], bins = freedman_diaconis(df[i]), alpha=0.5, label=i)
plt.legend()
plt.xlabel('Log10(TPM)')
plt.ylabel('counts')
plt.title(file_name+", No Bias Correction")
plt.show()

#%%
for i in list(df.columns.values):
    plt.scatter(np.arange(0, len(df)),df[i], alpha=0.5, label=i)
plt.legend(title='Samples', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('Gene #')
plt.ylabel('Log10(TPM)')
plt.title(file_name+", No Bias Correction")
plt.show()

