# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 14:37:19 2020

@author: holtj

Python script for differential gene expression analysis. 
"""

import pandas as pd
import numpy as  np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
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
    
result = {}

for key,value in gene_keys.items():
    if value not in result.values():
        result[key] = value

count_file = "PAO1"
metadata_file = "SraRunTable"

# tell pandas to make a new DataFrame with the contents of `brain_counts.csv`. This might take a minute.
count_dataframe = pd.read_csv(count_file+'.csv', index_col=0).sort_index(axis=0)
metadata_dataframe = pd.read_csv(metadata_file+'.csv')
metadata_dataframe["Treatment"]=metadata_dataframe['BioSample']

seq_id_to_gene_id_pa14 = dict_gene_num_to_ids("PAO1.gz")
count_dataframe.rename(mapper=seq_id_to_gene_id_pa14, axis="columns", inplace=True)

#drop emptys sets from TPM data
count_dataframe = count_dataframe[(count_dataframe.T != 0).any()]
count_dataframe = np.log10(count_dataframe + 1)

#Set index as SRR number for both data frames. Make sure they have the same index. 

metadata_dataframe = metadata_dataframe[metadata_dataframe["Index"].isin(count_dataframe.index)]
metadata_dataframe = metadata_dataframe.set_index("Index").sort_index(axis=0)
metadata_dataframe.to_csv("Subset.csv", header=True)

for name in list(count_dataframe.columns.values):
    if name in result:
        count_dataframe.rename(columns={name:result[name]}, inplace=True)

metadata_dataframe['Study'] = metadata_dataframe['BioSample'].map(lambda x: str(x)[:-1])

# GeorgiaList = metadata_dataframe[~metadata_dataframe["SAMPLE_TYPE"].str.contains('clinical', na=False)].sort_index(axis=0)
# count_dataframe = count_dataframe[count_dataframe.index.isin(GeorgiaList.index)].sort_index(axis=0)

adata = sc.AnnData(X = count_dataframe, obs = metadata_dataframe)

#%%
# """
# PCA
# Normalize was done earlier by doing log10(TPM+1)
# """
# sc.pp.pca(adata)
# sc.pl.pca_overview(adata, color='Study')

#%%
"""
Dimensionality Reduction
"""
# sc.tl.tsne(adata, perplexity=30, learning_rate=1000, random_state=0)
# sc.pl.tsne(adata, color='souce_name')

#%%

#sc.pl.highest_expr_genes(adata, n_top=20)

cake_lis = ['lasR', 'phzM', 'phoB','source_name']
gene_lis = ['PA0143']

sc.pp.neighbors(adata) # UMAP is based on the neighbor graph; we'll compute this first
sc.tl.umap(adata, min_dist=0.5, spread=1.0, random_state=1, n_components=2)

sc.pl.umap(adata, color=gene_lis)

# #%%
# """
# Clustering
# """
# umap_coordinates = adata.obsm['X_umap'] # extract the UMAP coordinates for each cell
# kmeans = KMeans(n_clusters=4, random_state=0).fit(umap_coordinates) # fix the random state for reproducibility

# adata.obs['kmeans'] = kmeans.labels_ # retrieve the labels and add them as a metadata column in our AnnData object
# adata.obs['kmeans'] = adata.obs['kmeans'].astype(str)

# sc.pl.umap(adata, color='kmeans') # plot the results

# # rand_index = adjusted_rand_score(labels_true = adata.obs['source_name'], labels_pred = adata.obs['kmeans'])
# # print('The Rand index is', round(rand_index, 2))

# # sc.tl.louvain(adata, resolution=0.1)

# #%%
# """
# Differential gene expression analysis
# """
# sc.tl.rank_genes_groups(adata, groupby='infection_side', use_raw=True, 
#                         method='t-test_overestim_var', n_genes=10) # compute differential expression
# sc.pl.rank_genes_groups_tracksplot(adata, groupby='Study') # plot the result
