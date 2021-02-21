import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD

import os
import gc

celltype = [
    'mega', 'Mono-CD14-HLA', 'Mono-CD14-S100A8-RETN', 'Mono-CD14-S100A8-CD163',
    'Mono-CD14-CCL3', 'Neu', 'Mono−CD16', 'T-CD4-FOS', 'T-CD4-LTB-S100A4',
    'T-CD8-SLC4A10'
]

pbmc_anndata = sc.read('/home/liunianping/workspace/projects/Infection_Virus/data/Integration/pbmc.virus.refine.6.h5ad')
pbmc_anndata = pbmc_anndata[(pbmc_anndata.obs['source'] == 'Covid19') |
                  (pbmc_anndata.obs['source'] == 'Control')]
pbmc_anndata = pbmc_anndata[pbmc_anndata.obs['dataset'] != 'd01']
pbmc_anndata = pbmc_anndata[pbmc_anndata.obs['dataset'] != 'd12']
pbmc_anndata = pbmc_anndata[pbmc_anndata.obs['dataset'] != 'd14']


pbmc_anndata = pbmc_anndata[pbmc_anndata.obs['louvain_celltype'].apply(
        lambda x: x in celltype).astype(np.bool)]

X = pbmc_anndata.X.todense()

group = pbmc_anndata.obs[['sampleID', 'dataset', 'state', 'inflammatory_cytokine']].groupby('sampleID')

ann = pd.DataFrame()
ann['sampleID'] = list(group.first().index)
ann['dataset'] = list(group.first()['dataset'])
ann['state'] = list(group.first()['state'])
ann['inflammatory_cytokine'] = list(group['inflammatory_cytokine'].mean())

data = pd.DataFrame(X, index=pbmc_anndata.obs['sampleID'], columns=pbmc_anndata.var.index)
data = data.groupby('sampleID').mean()
data = sc.AnnData(data)
data.var.index = list(pbmc_anndata.var.index)
data.var['gene'] = pbmc_anndata.var['gene_name-covid19']
data.obs['dataset'] = list(group.first()['dataset'])

geneLists = []
for i in data.obs['dataset'].unique():
    temp = data[data.obs['dataset'] == i]
    sc.pp.highly_variable_genes(temp)
    genes = temp[:,temp.var['highly_variable']].var.index.values
    geneLists.extend(genes)

from collections import Counter
c = Counter(geneLists)

genelist = []
for i in c.most_common(2000):
    genelist.append(i[0])

temp = data[:,genelist]
pd.DataFrame(np.exp(temp.X)-1, index=list(ann['sampleID']), columns=temp.var.index).to_csv('cluster10_hvg_overlap2000_exp.csv')
ann.to_csv('ann.csv', index=False)