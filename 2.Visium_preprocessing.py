#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')


# In[ ]:


sample_data = pd.read_csv('/sample_annot.csv')


# In[ ]:


#set the path of the folders
sp_data_folder = '/Visium_dmd/'


# In[ ]:


def read_and_qc(sample_name, path=sp_data_folder):
    """
    reads the data into the anndata object and calculates QC metrics. 
    """

    adata = sc.read_visium(path + str(sample_name) + '/outs/',
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    adata.var.drop(columns='ENSEMBL', inplace=True)

    # just in case there are non-unique ENSEMBL IDs
    adata.var_names_make_unique()

    #Calculate QC metrics
    from scipy.sparse import csr_matrix
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    adata.var['MT'] = [gene.startswith('MT-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['MT'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    #add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata


# In[ ]:


def select_slide(adata, s, s_col='sample'):
    """ 
    selects the data for one slide from the spatial anndata object.
    """

    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]

    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}

    return slide


# In[ ]:


# Read the data into anndata objects
slides = []
for i in sample_data['Sample']:
    slides.append(read_and_qc(i, path=sp_data_folder))


# In[ ]:


# Combine anndata objects together
adata_vis = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_data['Sample'],
    index_unique=None
)


# In[ ]:


# mitochondria-encoded (MT) genes should be removed for spatial mapping
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT'].values]


# In[ ]:


sc.pp.calculate_qc_metrics(adata_vis, inplace=True)


# In[ ]:


# PLOT QC FOR EACH SAMPLE for DMD
fig, axs = plt.subplots(len(slides), 4, figsize=(15, 4*len(slides)-4))
for i, s in enumerate(adata_vis.obs['sample'].unique()):

    slide = select_slide(adata_vis, s)
    sns.distplot(slide.obs['total_counts'],
                 kde=False, ax = axs[i, 0])
    axs[i, 0].set_xlim(0, adata_vis.obs['total_counts'].max())
    axs[i, 0].set_xlabel(f'total_counts | {s}')

    sns.distplot(slide.obs['total_counts']\
                 [slide.obs['total_counts']<20000],
                 kde=False, bins=40, ax = axs[i, 1])
    axs[i, 1].set_xlim(0, 20000)
    axs[i, 1].set_xlabel(f'total_counts | {s}')

    sns.distplot(slide.obs['n_genes_by_counts'],
                 kde=False, bins=60, ax = axs[i, 2])
    axs[i, 2].set_xlim(0, adata_vis.obs['n_genes_by_counts'].max())
    axs[i, 2].set_xlabel(f'n_genes_by_counts | {s}')

    sns.distplot(slide.obs['n_genes_by_counts']\
                 [slide.obs['n_genes_by_counts']<6000],
                 kde=False, bins=60, ax = axs[i, 3])
    axs[i, 3].set_xlim(0, 6000)
    axs[i, 3].set_xlabel(f'n_genes_by_counts | {s}')

plt.tight_layout()


# In[ ]:


keep = (adata_vis.obs['total_counts'] >= 750) & (adata_vis.obs['n_genes_by_counts'] >= 200)
adata_vis = adata_vis[keep].copy()


# In[ ]:


#Add counts matrix as adata.raw
adata_vis.raw = adata_vis


# In[ ]:


#normalize and logtransform the data
sc.pp.normalize_total(adata_vis)
sc.pp.log1p(adata_vis)

