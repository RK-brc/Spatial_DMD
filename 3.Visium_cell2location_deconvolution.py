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

os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'

import cell2location
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')


# In[ ]:


#load the visium object
adata_vis = sc.read("~/Visium_objects/visium_dmd_final.h5ad")  


# In[]:


#read the data 
adata_ref = sc.read_h5ad('~/Spatial/dmd_withfaps.h5ad')


# In[ ]:


#Check raw counts (should be integers)
if adata_ref.raw is not None:
    print("Raw counts exist in adata.raw.X")
    print(f"adata.raw.X dtype: {adata_ref.raw.X.dtype}") 

#Check normalized values
print(f"adata.X dtype: {adata_ref.X.dtype}") 


# In[13]:


print(adata_ref.X[:5,:5])


# In[ ]:


dense_X = adata_ref.X.toarray() if hasattr(adata_ref.X, "toarray") else adata_ref.X
print(dense_X[:10, :10])


# In[ ]:


#make the data compatible for cell2location
adata_ref.X = adata_ref.raw.X.copy()
adata_ref.X = adata_ref.X.astype(np.int64)


# In[11]:


print("adata_ref.X dtype:", adata_ref.X.dtype)


# In[ ]:


#prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        batch_key='orig.ident',
                        labels_key='merge_type',
                       )


# In[ ]:


#create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

mod.view_anndata_setup()


# In[]:

mod.train(max_epochs=250, accelerator="gpu")

# In[]:

#check the ELBO loss 
mod.plot_history(20)

# In[]:
#export the estimated cell abundance
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1500, 'batch_size': 2500, 'accelerator': 'gpu'}
)


# In[]:
mod.plot_QC()


# In[ ]:

#export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']


# In[ ]:
dir= '~/DMD'


# In[]:
mod.save(dir+"/rsignatures", overwrite=True)

#Save anndata object with results
adata_ref.write(dir+"/rsignatures/sc.h5ad")


# Spatial mapping

# In[ ]:

#from reference signature, remove hgnc_sybmol and ensemble id columns
inf_aver = inf_aver.drop(columns=['hgnc_symbol','ensembl_gene_id' ])

# In[ ]:


#find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# In[ ]:


#prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")


# In[ ]:


#create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    N_cells_per_location=10,
    detection_alpha=20
)
mod.view_anndata_setup()


# In[ ]:


mod.train(max_epochs=30000,
          batch_size=None,
          train_size=1,
          accelerator="cpu",
         )

output_dir = '~/DMD_withFAPs/output'

#Save model
mod.save(output_dir+"/dmdmodel", overwrite=True)


# In[ ]:

#export the estimated cell abundance 
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'accelerator': 'cpu'}
)


# In[ ]:

mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_vis.write(output_dir+"/dmdmodel/sp.h5ad")


# In[ ]:

adata_vis.obsm['q05_cell_abundance_w_sf'].to_csv('~/DMD_withFAPs/output/q05_cell_abundance_w_sf.csv')


# In[ ]:

adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

#select slide
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'sample1_name')

# plot
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(slide, cmap='magma',
                  color=['Slow.Fibres', 'Endothelial.Cells', 'Fast.Fibres', 'Smooth.Muscle.Cells',
                         'Inflammatory.FAPs','MME..FAPs','Pre.Adipogenic','Mixed.Fibres' ,'Satellite.Cells', 
                         'Inflammatory.Cells','Adipocytes','Regenerative.Fibres'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2'
                 )


# Spot-level normalization

# In[]:

cell_type_cols = ['Slow.Fibres', 'Endothelial.Cells', 'Fast.Fibres', 'Smooth.Muscle.Cells','Mixed.Fibres' ,
'Satellite.Cells', 'Inflammatory.Cells','Adipocytes','Regenerative.Fibres', 'Inflammatory.FAPs','MME..FAPs','Pre.Adipogenic']


# In[]:

#Compute the total abundance for each spot
row_sum = adata_vis.obs[cell_type_cols].sum(axis=1)

#For each cell type, compute the percentage relative to the total abundance in that spot.
for cell_type in cell_type_cols:
    adata_vis.obs[cell_type + "_perspot"] = 100 * adata_vis.obs[cell_type] / row_sum.replace(0, np.nan)


