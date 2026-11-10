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

os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32,openmp=True,force_device=True'

import cell2location

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq


# In[ ]:


adata_vis = sc.read_h5ad("~/DMD_withFAPs/output/dmdmodel/sp.h5ad")  


# Run NMF

# In[ ]:


from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(5,8), 
      'sample_name_col': 'sample', 
      'n_restarts': 3 
    },
    model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
    export_args={'path': f'~/DMD_withFAPs/output/NMF'}
)


# In[ ]:


res_dict['n_fact7']['mod'].plot_cell_type_loadings() # check this for all the factor loadings
plt.show()


# In[ ]:


#make the data more readable

#loading the factor loading for factor7 file
df = pd.read_csv('~/DMD_withFAPs/output/NMF/location_factors_mean/n_fact7.csv', header=0)

#rename the factor columns to Niche1, Niche2 etc
df.columns = ['spot_id'] + [f'Niche{i}' for i in range(len(df.columns)-1)]

#find the Niche (column) with the highest value for each spot
df['Niche_NMF'] = df[[col for col in df.columns if col.startswith('Niche')]].idxmax(axis=1)

#extract the factor number from the Niche string
df['Niche_NMF'] = df['Niche_NMF'].str.extract(r'Niche(\d+)').astype(int)

#ensure spot_id is the index in the NMF factor DataFrame
df.set_index('spot_id', inplace=True)

#assign all niche columns (excluding 'Niche_NMF') to adata_vis.obs
niche_cols = df.columns.drop('Niche_NMF')  # or use list of niche col names

adata_vis.obs[niche_cols] = df[niche_cols]
adata_vis.obs['Niche_NMF'] = adata_vis.obs.index.map(df['Niche_NMF'])
adata_vis.obs['Niche_NMF'] = adata_vis.obs['Niche_NMF'].astype(int)
adata_vis.obs['Niche_NMF'] = adata_vis.obs['Niche_NMF'].replace({i: f'Niche{i}' for i in range(adata_vis.obs['Niche_NMF'].max() + 1)})


# To get the dotplot

# In[ ]:


factor_loadings_df = pd.read_csv('~/DMD_withFAPs/output/NMF/cell_type_fractions_mean/n_fact7.csv', skiprows=1, names=column_names)


# In[ ]:


factor_loadings_df_percentage = factor_loadings_df.set_index('cell_type') * 100
df_melted = factor_loadings_df_percentage.reset_index().melt(
    id_vars='cell_type',
    var_name='Niche',
    value_name='Percentage'
)

plt.figure(figsize=(6, 6))
scatter = plt.scatter(
    x=df_melted["Niche"],
    y=df_melted["cell_type"],
    s=df_melted["Percentage"],           
    c=df_melted["Percentage"],         
    cmap="RdPu",
    alpha=1
)

cbar = plt.colorbar(scatter)

#format plot
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.xlabel("")
plt.ylabel("")
plt.tight_layout()
plt.savefig("niche_celltype_dotplot.tiff", dpi=300, bbox_inches="tight", pad_inches=0.05)
plt.show()


# To get the % of cell type in each niche

# In[ ]:


#normalise each niche
df_normalized = factor_loadings_df.div(factor_loadings_df.sum(axis=1), axis=0) * 100 

plt.figure(figsize=(8, 4))

sns.heatmap(df_normalized, 
            cmap="RdPu",  
            annot=True,  
            fmt=".1f",  
            linewidths=0.5,  
            linecolor='gray', 
            cbar_kws={'label': 'Percentage'}, 
            xticklabels=True,
            yticklabels=True,  
            annot_kws={"size": 8})  

plt.yticks(rotation=0)  
plt.tight_layout()

plt.show()


# In[ ]:


#get unique niche labels
niches = adata_vis.obs['Niche_NMF'].unique()
palette = sns.color_palette(["#E69F00", "#56B4E9", "#009E73",
                             "#F0E442", "#0072B2", "#D55E00", "#CC79A7"])
niche_colors = dict(zip(niches, palette))
niche_colors = dict(sorted(niche_colors.items()))


# In[ ]:


adata_vis.obs['Niche_NMF'] = adata_vis.obs['Niche_NMF'].astype('category')
adata_vis.uns['Niche_NMF_colors'] = [niche_colors[n] for n in adata_vis.obs['Niche_NMF'].cat.categories]


# Visualise the niches

# In[ ]:


from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'slice_name')
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(
        slide,
        color='Niche_NMF',
        cmap='tab20',
        size=1.3,
        img_key='hires', 
        vmin=0, vmax='p99.2'
)


# In[ ]:


niches = ['Niche0', 'Niche1', 'Niche2', 'Niche3', 'Niche4', 'Niche5', 'Niche6']

with mpl.rc_context({'figure.figsize': [2.3, 2.5],   # per-panel size
                     'font.size': 6}):
    fig = sc.pl.spatial(
        slide,
        color=niches,
        cmap='magma_r',
        size=1.3,
        ncols=3,
        vmin=0,
        vmax='p99.2',
        show=False,
        return_fig=True,
    )

    print(fig.get_size_inches())  # check size before saving

    fig.savefig('~/Plots_paper/all_niches.png', dpi=300,
                bbox_inches='tight', pad_inches=0.05, facecolor='white') 


# ## Neighbourhood enrichment based on the NMF factors we have generated

# In[ ]:


adata_vis.obs['Niche_NMF'] = adata_vis.obs['Niche_NMF'].astype('category')


# one spot radius

# In[ ]:


sq.gr.spatial_neighbors(adata_vis, coord_type="grid", spatial_key="spatial",library_key="sample",n_rings=1)
sq.gr.nhood_enrichment(adata_vis, cluster_key="Niche_NMF")


# In[ ]:


onespot_matrix = pd.DataFrame(adata_vis.uns['Niche_NMF_nhood_enrichment']['zscore'])
plt.figure(figsize=(8, 4))
sb.heatmap(onespot_matrix, annot=True, cmap="RdBu", fmt=".2f" )
plt.yticks(rotation=0) 
plt.show()


# three spot radius

# In[ ]:


sq.gr.spatial_neighbors(adata_vis, coord_type="grid", spatial_key="spatial",library_key="sample",n_rings=3)
sq.gr.nhood_enrichment(adata_vis, cluster_key="Niche_NMF")


# In[ ]:


threespot_matrix= pd.DataFrame(adata_vis.uns['Niche_NMF_nhood_enrichment']['zscore'])
plt.figure(figsize=(8, 4))
sb.heatmap(threespot_matrix, annot=True, cmap="RdBu", fmt=".2f" )
plt.yticks(rotation=0) 
plt.show()


# five spot radius

# In[ ]:


sq.gr.spatial_neighbors(adata_vis, coord_type="grid", spatial_key="spatial",library_key="sample",n_rings=5)
sq.gr.nhood_enrichment(adata_vis, cluster_key="Niche_NMF")


# In[ ]:


fivespot_matrix = pd.DataFrame(adata_vis.uns['Niche_NMF_nhood_enrichment']['zscore'])
plt.figure(figsize=(8, 4))
sb.heatmap(fivespot_matrix, annot=True, cmap="RdBu", fmt=".2f" )
plt.yticks(rotation=0) 
plt.show()


# # Pathway enrichment

# In[ ]:


adata_vis.var_names = adata_vis.var['gene_symbols'].astype(str)
adata_vis.var_names_make_unique()


# In[ ]:


sc.tl.rank_genes_groups(adata_vis, groupby = 'Niche_NMF', groups = "all", use_raw = False, method = "wilcoxon", pts= True,)
sc.pl.rank_genes_groups(adata_vis)


# In[ ]:


#combine into one table
result = adata_vis.uns['rank_genes_groups']
markers = []
for cluster in result['names'].dtype.names:
    current = pd.DataFrame({"gene": result["names"][cluster], "score": result["scores"][cluster],
                            "logfoldchange": result["logfoldchanges"][cluster], "pval": result["pvals"][cluster],
                            "pval_adj": result["pvals_adj"][cluster],
                            "pct_within":result["pts"].loc[result["names"][cluster]][cluster],
                            "pct_outside":result["pts_rest"].loc[result["names"][cluster]][cluster],  "cluster": cluster})
    markers.append(current)
markers = pd.concat(markers)


# In[ ]:


# Filter for logfoldchange > 0.5 
top_markers = markers[(markers['logfoldchange'] > 0.5) & (markers['pval_adj'] < 0.05)]

# Sort by cluster and logfoldchange (descending)
top_markers.sort_values(by=["cluster", "logfoldchange"], ascending=[True, False], inplace=True)

#Save anndata object with results
adata_vis.write("adata_final.h5ad")

