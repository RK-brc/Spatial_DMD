#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import anndata
import os
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import liana as li


# In[ ]:


adata_vis=sc.read_h5ad("~/adata_final.h5ad")


# In[ ]:


adata_vis.var_names = adata_vis.var['gene_symbols'].astype(str)
adata_vis.var_names_make_unique()


# ## LIANA rank_aggregate per sample, grouping by niche

# In[ ]:


# import liana's rank_aggregate
from liana.mt import rank_aggregate
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean


# In[ ]:


#setting the parameters 
N_SAMPLES = adata_vis.obs['sample_ID'].nunique() 
SAMPLE_THRESHOLD = 4
SLICE_THRESHOLD = 1
RANK_THRESHOLD = 0.05
BANDWIDTH = 263      
CUTOFF = 0.1         
NZ_PROP = 0.2    
K_CV = 10      

slice_to_sample = (
    adata_vis.obs[['sample', 'sample_ID']]
    .drop_duplicates()
    .set_index('sample')['sample_ID']
)


# In[ ]:


#rank aggregate per sample
results = []

for sample in adata_vis.obs['sample_ID'].unique():
    adata_sample = adata_vis[adata_vis.obs['sample_ID'] == sample].copy()

    li.mt.rank_aggregate(
        adata_sample,
        groupby='Niche_NMF',
        resource_name='consensus',
        expr_prop=0.01,
        min_cells=1,
        use_raw=False,
        verbose=False
    )

    res = adata_sample.uns['liana_res'].copy()
    res['sample'] = sample
    results.append(res)
    print(f"Done: {sample}")

liana_res = pd.concat(results, ignore_index=True)


# In[ ]:


liana_res.to_csv("~/liana_rank_aggregate_all_samples.csv", index=False)


# In[ ]:


#count how many samples support each LR pair per niche 
sample_counts = (
    liana_res[liana_res['is_significant']]
    .groupby(['source', 'target', 'ligand_complex', 'receptor_complex'])
    ['sample'].nunique()
    .reset_index()
    .rename(columns={'sample': 'n_samples'})
)

robust_lr_pairs = sample_counts[sample_counts['n_samples'] >= 4].copy()
robust_lr_pairs['sample_consistency'] = robust_lr_pairs['n_samples'] / N_SAMPLES

print(f"Robust LR pairs (>={4}/6 samples): {len(robust_lr_pairs)}")


# In[ ]:


from liana.method import MistyData, genericMistyData, lrMistyData
from liana.method.sp import LinearModel, RandomForestModel


# In[ ]:


#preparing intra view for prediction using MISTy
cell_type_names = list(adata_vis.obsm['means_cell_abundance_w_sf'].columns) \
    if hasattr(adata_vis.obsm['means_cell_abundance_w_sf'], 'columns') \
    else [f'celltype_{i}' for i in range(
        adata_vis.obsm['means_cell_abundance_w_sf'].shape[1])]

print(f"Cell types: {cell_type_names}")

#convert to AnnData using LIANA utility
comps = li.ut.obsm_to_adata(adata_vis, 'means_cell_abundance_w_sf')
print(f"Cell type AnnData shape: {comps.shape}")
print(f"Cell types as vars: {comps.var_names.tolist()}")


# In[ ]:


#prepare ligand expression data
#these are the biologically relevant ligands in the tissue

#get unique ligands from robust pairs
robust_ligands = robust_lr_pairs['ligand_complex'].unique().tolist()

#filter to ligands that exist in adata_vis
robust_ligands_found = [l for l in robust_ligands if l in adata_vis.var_names]
print(f"Robust ligands from rank aggregate: {len(robust_ligands)}")
print(f"Found in adata_vis: {len(robust_ligands_found)}")

#create ligand expression AnnData
ligand_adata = adata_vis[:, robust_ligands_found].copy()
print(f"Ligand AnnData shape: {ligand_adata.shape}")


# In[ ]:

# ANALYSIS B: lrMistyData- Receptors in intra, ligands in para

lr_target_metrics_list = []
lr_interactions_list = []

for sample in adata_vis.obs['sample'].unique():
    print(f"\nProcessing lrMistyData: {sample}")

    slice_idx = adata_vis.obs['sample'] == sample
    adata_slice = adata_vis[slice_idx].copy()

    try:
        misty_lr = lrMistyData(
            adata_slice,
            bandwidth=BANDWIDTH,
            cutoff=CUTOFF,
            set_diag=False,
            nz_threshold=0.1,
            resource_name='consensus'
        )

        misty_lr(
            bypass_intra=True,
            model=LinearModel,
            k_cv=K_CV,
            seed=345,
            verbose=False
        )

        tm = misty_lr.uns['target_metrics'].copy()
        tm['sample'] = sample
        tm['donor'] = slice_to_sample[sample]
        lr_target_metrics_list.append(tm)

        interactions = misty_lr.uns['interactions'].copy()
        interactions['sample'] = sample
        interactions['donor'] = slice_to_sample[sample]
        lr_interactions_list.append(interactions)

        print(f"  Done: {sample}")

    except Exception as e:
        print(f"  Failed: {sample} — {e}")
        continue

lr_target_metrics_df = pd.concat(lr_target_metrics_list, ignore_index=True)
lr_interactions_df = pd.concat(lr_interactions_list, ignore_index=True)

print(f"lr_target_metrics_df shape: {lr_target_metrics_df.shape}")
print(f"lr_interactions_df shape: {lr_interactions_df.shape}")
print(lr_interactions_df.head())

#Aggregate lrMistyData results
lr_donor_interactions = (
    lr_interactions_df
    .groupby(['target', 'predictor', 'view', 'donor'])
    ['importances']
    .median()
    .reset_index()
)

lr_final_interactions = (
    lr_donor_interactions
    .groupby(['target', 'predictor', 'view'])
    .agg(
        mean_importance=('importances', 'mean'),
        std_importance=('importances', lambda x: x.std() if len(x) > 1 else 0),
        n_donors=('donor', 'nunique'),
        n_donors_positive=('importances', lambda x: (x > 0).sum()),
        n_donors_negative=('importances', lambda x: (x < 0).sum())
    )
    .reset_index()
)


# In[ ]:


#Filter to extra view 
lr_final_extra = lr_final_interactions[
    lr_final_interactions['view'] == 'extra'
].sort_values('mean_importance', key=abs, ascending=False)

#Cross-validate with rank aggregate
rank_agg_pairs = set(zip(
    robust_lr_pairs['ligand_complex'],
    robust_lr_pairs['receptor_complex']
))
misty_pairs = set(zip(
    lr_final_extra[lr_final_extra['n_donors'] >= SAMPLE_THRESHOLD]['predictor'],
    lr_final_extra[lr_final_extra['n_donors'] >= SAMPLE_THRESHOLD]['target']
))

overlap = rank_agg_pairs.intersection(misty_pairs)
print(f"\nLR pairs in both rank aggregate AND lrMistyData: {len(overlap)}")
print("These are highest confidence interactions:")
for pair in sorted(overlap):
    print(f"  {pair[0]} → {pair[1]}")

lr_robust = lr_final_extra[
    lr_final_extra['n_donors'] >= SAMPLE_THRESHOLD
].copy()


# In[ ]:


#ANALYSIS A: genericMistyData- Intra view=cell type abundances, Para view=ligand expression in neighbourhood


target_metrics_list = []
interactions_list_misty = []

for sample in adata_vis.obs['sample'].unique():
    print(f"\nProcessing: {sample}")

    #Subset to this slice
    slice_idx = adata_vis.obs['sample'] == sample

    #Cell type abundances for this slice
    comps_slice = comps[slice_idx].copy()
    comps_slice.obsm['spatial'] = adata_vis[slice_idx].obsm['spatial']

    #Ligand expression for this slice
    ligand_slice = ligand_adata[slice_idx].copy()
    ligand_slice.obsm['spatial'] = adata_vis[slice_idx].obsm['spatial']

    #Compute spatial neighbours on ligand slice
    li.ut.spatial_neighbors(
        ligand_slice,
        bandwidth=BANDWIDTH,
        cutoff=CUTOFF,
        kernel='gaussian',
        set_diag=False,    
        spatial_key='spatial'
    )

    #Check connectivity
    conn = ligand_slice.obsp['spatial_connectivities']
    neighbours_per_spot = np.diff(conn.indptr)
    print(f"  Mean effective neighbours: {neighbours_per_spot.mean():.1f}")

    #Build MISTy object
    try:
        misty = genericMistyData(
            intra=comps_slice,
            extra=ligand_slice,
            bandwidth=BANDWIDTH,
            cutoff=CUTOFF,
            set_diag=False
        )

        #Run MISTy with LinearModel
        #signed t-values give directionality
        misty(
            model=LinearModel,
            k_cv=K_CV,
            seed=345,
            bypass_intra=True,
            verbose=False
        )

        #Extract results
        tm = misty.uns['target_metrics'].copy()
        tm['sample'] = sample
        tm['donor'] = slice_to_sample[sample]
        target_metrics_list.append(tm)

        interactions = misty.uns['interactions'].copy()
        interactions['sample'] = sample
        interactions['donor'] = slice_to_sample[sample]
        interactions_list_misty.append(interactions)

        print(f"  Done: {sample}")

    except Exception as e:
        print(f"  Failed: {sample} — {e}")
        continue

#Concatenate results across slices
target_metrics_df = pd.concat(target_metrics_list, ignore_index=True)
interactions_df = pd.concat(interactions_list_misty, ignore_index=True)

print(f"\nTarget metrics shape: {target_metrics_df.shape}")
print(f"Interactions shape: {interactions_df.shape}")
print(f"\nTarget metrics columns: {target_metrics_df.columns.tolist()}")
print(f"Interactions columns: {interactions_df.columns.tolist()}")


# In[ ]:


# AGGREGATE RESULTS ACROSS SLICES 

#Summarise gain_R2 per cell type per sample (median across slices)
sample_target_metrics = (
    target_metrics_df
    .groupby(['target', 'sample'])
    .agg(
        gain_R2=('gain_R2', 'median'),
        intra_R2=('intra_R2', 'median'),
        multi_R2=('multi_R2', 'median')
    )
    .reset_index()
)

#final summary across samples
final_target_metrics = (
    sample_target_metrics
    .groupby('target')
    .agg(
        mean_gain_R2=('gain_R2', 'mean'),
        std_gain_R2=('gain_R2', 'std'),
        median_gain_R2=('gain_R2', 'median'),
        n_donors=('donor', 'nunique')
    )
    .reset_index()
    .sort_values('mean_gain_R2', ascending=False)
)

print("\nFinal target metrics (sorted by mean gain_R2):")


# In[ ]:


## Interactions aggregation
#For each target-predictor-view combination:
#median importance across slices within sample
#then mean across sample
#count how many samples show positive importance

sample_interactions = (
    interactions_df
    .groupby(['target', 'predictor', 'view', 'donor'])
    ['importances']
    .median()
    .reset_index()
)

final_interactions = (
    sample_interactions
    .groupby(['target', 'predictor', 'view'])
    .agg(
        mean_importance=('importances', 'mean'),
        std_importance=('importances', 'std'),
        median_importance=('importances', 'median'),
        n_donors=('donor', 'nunique'),
        n_donors_positive=('importances', lambda x: (x > 0).sum()),
        n_donors_negative=('importances', lambda x: (x < 0).sum())
    )
    .reset_index()
)

#Filter to para view- directionality signal
final_para = final_interactions[
    final_interactions['view'] == 'para'
].sort_values('mean_importance', ascending=False)

print(f"\nTop para interactions (ligand→cell type):")
print(final_para.head(20))


# In[ ]:


#Robust interactions- consistent across samples
robust_para = final_para[
    final_para['n_donors'] >= SAMPLE_THRESHOLD
].copy()
print(f"\nRobust para interactions (>={SAMPLE_THRESHOLD} samples): {len(robust_para)}")


# In[ ]:


#output dir
os.makedirs('~/misty_results', exist_ok=True)

# ================================================================
# 1. CORE DATAFRAMES (Analysis A)
# ================================================================
final_target_metrics.to_csv('~/misty_results/A_final_target_metrics.csv') #Ranks cell types by how much they are driven by neighbourhood ligands across samples
robust_para.to_csv('~/misty_results/A_robust_para_interactions.csv') #ligand–cell type relationships that are strong, reproducible across samples
target_metrics_df.to_csv('~/misty_results/A_target_metrics_all.csv')  #which cell types are most shaped by their neighbourhood (high gain_R2)
sample_target_metrics.to_csv('~/misty_results/A_sample_target_metrics.csv') #sample level summary. how consistently do the ligands in the spot’s neighbourhood explain the abundance of cell types?


# ================================================================
# 2. CORE DATAFRAMES (Analysis B)
# ================================================================
lr_final_extra.to_csv('~/misty_results/B_lr_final_extra.csv') #ligand→receptor spatial CCC table.
lr_robust.to_csv('~/misty_results/B_lr_robust_interactions.csv') #same as above, high‑confidence ligand→receptor spatial set
lr_target_metrics_df.to_csv('~/misty_results/B_lr_target_metrics_all.csv') #which receptors are most strongly shaped by spatial signalling, slice‑by‑slice.

# ================================================================
# 3. INPUTS & METADATA 
# ================================================================
comps.write('~/misty_results/celltype_abundances.h5ad')  #Intra view
ligand_adata.to_df().T.to_csv('~/misty_results/ligand_expression.csv')  #Ligands
pd.Series(robust_ligands_found).to_csv('~/misty_results/robust_ligands.csv')
slice_to_sample.to_csv('~/misty_results/slice_to_donor.csv')

# ================================================================
# 4. SUMMARY OVERLAPS & PIVOTS
# ================================================================
pd.DataFrame(overlap, columns=['ligand', 'receptor']).to_csv('~/misty_results/rankagg_misty_overlap.csv') #highest‑confidence interactions


