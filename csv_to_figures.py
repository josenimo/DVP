#!C:/Users/jnimoca/Anaconda3/envs/scimap/bin/python

import sys
import os
import anndata as ad
import pandas as pd
import scanpy as sc
import seaborn as sns; sns.set(color_codes=True)
import scimap as sm

#os.chdir ("/home/josenimo/DVP")

adata = sm.pp.mcmicro_to_scimap ("TMA_01_data/unmicst-TMA_01_Core_01_cell.csv")

sc.pl.highest_expr_genes(adata, n_top=20)

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10)

sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 1)
sc.pl.umap(adata, color='leiden')

sc.set_figure_params(dpi=500, color_map = 'viridis_r', figsize=(6,6))
sc.pl.scatter(adata, x='X_centroid', y='Y_centroid', color='leiden')

