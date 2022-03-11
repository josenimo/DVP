#!C:/Users/jnimoca/Anaconda3/envs/scimap/bin/python

import sys
import os
import anndata as ad
import pandas as pd
import scanpy as sc
import seaborn as sns; sns.set(color_codes=True)
import scimap as sm
print("Imports completed")

#os.chdir ("/home/josenimo/DVP")
print("Transforming MCMICRO csv to adata object..")
adata = sm.pp.mcmicro_to_scimap ("TMA_01_data/unmicst-TMA_01_Core_01_cell.csv")
print("Ranked expression data")
sc.pl.highest_expr_genes(adata, n_top=20, save=True)

print("Creating neighborhood graph")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10)
print("Creating UMAP")
sc.tl.umap(adata)
print("Creating Leiden Categorization, res=0.3")
sc.tl.leiden(adata, resolution = 0.3)

sc.set_figure_params(scanpy=True, dpi=200, dpi_save=500, frameon=True, 
                                vector_friendly=True, fontsize=14, figsize=(6,6), 
                                color_map='viridis_r', format='pdf', facecolor=None, 
                                transparent=False, ipython_format='png2x')

print("Plotting UMAP, coloured with Leiden categorization")
sc.pl.umap(adata, color='leiden', save=True)
print("Plotting centroid data in scatterplot, coloured with leiden categorization ")
sc.pl.scatter(adata, x='X_centroid', y='Y_centroid', color='leiden', save=True)
print("done")
