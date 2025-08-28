#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import scanpy as sc
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import cell2location
from cell2location.models import RegressionModel
from matplotlib import rcParams
from loguru import logger
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

#take in arguments
import argparse
parser = argparse.ArgumentParser(description='Cell2location deconvolution script')
parser.add_argument('--dataset', type=int, default=1, help='Dataset to deconvolve')

args = parser.parse_args()
n_dataset = args.dataset
logger.info(f"Running cell2location deconvolution for dataset {n_dataset}")


n_cells = 10
n_epochs = 5000
# create paths and names to results folders for reference regression and cell2location models
##cell type
#ref_run_name = f'/data/NCI_LP/04.single_cell/01.RNA/MasahiNomura_NatureGenetics_2025/cell2location/reference_signatures/'
#results_folder = './ref_Masahi_NG_2025/'
###malignant cell states
ref_run_name = f'/data/NCI_LP/04.single_cell/01.RNA/MasahiNomura_NatureGenetics_2025/cell2location_malignantcellstate/reference_signatures/'
results_folder = './ref_Masahi_NG_2025_cellstates/'
###
run_name = f"{results_folder}/cell2location_map_n{n_cells}cell/"

if not os.path.exists(run_name):
    os.makedirs(run_name)

##Load model
if(False): 
    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref = sc.read_h5ad(adata_file)
    mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
    cell_types = adata_ref.uns['mod']['factor_names']
    del adata_ref #required


sp_data_list_file = "/data/wuz6/project/25.pedHGG_spatial/01.preproess/01.cluster2025HE/00.data.folder.csv"
sample_data = pd.read_csv(sp_data_list_file)
#fset first column name as folder
sample_data.rename(columns={sample_data.columns[0]: 'folder'}, inplace=True)

sample_data.head()


def read_and_qc(sample_name, path):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """

    adata = sc.read_visium(path + '/outs/',
    #adata = sq.read.visium(path + '/outs/',
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    #return adata
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    #adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    #adata.var_names = adata.var['ENSEMBL']
    #adata.var.drop(columns='ENSEMBL', inplace=True)
    adata.var_names = adata.var['SYMBOL']

    # Calculate QC metrics
    #from scipy.sparse import csr_matrix
    #adata.X = adata.X.toarray()
    #sc.pp.calculate_qc_metrics(adata, inplace=True)
    #adata.X = csr_matrix(adata.X)
    #adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    #adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    sc.pp.filter_cells(adata, min_counts = 1000)
    sc.pp.filter_genes(adata, min_cells=5)
    print(adata.uns['spatial'].keys())
    return adata

def select_slide(adata, s, s_uns,s_col='sample'):
    r""" This function selects the data for one slide from the spatial anndata object.

    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param s_col: column in adata.obs listing experiment name for each location
    """

    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[s_uns in k for k in s_keys]][0]
    #uns['spatial'].keys some case has defaulted by different name, need to manualy select by read their names
    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]} 
    slide.obsm['q05_cell_abundance_w_sf'] = adata.obsm['q05_cell_abundance_w_sf'][adata.obs[s_col].isin([s])]
    #adata_vis.obsm

    return slide



#[s_name].isin(adata_vis.obsm['q05_cell_abundance_w_sf'].index)
#s_name
#adata_vis.obsm['q05_cell_abundance_w_sf'].iloc[0:10,]
#adata_vis
#i = 0
#s_name = sample_name_list[keep_idx[i]]
#adata_vis1 = select_slide(adata_vis, s_name, s_uns= uns_keys[i])
#adata_vis1.obs[cell_types] = adata_vis1.obsm['q05_cell_abundance_w_sf']
#sc.pl.spatial(adata_vis1, color= list(cell_types) )


#######################
sample_name_list = sample_data.Sample.values
folder_list = sample_data.folder.values
sample_source = sample_data.Source.values
sample_source


if n_dataset == 3:
    # select samples from Alissa et al. 2023
    keep_idx = [i for i in range(len(sample_source)) if sample_source[i] == "Alissa et al." ]
    out_prefix = 'data_Alissa'
elif n_dataset == 2:
    keep_idx = [i for i in range(len(sample_source)) if sample_source[i] == "Ravi et al." ]
    out_prefix = 'data_Ravi'
elif n_dataset == 1:
    keep_idx = [i for i in range(len(sample_source)) if sample_source[i] == "Aldape Lab" ]
    out_prefix = 'data_Aldape'
else:
    raise ValueError("Invalid dataset number. Please choose 1, 2, or 3.")

#i = keep_idx[0]
#path = folder_list[i]
#sample_name = sample_name_list[i]
#adata = sc.read_visium(path + '/outs/',
#                       count_file='filtered_feature_bc_matrix.h5', load_images=True)


#slides = []
#for i in keep_idx:
#    folder = folder_list[i]
#    sample_name = sample_name_list[i]
#    print("Reading ", sample_name, ":", folder)
#    slides.append(read_and_qc(sample_name, path = folder))


# Combine anndata objects together
#sample_names = sample_data['Sample'][keep_idx]
#adata_vis = slides[0].concatenate(
#    slides[1:],
#    batch_key="sample",
#    uns_merge="unique",
#    batch_categories = sample_names,
#    index_unique=None
#)
#######################

# find shared genes and subset both anndata and reference signatures
inf_aver = pd.read_csv(f"{ref_run_name}/inf_aver.csv", index_col = 0)


for i in keep_idx:
    folder = folder_list[i]
    sample_name = sample_name_list[i]
    #out_prefix1 = out_prefix + "_" +sample_name
    out_prefix1 = sample_name
    #check if file exists        
    if os.path.exists(f"./{run_name}/{out_prefix1}_decon_mod.train_history.png"):
        continue
    adata_vis = read_and_qc(sample_name, path = folder)
    # find mitochondria-encoded (MT) genes
    adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]
    #adata_vis.obs['sample']
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    print(f"Intersect genes: {len(intersect)}: {intersect[:10]}")
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver1 = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver1, 
        # the expected average cell abundance: tissue-dependent 
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location= n_cells,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
        ) 
    mod.view_anndata_setup()

    mod.train(max_epochs= n_epochs,
          # train using full data (batch_size=None)
          batch_size=None,  #None: use all spots, 
          # use all data points in training because 
          # we need to estimate cell abundance at all locations
          train_size=1          
         )

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history()
    plt.legend(labels=['full data training']);
    plt.savefig(f"./{run_name}/{out_prefix1}_decon_mod.train_history.png",bbox_inches='tight')
    plt.close()
    ### export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
    )

    # Save model
    mod.save(f"{run_name}/{sample_name}", overwrite=True)

    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
    # "Alissa et al." , "Ravi et al." "Aldape Lab" 
    #keep_idx1
    # Save anndata object with results
    adata_file = f"{run_name}/{out_prefix1}_sp.h5ad"
    adata_vis.write(adata_file)
    adata_file

    plt.savefig(f"{run_name}/{out_prefix1}.reconstruction_accuracy.png",
                    bbox_inches='tight')
    plt.close()

    #cell_types = adata_vis.uns['mod']['factor_names']
    #adata_vis.obs[cell_types] = adata_vis.obsm['q05_cell_abundance_w_sf']
    #adata_vis.obsm['means_cell_abundance_w_sf']
    cell_types = adata_vis.obsm['q05_cell_abundance_w_sf'].columns.to_list()
    cell_types = [x.replace("q05cell_abundance_w_sf_", "") for x in cell_types]
    #cell_types[0].replace("q05cell_abundance_w_sf_", "") #, cell_types)
    #cell_types
    #adata_vis.obsm['means_cell_abundance_w_sf']
    #5% quantile of the posterior distribution, representing the value of cell abundance that 
    #the model has high confidence in (aka 'at least this amount is present').
    cell_decon_q05 = adata_vis.obsm['q05_cell_abundance_w_sf']
    cell_decon_q05.columns = cell_types
    cell_decon_q05.to_csv(f"{run_name}/{out_prefix1}_decon_q05.csv")
    ###
    cell_decon_q95 = adata_vis.obsm['q95_cell_abundance_w_sf']
    cell_decon_q95.columns = cell_types
    cell_decon_q95.to_csv(f"{run_name}/{out_prefix1}_decon_q95.csv")
    ###
    cell_decon_means = adata_vis.obsm['means_cell_abundance_w_sf']
    cell_decon_means.columns = cell_types
    cell_decon_means.to_csv(f"{run_name}/{out_prefix1}_decon_means.csv")
    ###
    cell_decon_std = adata_vis.obsm['stds_cell_abundance_w_sf']
    cell_decon_std.columns = cell_types
    cell_decon_std.to_csv(f"{run_name}/{out_prefix1}_decon_std.csv")
    ###
    cell_decon_std = []
    cell_decon_means = []
    cell_decon_q95 = []
    cell_decon_q05 = []

#adata_file = f"{run_name}/{out_prefix}_sp.h5ad"
#adata_vis = sc.read(adata_file)

