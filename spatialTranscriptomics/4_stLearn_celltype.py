#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import pickle
import sys

i_sample = int(sys.argv[1])
n_cpu = 60

decon_result_folder = f'ref_Masahi_NG_2025/cell2location_map_n10cell/'
out_dir = f'./ref_Masahi_NG_2025_stLearn_celltype/'

info = pd.read_csv("/data/wuz6/project/25.pedHGG_spatial/03.deconvolution/cell2location/00.data.folder.csv")
#change first column name to folder
info.rename(columns={info.columns[0]: 'folder'}, inplace=True)

#exit if file exist (out_dir + sample_name+'-stLearn.h5ad')
if os.path.exists(out_dir + info['Sample'].values[i_sample] + '-stLearn.h5ad'):
    print(f"File already exists: {out_dir + info['Sample'].values[i_sample] + '-stLearn.h5ad'}")
    sys.exit(1)
else:
    print("Run ", i_sample)


import stlearn as st
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from loguru import logger


# In[3]:

#cell2location_res = pd.read_csv("01.ref_Masahi_NG_2025_cell2location.n10cell5kepoch.decon.csv.gz",index_col = 0, header=0)



# In[4]:


for i in [i_sample]:
    data_dir = info['folder'].values[i]
    sample_name = info['Sample'].values[i]
    logger.info (f"Processing: {data_dir}: {sample_name}")
    # Loading raw data #
    #sample_name = 'AG11'
    data = st.Read10X(data_dir + '/outs/')
    data.var_names_make_unique()
    #st.add.image(adata=data,
    #             imgpath=data_dir + "/outs/spatial/tissue_hires_image.png",
    #             library_id= str(sample_name), visium=True)

    # Basic normalisation #
    sc.pp.filter_cells(data, min_counts = 1000)
    st.pp.filter_genes(data, min_cells=5)
    st.pp.normalize_total(data) # NOTE: no log1p
    data.obs_names = sample_name + "_" + data.obs_names
    ###----------------------------------------------------------------------------------------
    #spot_mixtures = cell2location_res[cell2location_res['sample'] == sample_name]
    #spot_mixtures = spot_mixtures.drop('sample', axis = 1)    
    spot_mixtures = pd.read_csv(f'./{decon_result_folder}/{sample_name}_decon_q05.csv', index_col = 0, header=0)
    spot_mixtures.head()

    if not (np.all(spot_mixtures.index.values==data.obs_names.values)):
        print("Reordering spot_mixtures")
        spot_mixtures = spot_mixtures.loc[data.obs_names.to_list()]
    labels = spot_mixtures.idxmax(axis=1).to_list()
    #print(labels)
    # NOTE: using the same key in data.obs & data.uns
    data.obs['cell_type'] = labels # Adding the dominant cell type labels per spot
    data.obs['cell_type'] = data.obs['cell_type'].astype('category')
    data.uns['cell_type'] = spot_mixtures # Adding the cell type scores
    st.pl.cluster_plot(data, use_label='cell_type')
    plt.savefig(out_dir + sample_name + "-celltype.png")
    ###----------------------------------------------------------------------------------------
    ### Running the Ligand-Receptor Analysis
    # Loading the LR databases available within stlearn (from NATMI)
    lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
    #print(len(lrs))

    # Running the analysis #
    st.tl.cci.run(data, lrs,
                      min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                      distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                      n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                      n_cpus= n_cpu, # Number of CPUs for parallel. If None, detects & use all available.
                      )

    lr_info = data.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
    lr_info['sample'] = sample_name
    lr_info.to_csv(out_dir + sample_name + "-lr_info.csv")
    #print('\n', lr_info)
    #P-value adjustment
    st.tl.cci.adj_pvals(data, correct_axis='spot',
                       pval_adj_cutoff=0.05, adj_method='fdr_bh')
    ###----------------------------------------------------------------------------------------
    #Visualise the overall ranking of LRs by significant spots
    # Showing the rankings of the LR from a global and local perspective.
    # Ranking based on number of significant hotspots.
    st.pl.lr_summary(data, n_top=500, show=False)
    plt.savefig(out_dir + sample_name + "-LRsummaryRank.png")
    st.pl.lr_summary(data, n_top=50, figsize=(10,3),  show=False)
    plt.savefig(out_dir + sample_name + "-LRsummaryRankTop50.png")
    #Diagnostic plots
    st.pl.lr_diagnostics(data, figsize=(10,2.5),  show=False)
    plt.savefig(out_dir + sample_name + "-lr_diagnostics.png")
    st.pl.lr_n_spots(data, n_top=50, figsize=(11, 3),max_text=100, show=False)
    ##
    plt.savefig(out_dir + sample_name + "-sigVSnonsig_spotN_top50.png")
    st.pl.lr_n_spots(data, n_top=500, figsize=(11, 3), max_text=100, show=False)
    plt.savefig(out_dir + sample_name + "-sigVSnonsig_spotN.png")
    ##
    #LR Statistics Visualisations
    best_lr = lr_info.index.values[0] # Just choosing one of the top from lr_summary
    stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
    cmap = plt.cm.RdBu.reversed()
    reversed_cmap = cmap.reversed()
    fig, axes = plt.subplots(ncols=len(stats)+1, figsize=(30,8))
    for i, stat in enumerate(stats):
        if stat == 'lr_scores' or stat == '-log10(p_adjs)':
            st.pl.lr_result_plot(data, use_result=stat, use_lr=best_lr, cmap=cmap, show_color_bar=True, ax=axes[i])
        else:
            st.pl.lr_result_plot(data, use_result=stat, use_lr=best_lr, cmap=reversed_cmap, show_color_bar=True, ax=axes[i])
        axes[i].set_title(f'{best_lr} {stat}')
    st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=best_lr, show_color_bar=True, ax=axes[len(stats)])
    plt.savefig(out_dir + sample_name + "-bestLR-statVis.png")  
    
    ###Continuous LR coexpression for significant spots
    # Only significant spots #
    st.pl.lr_plot(data, best_lr,
              inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4,
              outer_mode='continuous', pt_scale=60,
              use_label=None, show_image=True,
              sig_spots=True)
    plt.savefig(out_dir + sample_name + "-bestLR-exp.png")  
    ###
    #Adding Arrows to show the Direction of Interaction
    #This is only useful when zooming in and want to display cell information and direction of interaction at the same time.
    st.pl.lr_plot(data, best_lr,
                  inner_size_prop=0.08, middle_size_prop=.3, outer_size_prop=.5,
                  outer_mode='binary', pt_scale=50,
                  show_image=True, arrow_width=5, arrow_head_width=2,
                  sig_spots=True, show_arrows=True)
    plt.savefig(out_dir + sample_name + "-bestLR-direction.png")  
    ##
    ###----------------------------------------------------------------------------------------
    ### Predicting significant CCIs
    #With the establishment of significant areas of LR interaction, can now determine the significantly interacting cell types.
    # Running the counting of co-occurence of cell types and LR expression hotspots #
    st.tl.cci.run_cci(data, 'cell_type', # Spot cell information either in data.obs or data.uns
                      n_cpus= n_cpu,
                      min_spots=3, # Minimum number of spots for LR to be tested.
                      spot_mixtures=True, # If True will use the label transfer scores,
                                          # so spots can have multiple cell types if score>cell_prop_cutoff
                      cell_prop_cutoff=0.2, # Spot considered to have cell type if score>0.2
                      sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                      n_perms=1000 # Permutations of cell information to get background, recommend ~1000
                     )
    ###plot CCI
    ###----------------------------------------------------------------------------------------
    #Diagnostic plot: check interaction and cell type frequency correlation
    st.pl.cci_check(data, 'cell_type', show=False)
    plt.savefig(out_dir + sample_name + "-CCI-frequence.png")  
    ### CCI Visualisations
    #With the celltype-celltype predictions completed, 
    #we implement a number of visualisations to explore the interaction landscape across LR pairs or for each independent pair.
    #CCI network
    pos_1 = st.pl.ccinet_plot(data, 'cell_type', return_pos=True)
    plt.savefig(out_dir + sample_name + "-CCI.png")
    #lrs = data.uns['lr_summary'].index.values[0:3]
    #n_spots_sig_pval: Number of spots which are significant for a given LR pair without any pval adjustment.
    #n_spots_sig: Number of spots significant for a given LR pair with pval adjustment;
        #this will be the same as n_spots_sig_pval if correct_axis=None.
    mask = lr_info['n_spots_sig'] >= 10
    lrs = lr_info.index[mask].values
    if(len(lrs) > 5):
        lrs = lrs[0:5]
    #lrs
    # Just examining the cell type interactions between selected pairs #
    #for best_lr in lrs[0:3]:
    for k in range(len(lrs)):
        best_lr = lrs[k]
        st.pl.ccinet_plot(data, 'cell_type', best_lr, min_counts=2, figsize=(6,6), pos=pos_1)
        plt.savefig(out_dir + sample_name + "-CCI_" + str(k) + "_" + best_lr + ".png")  
    #CCI chord-plot
    #The chord-plot is really useful when visualising interactions between few cell types
    st.pl.lr_chord_plot(data, 'cell_type')
    #for lr in lrs:
    for k in range(len(lrs)):
        lr = lrs[k]
        st.pl.lr_chord_plot(data, 'cell_type', lr, show=False)
        plt.savefig(out_dir + sample_name + "-CCI-chord_" + str(k) + "_" + best_lr + ".png")
    #Heatmap Visualisations
    #LR-CCI-Map
    # This will automatically select the top interacting CCIs and their respective LRs #
    #st.pl.lr_cci_map(data, 'cell_type', lrs=None, min_total=100, square_scaler= 50000, figsize=(20,10)) #show=False
    # You can also put in your own LR pairs of interest #
    st.pl.lr_cci_map(data, 'cell_type', lrs=lrs, min_total=100, square_scaler= 50000, figsize=(20,3*len(lrs))) #show=False
    plt.tight_layout()
    plt.savefig(out_dir + sample_name + "-CCI-heatmap.png")
    #CCI Maps
    #This is a heatmap equivalent to the network diagrams and chordplots, it has more quantitative benefits.
    #The # of interactions refers to the number of times a spot with the reciever cell type expressed the ligand and the source cell type expressed the receptor in the same neighbourhood.
    st.pl.cci_map(data, 'cell_type')
    #lrs = data.uns['lr_summary'].index.values[0:3]
    #for lr in lrs[0:3]:
    for k in range(len(lrs)):
        lr = lrs[k]
        st.pl.cci_map(data, 'cell_type', lr)
        plt.savefig(out_dir + sample_name + "-CCI-"+ str(k) + lr + ".png")
    ###
    ###----------------------------------------------------------------------------------------
    ###output csv
    interaction_n = []
    for k in data.uns['per_lr_cci_cell_type'].keys():
        n = data.uns['per_lr_cci_cell_type'][k]
        n['key'] = k
        interaction_n.append(n)
    df_combined = pd.concat(interaction_n, ignore_index=False)
    df_combined['sample'] = sample_name
    df_combined.to_csv(out_dir + sample_name+'-stLearn_LRcnt.csv')
    #data.write_h5ad(out_dir + sample_name+'-stLearn.h5ad') #some result cannot be saved 
    
    file_path = f'{out_dir}{sample_name}_stLearn.pkl'
    with open(file_path, "wb") as fh:    
        pickle.dump(data, fh)
    print("Done")
# In[ ]:




