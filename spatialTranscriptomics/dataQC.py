import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData

def plot_QC(adata, out_prefix):
    fig, axs = plt.subplots(2, 3, figsize=(15, 8))  # Adjusted height for better spacing

    # Histogram: Total Counts
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0, 0])
    axs[0, 0].set_title("Total Counts")

    # Histogram: Total Counts < 10,000
    sns.histplot(
        adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
        kde=False,
        bins=40,
        ax=axs[0, 1],
    )
    axs[0, 1].set_title("Total Counts < 10,000")

    # Histogram: Number of Genes by Counts
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[0, 2])
    axs[0, 2].set_title("Number of Genes by Counts")

    # Histogram: Number of Genes by Counts < 4,000
    sns.histplot(
        adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
        kde=False,
        bins=60,
        ax=axs[1, 0],
    )
    axs[1, 0].set_title("Number of Genes by Counts < 4,000")

    # Scatter Plot: Pct Counts MT vs Total Counts
    axs[1, 1].scatter(
        adata.obs["pct_counts_mt"],
        adata.obs["total_counts"],
        label="Total Counts",
        alpha=0.7,
        color="blue",
    )
    axs[1, 1].set_title("Pct MT vs Total Counts")
    axs[1, 1].set_xlabel("Pct Counts MT")
    axs[1, 1].set_ylabel("Total Counts")
    axs[1, 1].legend(loc="upper right")

    # Scatter Plot: Pct Counts MT vs Number of Genes by Counts
    axs[1, 2].scatter(
        adata.obs["pct_counts_mt"],
        adata.obs["n_genes_by_counts"],
        label="Number of Genes by Counts",
        alpha=0.7,
        color="orange",
    )
    axs[1, 2].set_title("Pct MT vs Number of Genes")
    axs[1, 2].set_xlabel("Pct Counts MT")
    axs[1, 2].set_ylabel("Number of Genes")
    axs[1, 2].legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(f"{out_prefix}_QC.png", dpi=300)  # Save as PNG with 300 DPI
    plt.show()
    plt.close(fig)
    ###


#scanpy need to read "tissue_positions_list.csv" but 10x spaceranger has different file name

def manage_soft_link(path, action="create"):
    """
    Manages a soft link for 'tissue_positions.csv' and 'tissue_positions_list.csv'.

    Args:
        path (str): The directory containing the file or link.
        action (str): The action to perform: "create" to create the link, "remove" to remove it.
    """
    source = os.path.join(path, "tissue_positions.csv")
    destination = os.path.join(path, "tissue_positions_list.csv")
    
    try:
        if action == "create":
            # Remove the destination link if it already exists
            if os.path.islink(destination) or os.path.exists(destination):
                os.unlink(destination)
            
            # Create the soft link
            os.symlink(source, destination)
            print(f"Soft link created: {destination} -> {source}")
        
        elif action == "remove":
            if os.path.islink(destination):
                os.unlink(destination)
                print(f"Soft link removed: {destination}")
            else:
                print(f"No soft link found at: {destination}")
        
        else:
            print(f"Invalid action: {action}. Use 'create' or 'remove'.")
    
    except OSError as e:
        print(f"Error managing soft link: {e}")

# Example usage
# manage_soft_link("/path/to/directory", action="create")
# manage_soft_link("/path/to/directory", action="remove")

def read_sc(sample_name, data_home = f"/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/"):

    manage_soft_link(f"{data_home}/{sample_name}/outs/spatial/")
    adata = sc.read_visium(f"{data_home}/{sample_name}/outs/")
    manage_soft_link(f"{data_home}/{sample_name}/outs/spatial/", "remove")
    
    #img = cv2.imread(f"{data_home}/{sample_name}/outs/spatial/tissue_hires_image.png")    
    #print(f"Reading {sample_images[sample_name]}")
    #img = cv2.imread(sample_images[sample_name])
    
    spatial_coord_file = f"{data_home}/{sample_name}/outs/spatial/tissue_positions.csv"
                                     #,sep=",",header=None,na_filter=False,index_col=0)
    spatial = pd.read_csv(spatial_coord_file, header = 0, index_col = 0)
    #barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres
    adata.obs["x1"]=spatial['in_tissue']
    adata.obs["x2"]=spatial['array_row']
    adata.obs["x3"]=spatial['array_col']
    adata.obs["x4"]=spatial['pxl_row_in_fullres']
    adata.obs["x5"]=spatial['pxl_col_in_fullres']
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    plot_QC(adata, sample_name + ".raw_data")
    
    sc.pp.filter_cells(adata, min_counts=1000)
    #sc.pp.filter_genes(adata, min_cells=10)
    #sc.pp.filter_cells(adata, max_counts=35000)
    plot_QC(adata, sample_name + ".filtered_data")
    print("Normalizing...")
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    return adata #, img

####

##read sample names from Data_locations_v20250516.txt
meta_data = pd.read_csv("../Data_locations_v20250516.txt", sep="\t")
sample_names = meta_data['Name'].tolist()


summary_stat = []
meta_data = []

for ss in sample_names:
    print(ss)
    adata = read_sc(ss)
    adata.obs["sample"] = ss
    n_cells = adata.shape[0]
    #
    mean_genes = np.mean(adata.obs["n_genes_by_counts"]) 
    sd_genes = np.std(adata.obs["n_genes_by_counts"])
    mean_counts = np.mean(adata.obs["total_counts"]) 
    sd_counts = np.std(adata.obs["total_counts"])  
    summary_stat.append([ss, n_cells, mean_genes, sd_genes, mean_counts, sd_counts])
    
    meta_data.append(adata.obs)
    

summary_stat = pd.DataFrame(summary_stat, columns=["Sample", "n_spots", "mean_genes", "sd_genes", "mean_counts", "sd_counts"])
summary_stat.to_csv("summary_stat.csv", index=False)

meta_data = pd.concat(meta_data, axis=0)
meta_data.to_csv("meta_data.csv", index=True)

print("Summary statistics saved to summary_stat.csv")

###plot pct_counts_mt across samples
plt.figure(figsize=(10, 6))
sns.boxplot(x="sample", y="pct_counts_mt", data=meta_data)
plt.xticks(rotation=90)
plt.xlabel("Sample")
plt.ylabel("Percentage of Counts from Mitochondrial Genes")
plt.tight_layout()
plt.savefig("pct_counts_mt_across_samples.png", dpi=300)
plt.close()


###plot number of genes per spot across samples
plt.figure(figsize=(10, 6))
sns.boxplot(x="sample", y="n_genes_by_counts", data=meta_data)
plt.xticks(rotation=90)
plt.xlabel("Sample")
plt.ylabel("Number of Genes by Counts")
plt.tight_layout()
plt.savefig("n_genes_by_counts_across_samples.png", dpi=300)
plt.close()

###plot number of spots bar plot across samples
plt.figure(figsize=(10, 6))
#reorder barplot by n_spots
summary_stat = summary_stat.sort_values(by="n_spots", ascending=False)
sns.barplot(x="Sample", y="n_spots", data=summary_stat)
plt.xticks(rotation=90)
plt.xlabel("Sample")
plt.ylabel("Number of Spots")
plt.tight_layout()
plt.savefig("n_spots_across_samples.png", dpi=300)
plt.close()