# Mouse single-cell RNA-seq Analysis -- step 1
# Single-cell RNA-seq analysis using Scanpy, including quality control, doublet detection, clustering, marker visualization, and annotation.

# Import all necessary Python libraries for single-cell analysis, plotting, and logging.
import os
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import colorcet as cc
from loguru import logger
import warnings
warnings.filterwarnings('ignore')

# Analysis parameters
n_pc = 20
n_neighbor = 15 
'''
if len(sys.argv) < 3:
    print("Usage: python 01.scRNA_analysis.py <n_pc> <n_neighbor>")
    exit()
n_pc = 20 #int(sys.argv[1])
n_neighbor = 30 #int(sys.argv[2])
'''
logger.info(f"Using n_pc={n_pc}, n_neighbor={n_neighbor}")

logger.info(f"Using n_pc={n_pc}, n_neighbor={n_neighbor}")

# Paths and sample/tissue info
data_home = "/data/CaspiMicrobiomeData/5prime_Amy_analysis/01-Cellranger/multi/"
dataset_names = ["AZ1", "AZ2", "AZ3","AZ4", "AZ5", "AZ6", "AZ7", "AZ8"]
tissue_map = {"AZ1": "Eye", "AZ2":"SPL", "AZ3": "CLP", "AZ4": "Eye", "AZ5": "CLP", "AZ6": "Eye", "AZ7": "CLP", "AZ8": "SPF"}

def run_doublet_detection(adata, solo_threshold=0.5):
    try:
        import scrublet as scr
    except Exception:
        logger.error("Scrublet is not installed. Install it with: pip install scrublet")
        raise
    counts = adata.X.toarray() if not isinstance(adata.X, np.ndarray) and hasattr(adata.X, "toarray") else adata.X
    scrub = scr.Scrublet(counts)
    scores, preds = scrub.scrub_doublets()
    adata.obs["doublet_score"] = scores
    adata.obs["is_doublet"] = preds
    return adata

def _gene_vector(adata, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        logger.warning(f"{gene} not found in adata.var_names; filling zeros.")
        return np.zeros(adata.n_obs, dtype=float)
    vals = adata[:, gene].X
    if hasattr(vals, "toarray"):
        vals = vals.toarray()
    return np.ravel(vals)

def _violin_with_rot(adata_obj, keys=None, groupby=None, outpath='', fig_width=14, fig_height=5, rot=45, fs=8, hlines=None, hline_kwargs=None):
    sc.pl.violin(
        adata_obj,
        keys=keys,
        jitter=0.4,
        groupby=groupby,
        show=False,
    )
    fig = plt.gcf()
    fig.set_size_inches(fig_width, fig_height)
    if hlines is not None:
        defaults = dict(color="blue", linestyle="--", linewidth=1)
        if hline_kwargs:
            defaults.update(hline_kwargs)
        ys = hlines if isinstance(hlines, (list, tuple)) else [hlines]
        for ax in fig.axes:
            for y in ys:
                ax.axhline(y, **defaults)
    for ax in fig.axes:
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(rot)
            lbl.set_horizontalalignment("right")
            lbl.set_fontsize(fs)
    fig.tight_layout()
    plt.show()
    fig.savefig(outpath, bbox_inches="tight", dpi=150)    
    plt.close(fig)

def scatterplot_1(x, y, cat, xlab = '', ylab = '', title = ''):
    fig, ax = plt.subplots(figsize=(6,4))
    cat = cat.astype('category')
    codes = cat.cat.codes
    labels = list(cat.cat.categories)
    scatter = ax.scatter(
        x,
        y,
        s=3, marker='o',
        c=codes, cmap="tab10", alpha=0.3
    )
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    cmap = plt.get_cmap("tab10")
    norm = plt.Normalize(vmin=codes.min(), vmax=codes.max())
    colors = [cmap(norm(i)) for i in range(len(labels))]
    handles = [
            Line2D([0], [0], marker='o', linestyle='',
                markerfacecolor=colors[i], markeredgecolor='none',
                label=labels[i])
            for i in range(len(labels))
    ]
    ax.legend(handles=handles, title= title,
            loc='upper left', bbox_to_anchor=(1.02, 1),
            frameon=True, fontsize=7, markerscale=1.0,
            handlelength=0.6, handletextpad=0.4)
    fig.subplots_adjust(right=0.75)
    return fig, ax, scatter

# 1. Load, merge data and QC
if os.path.exists("01.scRNA_merged.h5ad"):
    adata = sc.read("01.scRNA_merged.h5ad")    
else:
    rna_list, antibody_list = [], []
    for dataset in dataset_names:
        use_folder = os.path.join(data_home, dataset, "outs/per_sample_outs/")
        sample_names = os.listdir(use_folder)
        logger.info(f"Reading samples from {dataset}, samples: {sample_names}")
        rna_list1 = []
        for sample in sample_names:
            logger.info(f"Processing sample {sample} from dataset {dataset}")
            h5_path = os.path.join(use_folder, sample, "count", "sample_filtered_feature_bc_matrix.h5")
            ad_all = sc.read_10x_h5(h5_path, gex_only=False)
            ad_all.var_names_make_unique()
            logger.info(ad_all)
            ft = ad_all.var["feature_types"].astype(str)
            rna = ad_all[:, ft == "Gene Expression"].copy()
            adt = ad_all[:, ft == "Antibody Capture"].copy()
            for obj in (rna, adt):
                obj.obs["sample"] = sample
                obj.obs["batch"] = dataset
                obj.obs["tissue"] = tissue_map[dataset]
                obj.obs_names = [f"{sample}_{x}" for x in obj.obs_names]
            rna.var["mt"] = rna.var_names.str.startswith("mt-")
            sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
            rna_list1.append(rna)
            if adt.n_vars > 0:
                antibody_list.append(adt)
        adata = anndata.concat(rna_list1, join="outer", merge="first")
        adata = run_doublet_detection(adata)
        rna_list.append(adata)
    adata_rna = anndata.concat(rna_list, join="outer", merge="first")
    del rna_list
    if len(antibody_list) > 0:
        adata_antibody = anndata.concat(antibody_list, join="outer", merge="first")
        from scipy import sparse
        Xp = adata_antibody.X.toarray() if sparse.issparse(adata_antibody.X) else adata_antibody.X
        adata_rna.obsm["Antibody Capture"] = pd.DataFrame(Xp, index=adata_antibody.obs_names, columns=adata_antibody.var_names)
        del antibody_list
    del adata_antibody
    adata = adata_rna
    plt.figure(figsize=(8, 4))
    sns.kdeplot(adata.obs['pct_counts_mt'], fill=True, alpha=0.3, cut=0, clip=(0, None), color="gray")
    plt.axvline(5, color='blue', linestyle='--')
    plt.axvline(8, color='black', linestyle='--')
    plt.axvline(10, color='red', linestyle='--')
    plt.xlabel('Percent Mitochondrial Counts')
    plt.ylabel('Density')
    plt.legend(title='Sample')
    plt.xlim(left=0)
    plt.tight_layout()
    plt.savefig("figures/01_pct_counts_mt_density.png", dpi=150)
    plt.show()
    plt.close()
    _violin_with_rot(adata, keys=["n_genes_by_counts"], groupby="sample", outpath="figures/01_filtered_n_genes_violin.png", rot=90, hlines=None)
    _violin_with_rot(adata, keys=["total_counts"], groupby="sample", outpath="figures/01_filtered_total_counts_violin.png", rot=90, hlines=[500, 8000])
    _violin_with_rot(adata, keys=["pct_counts_mt"], groupby="sample", outpath="figures/01_filtered_qc_violin.png", rot=90, hlines=5)
    fig, ax, scatter = scatterplot_1(adata.obs["total_counts"], adata.obs["n_genes_by_counts"], adata.obs["sample"], xlab='Total Counts', ylab='Number of Genes', title='Sample')
    ax.axhline(8000, color='red', linestyle='--')    
    fig.show()
    fig.savefig(f"figures/01_QC_totalCounts-nGene.png")
    fig, ax, scatter = scatterplot_1(adata.obs["n_genes_by_counts"], adata.obs["pct_counts_mt"], adata.obs["sample"], xlab = 'Total # Genes', ylab = 'Percent Mitochondrial Counts', title = 'Sample')
    ax.axhline(5, color='blue', linestyle='--')
    ax.axvline(500, color='red', linestyle='--')
    ax.axvline(8000, color='red', linestyle='--')
    fig.show()
    fig.savefig(f"figures/01_QC_MT-nGene.png")
    cell_counts_raw = adata.obs['sample'].value_counts()
    adata = adata[adata.obs["is_doublet"] == False, :]
    adata = adata[(adata.obs["pct_counts_mt"] < 5) & (adata.obs["n_genes_by_counts"] > 500) & (adata.obs["n_genes_by_counts"] < 8000)]
    cell_counts = adata.obs['sample'].value_counts()
    cell_counts = pd.DataFrame({
        "sample": cell_counts_raw.index,
        "raw": cell_counts_raw.values,
        "filtered": cell_counts.reindex(cell_counts_raw.index).values
    })
    cell_counts.to_csv("figures/01_cell_counts_per_sample.csv", index=False)
    print(cell_counts)
    sc.pl.violin(adata, keys="n_genes_by_counts", groupby="sample", rotation=90, jitter=0.4, save="01_filtered_genes_per_cell_violin.png")
    adata.write("01.scRNA_merged.h5ad")

fig, ax, scatter = scatterplot_1(adata.obs["total_counts"], adata.obs["n_genes_by_counts"], adata.obs["sample"], xlab='Total Counts', ylab='Number of Genes', title='Sample')
fig.show()
fig.savefig(f"figures/01_filtered_totalCounts-nGene.png")
fig, ax, scatter = scatterplot_1(adata.obs["n_genes_by_counts"], adata.obs["pct_counts_mt"], adata.obs["sample"], xlab='Total # Genes', ylab='Percent Mitochondrial Counts', title='Sample')
fig.show()
fig.savefig(f"figures/01_filtered_MT-nGene.png")

# Normalization and Variable Gene Selection
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, batch_key="batch")
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

sc.pl.violin(adata, keys=["Cd4", "Cd8a", "Cd3e"], jitter=0.4, save="01_all_cell_check_Cd4cell_violin.png")

adata.obs['group'] = adata.obs['sample'].str[:3]
adata.obs['group2'] = adata.obs["group"].astype(str) + "-" + adata.obs["tissue"].astype(str)

if True:
    adata1 = adata[:, adata.var["highly_variable"]]
    sc.pp.scale(adata1)
    sc.tl.pca(adata1, svd_solver="arpack")
    sc.pl.pca_variance_ratio(adata1, n_pcs=50, save="_elbowplot.png")
    adata.obsm["X_pca"] = adata1.obsm["X_pca"]
    adata.uns["pca"] = adata1.uns["pca"]
    del adata1
else:
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, save="_all_genes_elbowplot.png")

# Add SingleR Annotation generated by R
singleR_results = pd.read_csv("01.allSample_scRNA_singleRImmGenAnn.csv", index_col=0)
adata.obs = adata.obs.join(singleR_results["pruned.labels"], how="left")
adata.obs["sample2"] = adata.obs["batch"].astype(str) + "-" + adata.obs["sample"].astype(str)

tableau10 = plt.get_cmap("tab10").colors
tableau20 = plt.get_cmap("tab20").colors
color35 = palette = cc.glasbey_light[:35]
do_batch_correction = True

if False:
    for n_pc in [20, 30]:
        #n_pc = 30 #testing 20,30
        #n_neighbor = 10  #testing 3, 5,10, 15,30, 50 #for bbkn, large n_neighbor is very slow
        n_neighbor_to_test = [15,30,50]
        if do_batch_correction:
            n_neighbor_to_test = [3,5,10,15]
            #With BBKNN, each cell will have
            # neighbors_within_batch × number_of_batches (8)
            # neighbors in total
        for n_neighbor in n_neighbor_to_test:
            logger.info(f"\n\nRunning do_batch_correction={do_batch_correction}, n_pc={n_pc}, n_neighbor={n_neighbor}")
            if do_batch_correction:
                #neighbors_within_batch default is 3
                sc.external.pp.bbknn(adata, batch_key="batch", n_pcs = n_pc, neighbors_within_batch = n_neighbor)
            else:
                sc.pp.neighbors(adata, n_neighbors=n_neighbor, n_pcs=n_pc)
            sc.tl.leiden(adata, resolution=1)
            sc.tl.umap(adata, random_state=42)

            fig, axes = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)
            sc.pl.umap(adata, color="sample", palette=tableau20, ax=axes[0, 0], show=False, title="sample")
            sc.pl.umap(adata, color="tissue", palette=tableau10, ax=axes[0, 1], show=False, title="tissue")
            sc.pl.umap(adata, color="leiden", palette=color35, ax=axes[1, 0], show=False, title="cluster", legend_loc='on data')
            sc.pl.umap(adata, color="leiden", palette=color35, ax=axes[1, 1], show=False, title="cluster")
            sc.pl.umap(adata, color="pruned.labels", palette=tableau20, ax=axes[2, 0], show=False, title="singleR", legend_loc='on data')
            sc.pl.umap(adata, color="pruned.labels", palette=tableau20, ax=axes[2, 1], show=False, title="singleR")
            for ax in axes.flatten():
                ax.set_box_aspect(aspect=1)
                legend = ax.get_legend()
                if legend is not None:
                    legend.set_bbox_to_anchor((1.0, 0.5))
                    legend._legend_box.align = "left"
            if do_batch_correction:
                outname = f"figures/01-0_BBKNumapPC{n_pc}_neighbor{n_neighbor}_0sample_batch_tissue.png"
            else:
                outname = f"figures/01-0_umapPC{n_pc}_neighbor{n_neighbor}_0sample_batch_tissue.png"
            plt.show()
            fig.savefig(outname, bbox_inches="tight", dpi=150)
            plt.close(fig)



if True:
    n_pc = 30
    n_neighbor = 5
    do_batch_correction = True
    logger.info(f"Running do_batch_correction={do_batch_correction}, n_pc={n_pc}, n_neighbor={n_neighbor}")
    if do_batch_correction:
        print(adata.obs["batch"].value_counts())
        sc.external.pp.bbknn(adata, batch_key="batch", n_pcs = n_pc, neighbors_within_batch = n_neighbor)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbor, n_pcs=n_pc)
    sc.tl.leiden(adata, resolution=1)
    sc.tl.umap(adata, random_state=42)
    fig, axes = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)
    sc.pl.umap(adata, color="sample", palette=tableau20, ax=axes[0, 0], show=False, title="sample")
    sc.pl.umap(adata, color="tissue", palette=tableau10, ax=axes[0, 1], show=False, title="tissue")
    sc.pl.umap(adata, color="leiden", palette=color35, ax=axes[1, 0], show=False, title="cluster", legend_loc='on data')
    sc.pl.umap(adata, color="leiden", palette=color35, ax=axes[1, 1], show=False, title="cluster")
    sc.pl.umap(adata, color="pruned.labels", palette=tableau20, ax=axes[2, 0], show=False, title="singleR", legend_loc='on data')
    sc.pl.umap(adata, color="pruned.labels", palette=tableau20, ax=axes[2, 1], show=False, title="singleR")
    for ax in axes.flatten():
        ax.set_box_aspect(aspect=1)
        legend = ax.get_legend()
        if legend is not None:
            legend.set_bbox_to_anchor((1.0, 0.5))
            legend._legend_box.align = "left"
    if do_batch_correction:
        outname = f"figures/01-0_BBKNumapPC{n_pc}_neighbor{n_neighbor}_0sample_batch_tissue.png"
    else:
        outname = f"figures/01-0_umapPC{n_pc}_neighbor{n_neighbor}_0sample_batch_tissue.png"
    plt.show()
    fig.savefig(outname, bbox_inches="tight", dpi=150)
    plt.close(fig)

# Antibody Capture Visualization
def plot_antibody_capture(adata, n_pc=30, n_neighbor=50, output_figure=""):
    """Plot UMAPs for antibody capture features if present in AnnData.obsm."""
    if output_figure == "":
        logger.warning("No output figure name provided.")
        return
    if "Antibody Capture" not in adata.obsm_keys():
        logger.warning("No Antibody Capture data found in adata.obsm.")
        return
    prot_raw = adata.obsm["Antibody Capture"]
    prot_raw = prot_raw.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(np.float32)
    log1p = np.log1p(prot_raw)
    prot_clr = log1p.sub(log1p.mean(axis=1), axis=0)
    adata.obsm["Antibody Capture (CLR)"] = prot_clr
    ab_names = list(prot_clr.columns)
    adata_adt = anndata.AnnData(
        X=prot_clr.values,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=ab_names),
    )
    adata_adt.var_names_make_unique()
    if "X_umap" in adata.obsm:
        adata_adt.obsm["X_umap"] = adata.obsm["X_umap"].copy()
    fig = sc.pl.umap(
        adata_adt,
        color=ab_names,
        ncols=6,
        show=False,
        return_fig=True,
        color_map="rocket",
    )    
    fig.savefig(output_figure, bbox_inches="tight", dpi=150)
    plt.close(fig)

# Marker genes (example: T cell, B cell, etc.)
t_cell_markers = ["Cd3e", "Cd3d", "Cd3g", "Cd2", "Cd247", "Tcrb"]
cd4_t_markers = ["Cd3e", "Cd4", "Il7r", "Cd27"]
cd8_t_markers = ["Cd3e", "Cd8a", "Cd8b1", "Gzmb", "Prf1"]
naive_t_markers = ["Cd3e", "Cd4", "Cd8a", "Sell", "Ccr7", "Lef1"]
effector_t_markers = ["Cd3e", "Cd4", "Cd8a", "Gzmb", "Prf1", "Ifng", "Cd44"]
memory_t_markers = ["Cd3e", "Cd4", "Cd8a", "Cd44", "Sell"]
treg_markers = ["Cd3e", "Cd4", "Il2ra", "Foxp3", "Ctla4", "Ikzf2"]
tex_markers = ["Pdcd1", "Lag3", "Havcr2", "Tox", "Ctla4"]
cytotoxic_t_markers = ["Cd8a", "Gzmb", "Gzma", "Prf1", "Ifng"]
nk_cell_markers = ["Ncr1", "Klra1", "Klrk1", "Nkg7", "Gnly", "Gzmb", "Klrb1a","Klrb1b", "Klrb1c", "Klrb1f", "Klrb1"]
ilc_markers = ["Id2", "Il7r", "Gata3", "Rora", "Tcf7"]
ilc1_markers = ["Tbx21", "Ifng", "Ncr1", "Cd122", "Cd127"]
ilc2_markers = ["Gata3", "Il1rl1", "Il5", "Il13", "Areg", "Rora"]
ilc3_markers = ["Rorc", "Il22", "Il17a", "Cd117", "Ahr"]
b_cell_markers = ["Cd19", "Cd79a", "Cd79b", "Ms4a1", "Ighm", "Ighd"]
b_early_immature = ["Cd34", "Vpreb1", "Vpreb3", "Igll1", "Dntt"]
b_mature_naive = ["Ighm", "Igkc", "Iglc1", "Cd24a", "Fcrl1"]
b_germinal_center = ["Bcl6", "Mki67", "Aicda"]
b_plasma = ["Prdm1", "Irf4", "Xbp1", "Sdc1", "Jchain"]
b_memory = ["Cd27", "Tnfrsf13b", "Tnfrsf13c"]
b_negative_selection = ["Cd3e", "Cd14", "Itgam", "Adgre1", "Nkg7", "Klrb1c", "Epcam", "Krt8"]
plasma_cell_markers = ["Sdc1", "Prdm1", "Xbp1", "Ighg1", "Ighg2b"]
mait_markers = ["Trav1", "Mr1", "Cd3e", "Cd8a", "Cd4", "Zbtb16", "Il18r1", "Il7r", "Rorc", "Tbx21", "Ifng", "Il17a", "Klrk1", "Cd69", "Maf", "Ccr6", "Ccr5"]
iel_markers = ["Cd8a", "Cd8b1", "Cd3e", "Itgae", "Cd103", "Tcrd", "Tcrg", "Cd69", "Cd122"]
monocyte_markers = ["Ly6c2", "Cd14", "Fcgr3", "Ccr2", "Itgam", "Cd11b", "Lyz2"]
monocyte_nonclassical = ["Ly6c2", "Cx3cr1", "Nr4a1"]
macrophage_markers = ["Adgre1", "Cd68", "Cd14", "Cd64", "Lyz2", "Mertk", "Maf", "Csf1r"]
m1_macrophage_markers = ["Nos2", "Tnf", "Il6", "Cd86"]
m2_macrophage_markers = ["Arg1", "Mrc1", "Retnla", "Cd206", "Chil3"]
dc_common_markers = ["Itgax", "Cd11c", "H2-Ab1", "Cd74"]
cDC1_markers = ["Xcr1", "Batf3", "Clec9a", "Cd8a", "Irf8"]
cDC2_markers = ["Cd11b", "Sirpa", "Irf4", "Clec10a"]
pDC_markers = ["Siglech", "Bst2", "Ly6d", "Tcf4", "Irf7", "B220"]
neutrophil_markers = ["Ly6g", "S100a8", "S100a9", "Cxcr2", "Itgam", "Fcgr3", "Elane", "Mpo", "Csf3r"]
eosinophil_markers = ["SiglecF", "Ccr3", "Prg2", "Il5ra", "Gata1", "Epx", "Ccl24"]
basophil_markers = ["Mcpt8", "Cd200r3", "Fcer1a", "Il4", "Il13", "Prss34", "Kit"]
hsc_markers = ["Kit", "Cd34", "Flt3", "Ly6a", "Cd150", "Cd48", "Procr", "Mpl", "Mecom", "Runx1", "Gata2"]
fibroblast_markers = ["Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra", "Fbln1", "Thy1", "Mcam"]
epithelial_markers = ["Epcam", "Krt8", "Krt18", "Krt19", "Krt7", "Ocln", "Cldn3", "Cldn4", "Cdh1"]
endothelial_markers = ["Pecam1", "Cdh5", "Kdr", "Flt1", "Tie1", "Tek", "Vwf"]
mouse_immune_markers = {
    "T_Cell": t_cell_markers,
    "CD4_T": cd4_t_markers,
    "CD8_T": cd8_t_markers,
    "Naive_T": naive_t_markers,
    "Effector_T": effector_t_markers,
    "Memory_T": memory_t_markers,
    "Treg": treg_markers,
    "Tex": tex_markers,
    "Cytotoxic_T": cytotoxic_t_markers,
    "B_Cell": b_cell_markers,
    "B_Early_Immature": b_early_immature,
    "B_Mature_Naive": b_mature_naive,
    "B_Germinal_Center": b_germinal_center,
    "B_Plasma": b_plasma,
    "B_Memory": b_memory,
    "B_Negative_Selection": b_negative_selection,
    "Plasma_Cell": plasma_cell_markers,
    "NK_Cell": nk_cell_markers,
    "ILC": ilc_markers,
    "ILC1": ilc1_markers,
    "ILC2": ilc2_markers,
    "ILC3": ilc3_markers,
    "MAIT": mait_markers,
    "IEL": iel_markers,
    "Monocyte": monocyte_markers,
    "Monocyte_Nonclassical": monocyte_nonclassical,
    "Macrophage": macrophage_markers,
    "M1_Macrophage": m1_macrophage_markers,
    "M2_Macrophage": m2_macrophage_markers,
    "DC_Common": dc_common_markers,
    "cDC1": cDC1_markers,
    "cDC2": cDC2_markers,
    "pDC": pDC_markers,
    "Neutrophil": neutrophil_markers,
    "Eosinophil": eosinophil_markers,
    "Basophil": basophil_markers,
    "HSC": hsc_markers,
    'Fibroblast': fibroblast_markers,
    'Epithelial': epithelial_markers,
    'Endothelial': endothelial_markers
}
immune_markers = [marker for sublist in mouse_immune_markers.values() for marker in sublist]
immune_markers = list(dict.fromkeys(immune_markers))
immune_markers_present = [g for g in immune_markers if g in adata.var_names]
missing = sorted(set(immune_markers) - set(immune_markers_present))
if missing:
    logger.warning(f"{len(missing)} immune markers not found and will be skipped (e.g., {', '.join(missing)})")

fig = sc.pl.umap(adata, color=immune_markers_present, ncols=15, color_map='rocket', show=False, return_fig=True)
outname = f"figures/01-4_umapPC{n_pc}_neighbor{n_neighbor}_immuneMarkers.png"
fig.savefig(outname, bbox_inches="tight", dpi=150)
plt.close(fig)

for cell_type_marker in mouse_immune_markers.keys():
    markers = mouse_immune_markers[cell_type_marker]
    markers_present = [g for g in markers if g in adata.var_names]
    fig = sc.pl.umap(adata, color=markers_present, ncols=5, color_map='rocket', show=False, return_fig=True)
    outname = f"figures/01-3_umapPC{n_pc}_neighbor{n_neighbor}_{cell_type_marker}.png"
    fig.savefig(outname, bbox_inches="tight", dpi=150)
    plt.close(fig)

x_table = pd.crosstab(adata.obs["leiden"], adata.obs["pruned.labels"])
import seaborn as sns
sns.set(font_scale=0.7)
sns.heatmap(x_table, cmap="viridis", annot=True, fmt="d")
plt.savefig("figures/01-5_cluster_singleR-heatmap.png", dpi=300, bbox_inches="tight")

adata.write(f"01.scRNA_annotated_step1_PC{n_pc}_neighbor{n_neighbor}.h5ad")

# Assign cell types to each cluster
cluster_to_cell_type = {
    # Fill in cluster to cell type mapping
}
missing_clusters = set(adata.obs["leiden"].cat.categories) - set(cluster_to_cell_type.keys())
if missing_clusters:
    logger.warning(f"Missing cluster annotations for clusters: {missing_clusters}")
adata.obs["cell_type1"] = adata.obs["leiden"].map(cluster_to_cell_type).astype("category")
adata.obs["cell_type1"].value_counts()

plt.rcdefaults()
fig = sc.pl.umap(adata, color="cell_type1", palette=tableau20, legend_loc='on data', title='cell_type', show=False, return_fig=True)
plt.savefig(f"figures/01-6_umapPC{n_pc}_neighbor{n_neighbor}_cell_type1.png", bbox_inches="tight", dpi=150)
plt.show()
plt.close(fig)

fig = sc.pl.umap(adata, color="total_counts", show=False, return_fig=True)
plt.savefig(f"figures/01-7_umapPC{n_pc}_neighbor{n_neighbor}_n_counts.png", bbox_inches="tight", dpi=150)
plt.show()
plt.close(fig)

print("scanpy:", sc.__version__)
print("anndata:", anndata.__version__)
print("pandas:", pd.__version__)
print("numpy:", np.__version__)

exit()
