#library(celldex)
#library(SingleR)
#library(Seurat)
#library(scmap)
#library(scater)
library(SingleCellExperiment)
#output = "single"

meta <- read.table("/data/wuz6/data/single_cell/expression/IDHmut_A_glioma.Venteicher.Science.2017.SmartSeq/IDH_A_cell_type_assignment_portal_v2.txt", head =T, sep = "\t");
meta$NAME <- gsub(" ", "", meta$NAME);
cnt <- read.table("/data/wuz6/data/single_cell/expression/IDHmut_A_glioma.Venteicher.Science.2017.SmartSeq/IDH_A_processed_data_portal.txt.gz", head =T, row.names = 1, check.names = F)
cnt <- cnt[,colnames(cnt) %in% meta$NAME]
meta = meta[match(colnames(cnt), meta$NAME),]
rownames(meta) <- meta$NAME
#needs to be integer
ref <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=round(cnt)), colData = DataFrame(label = meta$CLUSTER), rowData = DataFrame(name = rownames(cnt)))
ref$cell.type <- meta$CLUSTER
colData(ref)$cell_type1 <- meta$CLUSTER
set.seed(12345)
#O IDH: FNF, GPU, RPO, 
#A IDH: CEH, JAS, JWS, SAB, GP29, GP30 #DBG, VNA
#ref <- celldex::HumanPrimaryCellAtlasData()
#colData(ref)$cell_type1 <- colData(ref)$label.fine
#colData(ref)$cell_type1 <- colData(ref)$label.main
#ref <- scRNAseq::DarmanisBrainData()
#colData(ref)$cell_type1 <- colData(ref)$cell.type


rowData(ref)$feature_symbol <- rownames(ref)
#ref <- scater::logNormCounts(ref)
#ref_sce <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=Matrix::Matrix(assays(ref)$logcounts)), 
#                                                      colData=colData(ref), rowData=rowData(ref))
ref_sce <- ref
png("ref-IDH.scmap.featureSel.png", width = 7, height = 7, units = "in", res = 300) 
ref_sce <- scmap::selectFeatures(ref_sce, suppress_plot=F)
dev.off()

#ref.seurat <- Seurat::CreateSeuratObject(counts = assays(ref)$logcounts)
#ref.seurat@assays$data = ref.seurat@assays$RNA
#ref.seurat <- Seurat::FindVariableFeatures(ref.seurat, selection.method = "vst", nfeatures = 2000)
#features.sel <- Seurat::VariableFeatures(ref.seurat);rm(ref.seurat)
#Seurat::Idents(ref.seurat) <-  colData(ref)$label.main
#mk <- Seurat::FindAllMarkers(ref.seurat, only.pos = TRUE,min.cells.group = 1)
#mk1 <- mk %>% dplyr::group_by(cluster) %>% dplyr::slice_max(n = 40, order_by = avg_log2FC)
#scmap::selectFeatures work with only data with missing (==0)

#rowData(ref_sce)$scmap_features <- FALSE
#rowData(ref_sce)$scmap_features[rownames(ref_sce) %in% features.sel] <- TRUE

# Inspect the first 50 genes selected by scmap
#rownames(ref_sce)[which(rowData(ref_sce)$scmap_features)][1:50]

# Create a list of key markers that you want to use
#my_key_markers = c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "IGKC")
# Ensure markers are in the list of features used by scmap
#rowData(ref_sce)$scmap_features[rownames(ref_sce) %in% my_key_markers] <- TRUE
# You can check and see if this added any genes by checking the length 
# of the vector of gene names again
#length(rownames(ref_sce)[which(rowData(ref_sce)$scmap_features)])
mt_genes <- rownames(ref_sce)[grep("^MT-", rownames(ref_sce))]
# Remove these genes from the features used by scmap
rowData(ref_sce)$scmap_features[rownames(ref_sce) %in% mt_genes] <- FALSE
# Check how many genes this is
#length(rownames(ref_sce)[which(rowData(ref_sce)$scmap_features)])
# Extract the features and assign them to a new variable, "scmap_feature_genes"
scmap_feature_genes <- rownames(ref_sce)[which(rowData(ref_sce)$scmap_features)]
# Note that the number of genes/features is identical to what we just checked
message("** N=", length(scmap_feature_genes), " genes/features selected.")
ref_sce <- scmap::indexCluster(ref_sce)

# Store expression information as a variable
scmap_cluster_reference <- metadata(ref_sce)$scmap_cluster_index
rm(ref_sce)
# Update the previous reference to also contain the 'scmap-cell' reference
#ref_sce <- scmap::indexCell(ref_sce)
# Extract the scmap index from the reference and store as a variable
#scmap_cell_reference <- metadata(ref_sce)$scmap_cell_index
# Extract the associated cell IDs from the reference and save as a variable
#scmap_cell_metadata <- colData(ref_sce)




cell_annotate <- function(id, data.path, ref, scmap_cluster_reference){
    message("Annotate ", id)
    sc.data = Seurat::Read10X(data.dir = paste0(data.path, id, "/outs/filtered_feature_bc_matrix/"))
    scRNA <- Seurat::CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 500, project = id); rm(sc.data);
    scRNA[["percent.mt"]] <- Seurat::PercentageFeatureSet(scRNA, pattern = "^MT-")
    scRNA = subset(scRNA, features = which(Matrix::rowSums(scRNA@assays$RNA@counts) > 10));
    scRNA <- subset(scRNA, subset = nFeature_RNA > 500 & percent.mt < 10)
    scRNA <- Seurat::NormalizeData(scRNA,  verbose = F)
    
    
    #run seurat umap
    scRNA <- Seurat::FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000,  verbose = F)
    scRNA <- Seurat::ScaleData(scRNA,  verbose = F)
    scRNA <- Seurat::RunPCA(scRNA, features = Seurat::VariableFeatures(object = scRNA), verbose = F)
    scRNA <- Seurat::RunUMAP(scRNA, dims = 1:15, verbose = F)
    
    res = cbind(scRNA[[]], scRNA@reductions$umap@cell.embeddings)
    #print(head(res))
    query_sce <- Seurat::as.SingleCellExperiment(scRNA, assay = 'RNA')#
    # normalize the data using the scater package
    #query_sce <- scater::logNormCounts(query_sce)
    rowData(query_sce)$feature_symbol <- rownames(scRNA)
    rm(scRNA)
    # Run scmapCluster
    scmap_cluster_res <- scmap::scmapCluster(projection=query_sce, 
                                         index_list=list(hpca = scmap_cluster_reference), 
                                         threshold=0.1)

    # plot the results of our annotation
    png(paste0(id, ".scmap.barplot.png"), width = 10, height = 5, units = "in", res = 150)
    par(mar=c(13, 4, 1, 0))
    barplot(table(scmap_cluster_res$combined_labs), las=2)
    dev.off()


    # Store this annotation information within the query object
    colData(query_sce)$scmap_cluster <- scmap_cluster_res$combined_labs
    res$scmap_cluster <- scmap_cluster_res$combined_labs

    # Make a UMAP of the cells, labeled with the cell-type annotations from scmapCluster
    #query_sce <- scater::runUMAP(query_sce)
    p <- ggplot2::ggplot(res, ggplot2::aes(x= UMAP_1, y = UMAP_2, col = scmap_cluster)) + ggplot2::geom_point(size = 1) + ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1)
    png(paste0(id, ".scmap.umap.png"), width = 15, height = 5, units = "in", res = 150)
    #print(scater::plotReducedDim(query_sce, dimred="UMAP", colour_by="scmap_cluster")+ggplot2::theme(aspect.ratio = 1))
    print(p)
    dev.off()

    # Run SingleR on the query data and the reference to acquire
    # cell-type predictions for the cells in the query dataset
    #predictions <- SingleR::SingleR(test=query_sce, ref=ref, labels=ref$label.main)
    predictions <- SingleR::SingleR(test=query_sce, ref=ref, labels=ref$cell.type)
    
    # Change NAs to "ambiguous"
    predictions$pruned.labels[which(is.na(predictions$pruned.labels))] <- "unassigned"
    # Add singleR labels to query_sce
    colData(query_sce)$singleR <- predictions$pruned.labels
    res$singleR <- predictions$pruned.labels
    # Create a bar plot of number of cells per assigned cell ID
    png(paste0(id, ".singleR.barplot.png"), width = 10, height = 5, units = "in", res = 150)
    par(mar=c(13, 4, 2, 0))
    barplot(table(predictions$pruned.labels), las=2)
    dev.off()
    rm(p)
    p <- ggplot2::ggplot(res, ggplot2::aes(x= UMAP_1, y = UMAP_2, col = singleR)) + ggplot2::geom_point(size = 1) + ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1)
    png(paste0(id, ".singleR.umap.png"), width = 15, height = 5, units = "in", res = 150)
    #print(scater::plotReducedDim(query_sce, dimred="UMAP", colour_by="singleR") + ggplot2::theme(aspect.ratio = 1))
    print(p)
    dev.off()
    #query_sce;
    rownames(res) <- paste0(res$orig.ident, ".", rownames(res))
    res
}


data.folder = "/data/wuz6/data/aldape_lab/scRNA/cellranger/"
patient.list = c("CEH", "JAS", "JWS", "SAB", "GP29", "GP30", "GP32" );

all.ann <- lapply(patient.list, function(ID){
    res <- cell_annotate(id = ID, data.path = data.folder, ref = ref, scmap_cluster_reference = scmap_cluster_reference)
    res;
})


save(all.ann, file = "IDH-A.ann.rda")
