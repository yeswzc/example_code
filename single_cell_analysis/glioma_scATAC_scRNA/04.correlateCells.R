library(Seurat)


output = "04"
patient.list = c("ABS", "CCH", "MBF", "MBT", "RCA", "ROB", "STB", "STR", "GP27", "GP28")

feature.sel.by.patient = TRUE;
n.feature.sel.by.patient = 3000
do.sample.wise.cluster = FALSE;

data.list = lapply(patient.list, function(id){
    cat("Reading", id, "\n")
    #sc.data = Read10X(data.dir = file.path(data.home, id, "/outs/filtered_feature_bc_matrix/"))
    sc.data = t(read.table(paste0("02.macrophage.cNMF-", id, ".tsv.gz"), head =T, check.names = F, row.names =1))
    #seurat <- CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 500, project = id); rm(sc.data);
    seurat <- CreateSeuratObject(counts = sc.data, min.cells = 0, min.features = 0, project = id); rm(sc.data);
    print(seurat);
    return(seurat);
})
##
cat("Combine all samples\n")
scRNA = merge(data.list[[1]], y = data.list[-1], add.cell.ids = patient.list, project = "GBM")
#scRNA = subset(scRNA, features = which(Matrix::rowSums(scRNA@assays$RNA@counts >0) > 6));
rm(data.list);
scRNA;

scRNA <- NormalizeData(scRNA);
#scRNA <- NormalizeData(scRNA, normalization.method = "RC");
if( feature.sel.by.patient){
    cat("Sample-wise feature selection...\n");
    #scRNA <- NormalizeData(scRNA);
    feature.sel = lapply(unique(scRNA[[]][,1]), function(patient){
            cat("Find features for ", patient, "\n")
            idx = which(scRNA[[]][,1] == patient)
            tmp = subset(scRNA, cells = idx);
            tmp = FindVariableFeatures(tmp, selection.method = "vst", nfeatures = n.feature.sel.by.patient, verbose =F)
            if(do.sample.wise.cluster){
                #tmp <- SCTransform(tmp, verbose = T, variable.features.n = 2000);
                tmp <- ScaleData(tmp)
                tmp <- RunPCA(tmp, verbose = FALSE);
                tmp <- RunTSNE(tmp, dims = 1:dim.use, verbose = FALSE)
                tmp <- RunUMAP(tmp, dims = 1:dim.use, verbose = FALSE)
                tmp <- FindNeighbors(tmp, dims = 1:dim.use, verbose = FALSE)
                tmp <- FindClusters(tmp, verbose = FALSE, resolution = seq(0.5, 1.5, 0.1));
                res = cbind(tmp[[]], tmp@reductions$tsne@cell.embeddings, tmp@reductions$umap@cell.embeddings, tmp@reductions$pca@cell.embeddings[,1:30])
                saveRDS(res, file = paste0(output, "-", patient, ".clustering.rds"));
            }
            #tmp = SCTransform(scRNA, verbose = FALSE, variable.features.n = 2000)
            print(length(VariableFeatures(tmp)));
            VariableFeatures(tmp);
        });
        #feature.selected = Reduce(intersect, feature.sel)
        feature.selected = unlist(feature.sel)
        feature.selected =  names(which(table(feature.selected) >= 3)); #2 or 3?
        cat("N = ", length(feature.selected), "features selected.\n")
        VariableFeatures(scRNA) = feature.selected;
        #scRNA <- ScaleData(scRNA);
        #scRNA <- RunPCA(scRNA, verbose = FALSE, features = feature.selected);
}else{
        cat("Non sample-wise feature selection using SCT...\n");
        #scRNA <- ScaleData(scRNA); #will only scale variable features
        scRNA <- SCTransform(scRNA, verbose = T, variable.features.n = 3000); #variable.features.n default is 3000
        #scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA));
};


cor.m = cor( as.matrix(scRNA[["RNA"]]@data), method = "spearman");
print(head(cor.m[,1:4]))
saveRDS(cor.m, file = "04.sc.cor.rds");
