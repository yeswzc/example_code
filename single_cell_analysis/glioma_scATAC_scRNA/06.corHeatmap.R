library(dendextend)
library(Seurat)
library(parallel)
library(pheatmap)

output = "06."
patient.list = c("ABS", "CCH", "MBF", "MBT", "RCA", "ROB", "STB", "STR", "GP27", "GP28")
#patient.list = c("ABS") #, "CCH", "MBF", "MBT", "RCA", "ROB", "STB", "STR", "GP27", "GP28")

xx = lapply(patient.list, function(id){
    cat("Reading", id, "\n")
    #sc.data = Read10X(data.dir = file.path(data.home, id, "/outs/filtered_feature_bc_matrix/"))
    sc.data = t(read.table(paste0("02.macrophage.cNMF-", id, ".tsv.gz"), head =T, check.names = F, row.names =1))
    #seurat <- CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 500, project = id); rm(sc.data);
    seurat <- CreateSeuratObject(counts = sc.data, min.cells = 0, min.features = 0, project = id); rm(sc.data);
    seurat = subset(seurat, features = which(Matrix::rowSums(seurat@assays$RNA@counts >0) > 0.05*ncol(seurat)));
    seurat <- NormalizeData(seurat);
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 1000)
    sel.genes = VariableFeatures(seurat)
    #seurat <- ScaleData(seurat)
    print(seurat);
    print(head(seurat[["RNA"]]@data[,1:4]))
    #cor1 = cor(as.matrix(seurat[["RNA"]]@data), method = "spearman");
    #cor2 = cor(as.matrix(seurat[["RNA"]]@data), method = "pearson");
    cor1 = cor(as.matrix(seurat[["RNA"]]@data[sel.genes,]), method = "spearman");
    cor2 = cor(as.matrix(seurat[["RNA"]]@data[sel.genes,]), method = "pearson");
    #cor1 = cor(as.matrix(seurat[["RNA"]]@scale.data), method = "spearman");
    #cor2 = cor(as.matrix(seurat[["RNA"]]@scale.data), method = "pearson");
    #saveRDS(cor.m, file = paste0(output, id, ".spearman.rds"))
    h1 = hclust(as.dist(1-cor1), method = "average");
    h2 = hclust(as.dist(1-cor2), method = "average");
    ### split cluster by k 2- 8 and identify over expressed genes
    png(paste0(output, id, ".vst1k-spearman.png"), width = 10, height = 10, unit = "in", res = 200)
    pheatmap(cor1, cluster_rows = h1, cluster_cols = h1, show_colnames = F, show_rownames = F)
    dev.off()
    png(paste0(output, id, "-vst1k.pearson.png"), width = 10, height = 10, unit = "in", res = 200)
    pheatmap(cor2, cluster_rows = h2, cluster_cols = h2, show_colnames = F, show_rownames = F)
    dev.off()
})

##
