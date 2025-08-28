library(dendextend)
library(Seurat)
library(parallel)

output = "06."
#patient.list = c("ABS", "CCH", "MBF", "MBT", "RCA", "ROB", "STB", "STR", "GP27", "GP28")
patient.list = c("ABS") #, "CCH", "MBF", "MBT", "RCA", "ROB", "STB", "STR", "GP27", "GP28")

xx = lapply(patient.list, function(id){
    cat("Reading", id, "\n")
    #sc.data = Read10X(data.dir = file.path(data.home, id, "/outs/filtered_feature_bc_matrix/"))
    sc.data = t(read.table(paste0("02.macrophage.cNMF-", id, ".tsv.gz"), head =T, check.names = F, row.names =1))
    #seurat <- CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 500, project = id); rm(sc.data);
    seurat <- CreateSeuratObject(counts = sc.data, min.cells = 0, min.features = 0, project = id); rm(sc.data);
    seurat = subset(seurat, features = which(Matrix::rowSums(seurat@assays$RNA@counts >0) > 0.05*ncol(seurat)));
    seurat <- NormalizeData(seurat);
    print(seurat);
    #cor.m = cor(as.matrix(seurat[["RNA"]]@data), method = "spearman");
    cor.m = cor(as.matrix(seurat[["RNA"]]@data), method = "pearson");
    #saveRDS(cor.m, file = paste0(output, id, ".spearman.rds"))
    h = hclust(as.dist(1-cor.m), method = "average");
    ### split cluster by k 2- 8 and identify over expressed genes
    min.cluster.size = 10;
    all.clusters.n.markers = mclapply(2:ncol(cor.m), function(k){
        message(k);
        k.cut = cutree(h, k = k);
        #print(table(k.cut));
        large.cluster = which(table(k.cut) > 0.8*length(k.cut));
        small.cluster = which(table(k.cut) < min.cluster.size);
        if(max(table(k.cut)) < min.cluster.size) return(NULL)
        if(length(large.cluster) >0) k.cut[which(k.cut %in% large.cluster)] = NA;
        if(length(small.cluster) >0) k.cut[which(k.cut %in% small.cluster)] = NA;
        uniq.cluster = unique(k.cut)
        uniq.cluster = uniq.cluster[!is.na(uniq.cluster)];
        if(length(uniq.cluster) < 2) return(NULL);
        message("Find over expressed genes for clusters:", paste0(uniq.cluster, collapse = "&"))
        cluster.markers = lapply(uniq.cluster, function(k){
            k.cluster = ifelse(k.cut %in% k, k, "others")
            names(k.cluster) = names(k.cut);
            if(! identical(names(k.cluster), colnames(seurat)) ) stop("Error: sample order not identical.")
            Idents(seurat) = factor(k.cluster);
            mk =  FindMarkers(seurat, test.use = "t", only.pos = T, ident.1 = k);
            message("Cluster ", k, nrow(mk), " markers found!");
            mk;
        })
        names(cluster.markers) = paste0(k, "-", uniq.cluster);
        cluster.markers <- cluster.markers[lapply(cluster.markers,length)>0]
        res = list(k = k, cluster = k.cut, over.genes = cluster.markers);
        res;
    }, mc.cores = 10)
    all.clusters.n.markers <- all.clusters.n.markers[lapply(all.clusters.n.markers,length)>0]
    save(h, all.clusters.n.markers, file = paste0(output, id, ".rda"));
})

##
