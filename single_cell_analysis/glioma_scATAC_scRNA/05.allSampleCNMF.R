library(Seurat)
#library(SeuratDisk)

n.core = 10
n.genes = 1000
nmf.cmd = "python /data/wuz6/software/cNMF/cnmf.py"
#
input = "01.inhouse.GBM.scRNA.feature1500.SeuratObjCluster.rda" 
#"02.inhouse.GBM.scRNA.macrophage.seuratObj.rda"
output = "05.cNMF.selFeatures"
k.nmf = 2:12
k.nmf = paste(k.nmf, collapse = " ")

macrophage.markers = c("AIF1", "CD14", "CD16", "CD64", "CD68", "CD71", "CCR5", "CD163", "ITGAM");
macroglia.vs.macrophage = c("ITGA4") #ITGA4 is CD49D, high CD49D is macrophage, low CD49D is microglia

sample.list = c("ABS", "CCH", "MBF", "MBT", "RCA", "ROB", "STB", "STR", "GP27", "GP28")

macrophage.cluster = list("ABS" = c(2,3,10,12),
                          "CCH" = c(4),
                          "MBF" = c(0, 1, 3, 4, 9, 11),
                          "MBT" = c(1, 4, 5),
                          "RCA" = c(3),
                          "ROB" = c(3,5),
                          "STB" = c(0,1,2,7,14),
                          "STR" = c(0, 1, 16),
                          "GP27" = c(2,8),
                          "GP28" = c(2, 5, 12))
if(1){
file.list = lapply(sample.list, function(ID){

    keep.clusters = macrophage.cluster[[ID]];
    load(paste0("01.GBM-", ID, ".SeuratObjCluster.rda"));
    ##select cells of interest
    Idents(scRNA) = factor(scRNA[[]]$SCT_snn_res.0.6) 
    DefaultAssay(scRNA) <- "RNA"
    scRNA[['SCT']] <- NULL
    #Idents(scRNA) = scRNA[[]]$RNA_snn_res.0.8
    idx = which(Idents(scRNA) %in% keep.clusters); # & (!grepl("TRCnormal", colnames(scRNA))))
    scRNA = subset(scRNA, cells = idx);
    scRNA = subset(scRNA, features = which(Matrix::rowSums(scRNA@assays$RNA@counts >0) > 0.05*ncol(scRNA)));
    print(ID);
    print(scRNA);
    if(0){
        x = t(as.matrix(scRNA@assays$RNA@counts))
        file = paste0(output, "-", ID, ".tsv.gz")
        write.table(sep = "\t", x, file =gzfile(file), quote = F, col.names=NA)
    }
    if(0){
        #cluster again and check marker expression;
        scRNA <- NormalizeData(scRNA)
        scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
        scRNA <- ScaleData(scRNA)
        scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
        
        scRNA <- RunUMAP(scRNA, dims = 1:10)
        p = FeaturePlot(scRNA, features = c(macrophage.markers, macroglia.vs.macrophage) ) ;
        png(paste0(output, "-", ID, ".mk.png"), width = 10, height = 10, unit = "in", res = 150)
        print(p);
        dev.off();
        cluster = cbind(scRNA[[]], scRNA@reductions$umap@cell.embeddings, scRNA@reductions$pca@cell.embeddings[,1:20]);
        saveRDS(cluster, file = paste0(output,".cluster.rds"));
    }
    #file;
    return(scRNA);
})

###feature selection
features = lapply(file.list, function(scRNA){
    scRNA <- NormalizeData(scRNA)
    scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = n.genes)
    VariableFeatures(scRNA)
})
sel.features = names(which(table(unlist(features)) > 3));

message("Merge scRNA datas...")
scRNA = merge(file.list[[1]], y = file.list[-1], add.cell.ids = sample.list, project = "GBM")
#scRNA = subset(scRNA, features = which(Matrix::rowSums(scRNA@assays$RNA@counts >0) > 0.05*ncol(scRNA)));
VariableFeatures(scRNA) = sel.features;
scRNA <- NormalizeData(scRNA)
#scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(scRNA)
scRNA <- RunUMAP(scRNA, dims = 1:10)
p = FeaturePlot(scRNA, features = c(macrophage.markers, macroglia.vs.macrophage) ) ;
png(paste0(output, ".mk.png"), width = 10, height = 10, unit = "in", res = 150)
print(p);
dev.off();
cluster = cbind(scRNA[[]], scRNA@reductions$umap@cell.embeddings, scRNA@reductions$pca@cell.embeddings[,1:20]);
saveRDS(cluster, file = paste0(output,".cluster.rds"));
save(scRNA, file = paste0(output, ".seurat.rda"));

###
x = t(as.matrix(scRNA@assays$RNA@counts[sel.features,])); rm(scRNA);
file = paste0(output, ".tsv.gz")
write.table(sep = "\t", x, file =gzfile(file), quote = F, col.names=NA)
}

#q(save = "no")

file = paste0(output, ".tsv.gz")
#run cNMF
parallelize_numworkers = paste(0:(n.core-1), collapse = " ")

    #file = file.list[[k]]
    k.output = output;
### 1 cNMF prepare
#cmd1 = paste(sep=" ", nmf.cmd, "prepare -c ", paste0(output, ".tsv.gz"), "--output-dir ", output, "--name", name, "-k", k.nmf, "--n-iter 50 --seed 123 --numgenes ", n.genes, " --total-workers", n.core)
    cmd1 = paste(sep=" ", nmf.cmd, "prepare -c ", file, "--output-dir ", k.output, "-k", k.nmf, "--n-iter 50 --seed 14 --numgenes ", n.genes, " --total-workers", n.core)
    cat(cmd1, "\n")
    system(cmd1);

### 2 cNMF factorize
#conda activate cnmf_env
#python /data/wuz6/software/cNMF/cnmf.py prepare -c 02.GBM.scRNA.macrophage.tsv.gz --output-dir macrophage_cNMF --name macrophage -k 5 6 7 8 9 10 --n-iter 50 --seed 14 --numgenes 2000 --total-workers 10#parallel python /data/wuz6/software/cNMF/cnmf.py factorize --output-dir macrophage_cNMF --name macrophage --worker-index {} ::: 0 1 2 3 4 5 6 7 8 9
#python /data/wuz6/software/cNMF/cnmf.py combine --output-dir macrophage_cNMF --name macrophage
#python /data/wuz6/software/cNMF/cnmf.py k_selection_plot --output-dir macrophage_cNMF --name macrophage

#cmd2 = paste(sep=" ", "parallel", nmf.cmd, "factorize --output-dir", output, "--name", name, "--worker-index {} :::", parallelize_numworkers)
    cmd2 = paste(sep=" ", "parallel", nmf.cmd, "factorize --output-dir", k.output, "--worker-index {} :::", parallelize_numworkers)
    cat(cmd2, "\n")
    system(cmd2);

### cNMF combine
#cmd3 = paste(sep = " ", nmf.cmd, "combine --output-dir", output, "--name", name)
    cmd3 = paste(sep = " ", nmf.cmd, "combine --output-dir", k.output)
    cat(cmd3, "\n")
    system(cmd3)

#rm ./example_data/example_cNMF/cnmf_tmp/example_cNMF.spectra.k_*.iter_*.df.npz

### Step 4 select an optimal K by considering the trade-off between stability and error
#cmd4 = paste(sep = " ", nmf.cmd, "k_selection_plot --output-dir", output, "--name", name);
    cmd4 = paste(sep = " ", nmf.cmd, "k_selection_plot --output-dir", k.output);
    cat(cmd4, "\n")
    system(cmd4)
    1;
#python /data/wuz6/software/cNMF/cnmf.py consensus --output-dir macrophage_cNMF --show-clustering --local-density-threshold 2.00 --components 10
cat("Done\n")

