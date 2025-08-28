#R/4.3
library(Seurat)
library(Banksy)
library(parallel)
library(SeuratWrappers)
library(clustree)
library(dplyr)


###
banksy.cluster <- function(seurat.obj = NULL, n.pc = 20){
    ###with spatial 
    seurat.obj <- RunBanksy(seurat.obj, lambda = 0.2, verbose=TRUE, dimx = "X", dimy = "Y",
        assay = 'SCT', slot = 'data', features = 'variable', k_geom = 15)
        #assay = 'Spatial', slot = 'data', features = 'variable', k_geom = 15)
    seurat.obj <- RunPCA(seurat.obj, assay = 'BANKSY', features = rownames(seurat.obj), npcs = n.pc)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:n.pc)
# Clustering
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:n.pc)
    seurat.obj <- FindClusters(seurat.obj, cluster_id = "BANKSY", resolution = seq(0.6, 1.2, 0.05), algorithm = 4)
    res = seurat.obj@meta.data
    res = cbind(seurat.obj@reductions$umap@cell.embeddings, res[,grep("BANKSY", colnames(res))])
    res
}
###
nFeature_Spatial.cut.off <- 800; #0
nCount_Spatial.cut.off <- 1.5*1e3; #1000

###
data.home = "/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/"
#all.sample.names <- read.table("00.sample.list.txt", head=F, stringsAsFactors = F)[,1]
#all.sample.names <- read.table("00.sample.list.txt.GBM", head=F, stringsAsFactors = F)[,1]
all.sample.names <- read.table("00.pHGG_GBM_IDH.sample.list", head=F, stringsAsFactors = F)[,1]
#remove U481PHE
all.sample.names <- c("AP67PHE","Y963_B11PHE","U529","X744","Y997","AD01","V203","X210","W908_1O","Y168","AE36","AL52","Y216","Z256","Z930","W396")
all.sample.names <- "W908_1O"
all.sample.names <- c('AE30', 'BK67', 'Y758', 'AR56_A3', 'Y597', 'CP86')
all.sample.names <- c("Y597")
n.cores = 1
###
all.sample.names <- lapply(all.sample.names, function(sample.name){
    if( file.exists(paste0("01.",sample.name, ".Banksycluster.csv")) ){message("Skip ", sample.name); return(NULL)}
    sample.name;
})
all.sample.names = unlist(all.sample.names)


message("left samples: ", all.sample.names)
mclapply(1:length(all.sample.names), function(k){
        ###0 load data
    #data.home <- info[k,1]
    #sample.name <- info$Sample[k]
    sample.name <- all.sample.names[k]
    spata.coord.file = paste0(data.home, sample.name,"/outs/spatial/tissue_positions.csv")
    spata.coord = read.csv(spata.coord.file, head = T, check.names = F, row.names = 1)
    colnames(spata.coord)[c(4:5)] <- c("spatialX", "spatialY")
    n.pc <- 20 #for clustering and umap

    sp.data = Load10X_Spatial(paste0(data.home,sample.name, "/outs/"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
    sp.data <- PercentageFeatureSet(sp.data, pattern = "^MT-", col.name = "percent.mt")
    sp.data <- subset(sp.data, subset = nFeature_Spatial > nFeature_Spatial.cut.off & nCount_Spatial > nCount_Spatial.cut.off & percent.mt < 20);

    if(1){
        sp.data = NormalizeData(sp.data) #need normalized for cell cycle score
        sp.data = CellCycleScoring(sp.data, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
        ###remove MT genes for clustering
        keep.nonMT.genes = grep(pattern = "^MT", rownames(sp.data), value = T, invert = T)
        sp.data = subset(sp.data, features = keep.nonMT.genes)
        ####
        ###1.1 seurat preprocess and cluster
        sp.data <- SCTransform(sp.data, assay = "Spatial", verbose = FALSE, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), variable.features.n = 7000)
        spata.coord.match = spata.coord[rownames(sp.data@meta.data),]
        sp.data@meta.data[,c("X", "Y")] = spata.coord.match[,c("spatialX", "spatialY")]
    }else{
        sp.data = NormalizeData(sp.data, scale.factor = 100, normalization.method = 'RC')
    }
      
    ###1.5
    message("Banksy")
    banksy.cluster.res <- banksy.cluster(sp.data, n.pc = n.pc)
    prefix = "BANKSY_snn_res."
    out = clustree(banksy.cluster.res, prefix = prefix)
    out.plot <- out$data %>% select(!!sym(prefix), sc3_stability) %>% 
            group_by(!!sym(prefix)) %>%   summarise(mean(sc3_stability), sd(sc3_stability)) %>% 
            as.data.frame() 
    best.stable = as.character(out.plot[,1][which.max(out.plot$`mean(sc3_stability)`)])
    banksy.cluster.res$best = best.stable
    banksy.cluster.res$banksy.cluster = banksy.cluster.res[, paste0(prefix, best.stable)]
    #pdf(paste0("01.",sample.name, ".RCnormBanksyClustree.pdf"), width = 10, height = 10)

    pdf(paste0("01.",sample.name, ".BanksyClustree.pdf"), width = 10, height = 10)
    print(out)
    dev.off()
    ###
    cluster.res = banksy.cluster.res
    #write.csv(cluster.res, file = paste0("01.",sample.name, ".RCnormBanksycluster.csv"), quote = F)
    write.csv(cluster.res, file = paste0("01.",sample.name, ".Banksycluster.csv"), quote = F)

}, mc.cores = n.cores)

##

