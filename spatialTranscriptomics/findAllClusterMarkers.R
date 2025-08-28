#R/4.3
library(Seurat)
library(Banksy)
library(parallel)
library(SeuratWrappers)
library(clustree)
library(dplyr)
Sys.setenv(RETICULATE_PYTHON = "/usr/local/Anaconda/envs/py3.10/bin/python")

###
nFeature_Spatial.cut.off <- 800;
nCount_Spatial.cut.off <- 1.5*1e3;

all.sample.names <- read.table("00.pHGG_GBM_IDH.sample.list", head=F, stringsAsFactors = F)[,1]
#remove U481PHE
all.sample.names <- c("AP67PHE","Y963_B11PHE","U529","X744","Y997","AD01","V203","X210","W908_1O","Y168","AE36","AL52","Y216","Z256","Z930","W396")
all.sample.names <- "W908_1O"
all.sample.names <- c('AE30', 'BK67', 'Y758', 'AR56_A3', 'Y597', 'CP86')
all.sample.names <- c("Y597")
data.home = "/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/"

all.folder.names = all.sample.names
n.cores = 1
###
###
###
###
if(length(all.sample.names) == 0){quit(save = "no")}

message("left samples: ", all.sample.names)
mclapply(1:length(all.sample.names), function(k){
        ###0 load data
    #spata.coord.file = paste0(data.home, sample.name, "/outs/spatial/tissue_positions.csv")
    sample.folder = all.folder.names[k]
    sample.name = all.sample.names[k]
    sample.name = gsub("#", "", sample.name)

    sp.data = Load10X_Spatial(paste0(data.home, sample.folder, "/outs/"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
    sp.data <- PercentageFeatureSet(sp.data, pattern = "^MT-", col.name = "percent.mt")
    sp.data <- subset(sp.data, subset = nFeature_Spatial > nFeature_Spatial.cut.off & nCount_Spatial > nCount_Spatial.cut.off & percent.mt < 20);
    #spatial.coord = GetTissueCoordinates(sp.data)
    #colnames(spata.coord) <- c("row", "col")

    sp.data = NormalizeData(sp.data) #need normalized for cell cycle score
    ###1.5
    message("Banksy")
    banksy.cluster.res <- read.csv( paste0("01.",sample.name, ".Banksycluster.csv"), row.names =1, check.names =F)
    use.clusters = grep("BANKSY_snn_res", colnames(banksy.cluster.res), value = T)
    #banksy.cluster.res <- read.csv( paste0("01.",sample.name, ".Seuratcluster.csv"), row.names =1, check.names =F)
    #use.clusters = grep("SCT_snn_res", colnames(banksy.cluster.res), value = T)
    sp.data = AddMetaData(sp.data, banksy.cluster.res)
    all.mk = lapply(use.clusters, function(kk){
        Idents(sp.data) = as.factor(sp.data@meta.data[[kk]])
        mk <- FindAllMarkers(sp.data, only.pos = TRUE, assay = "Spatial")
        mk$res = kk
        #mk$sample = sample.name
        mk
    })
    all.mk = do.call(rbind, all.mk)
    write.csv(all.mk, file = paste0("01.",sample.name, ".AllBanksyMarkers.csv"))
    return(1);;

    ###
    rm(banksy.cluster.res, use.clusters)
    
    banksy.cluster.res <- read.csv( paste0("01.",sample.name, ".Seuratcluster.csv"), row.names =1, check.names =F)
    use.clusters = grep("SCT_snn_res", colnames(banksy.cluster.res), value = T)
    sp.data = AddMetaData(sp.data, banksy.cluster.res)
    all.mk = lapply(use.clusters, function(kk){
        Idents(sp.data) = as.factor(sp.data@meta.data[[kk]])
        mk <- FindAllMarkers(sp.data, only.pos = TRUE, assay = "Spatial")
        mk$res = kk
        #mk$sample = sample.name
        mk
    })
    all.mk = do.call(rbind, all.mk)
    write.csv(all.mk, file = paste0("01.",sample.name, ".AllSeuratMarkers.csv"))
    1;

}, mc.cores = n.cores)

##

