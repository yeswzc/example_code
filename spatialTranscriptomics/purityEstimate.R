#R/4.3
library(Seurat)
library(parallel)
library(tidyestimate)


nCount_Spatial.cut.off <- 1000
nFeature_Spatial.cut.off <- 0
###
data.home = "/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/"
#all.sample.names <- read.table("00.sample.list.txt", head=F, stringsAsFactors = F)[,1]
#all.sample.names <- read.table("00.sample.list.txt.GBM", head=F, stringsAsFactors = F)[,1]
#all.sample.names <- read.table("00.pHGG_GBM_IDH.sample.list", head=F, stringsAsFactors = F)[,1]
#all.sample.names <- c("AP67PHE","U481PHE","Y963_B11PHE","U529","X744","Y997","AD01","V203","X210","W908_1O","Y168","AE36","AL52","Y216","Z256","Z930","W396")
all.sample.names <- read.table("/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/Data_locations_v20250716.txt", head =T, sep = "\t")[,2]
n.cores = 2
###

message("left samples: ", paste0(all.sample.names, collapse = ";"))

all.scores <- mclapply(1:length(all.sample.names), function(k){
        ###0 load data
    #data.home <- info[k,1]
    #sample.name <- info$Sample[k]
    sample.name <- all.sample.names[k]
    message(sample.name)
    
    sp.data = Load10X_Spatial(paste0(data.home,sample.name, "/outs/"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
    sp.data <- PercentageFeatureSet(sp.data, pattern = "^MT-", col.name = "percent.mt")
    sp.data <- subset(sp.data, subset = nFeature_Spatial > nFeature_Spatial.cut.off & nCount_Spatial > nCount_Spatial.cut.off & percent.mt < 20);
    sp.data <- NormalizeData(sp.data)
    sp.data <- data.matrix(sp.data[['Spatial']]$data)
    #print(sp.data[1:4,1:4])
    #message("scoring")
    scores <- sp.data |> filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> estimate_score(is_affymetrix = T) #is_affymetrix=T has purity, otherwise does not
    print(head(scores))
    scores$sampleName <- sample.name
    print(head(scores))

    return(scores)

}, mc.cores = n.cores)
names(all.scores) <- all.sample.names
saveRDS(all.scores, file = "01.inhouse.ESTIMATE.scores.rds")
##

