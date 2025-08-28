#R/4.3
library(Seurat)
library(parallel)
library("foreach")
library("doParallel")
library("doMC")
library(NMF)
registerDoMC(cores = 6) #required


args = commandArgs(trailingOnly=T)
i_sample = as.numeric(args[1])
cat("Working on ", i_sample, " sample\n")
###
nFeature_Spatial.cut.off <- 800;
nCount_Spatial.cut.off <- 1.5*1e3; #ZC initial
#nFeature_Spatial.cut.off <- 0;
#nCount_Spatial.cut.off <- 1.0*1e3;

rank_lb = 2
rank_ub = 11
out.path = "01.NMF/"
###
data.home = "/data/wuz6/data/aldape_lab/spatial10xg/spaceranger_v201_he/"

#all.sample.names <- read.table("00.pHGG_GBM_IDH.sample.list", head=F, stringsAsFactors = F)[,1]
#all.sample.names <- c("AP67PHE","U481PHE","Y963_B11PHE","U529","X744","Y997","AD01","V203","X210","W908_1O","Y168","AE36","AL52","Y216","Z256","Z930","W396")
all.sample.names <- c("AP67PHE","Y963_B11PHE","U529","X744","Y997","AD01","V203","X210","W908_1O","Y168","AE36","AL52","Y216","Z256","Z930","W396")
all.sample.names <- c('AE30', 'BK67', 'Y758', 'AR56_A3', 'Y597', 'CP86')
all.sample.names <- c('Y597')
sample.name <- all.sample.names[i_sample]
cat("Working on ", sample.name, "\n")

if(file.exists(paste0(out.path, sample.name, ".nmf.rds"))){return(NULL)}

###
#lapply(all.sample.names, function(sample.name){
    cat("Loading ", sample.name, "\n")
    message("Loading ", sample.name)
    #spata.coord.file = paste0(data.home, sample.name, "/outs/spatial/tissue_positions.csv")
    #spata.coord = read.csv(spata.coord.file, head = T, check.names = F, row.names = 1)
    #colnames(spata.coord)[c(4:5)] <- c("spatialX", "spatialY")
    sp.data = Load10X_Spatial(paste0(data.home, sample.name, "/outs/"), filename = "filtered_feature_bc_matrix.h5", assay = "Spatial")
    sp.data <- PercentageFeatureSet(sp.data, pattern = "^MT-", col.name = "percent.mt")
    sp.data <- subset(sp.data, subset = nFeature_Spatial > nFeature_Spatial.cut.off & nCount_Spatial > nCount_Spatial.cut.off & percent.mt < 20);
    nonMT.genes = grep(pattern = "^MT", rownames(sp.data), value = T, invert = T)
    sp.data = subset(sp.data, features = nonMT.genes)
    m = data.matrix(sp.data[["Spatial"]]$counts);rm(sp.data)
    ###process data based method provided by https://github.com/tiroshlab/Spatial_Glioma/blob/main/Module1_MP.R
    if(min(colSums(m)) == 0){m <- m[, colSums(m) != 0]}
    scaling_factor <- 1000000/colSums(m)
    m_CPM <- sweep(m, MARGIN = 2, STATS = scaling_factor, FUN = "*")
    m_loged <- log2(1 + (m_CPM/10))
    # removing genes with zero variance across all cells
    var_filter <- apply(m_loged, 1, var)
    m_proc <- m_loged[var_filter != 0, ]
    # filtering out lowly expressed genes
    exp_genes <- rownames(m_proc)[(rowMeans(m_proc) > 0.4)]
    m_proc <- m_proc[exp_genes, ]
    # centering data gene-wise
    count_mat_cent<- m_proc - rowMeans(m_proc)
    #nmf preprocessing
    count_mat_nmf<- count_mat_cent
    count_mat_nmf[count_mat_nmf<0] <- 0 # negative values should be set to 0 in initial matrix
    # output to a list of gene expression profiles (GEP)
    #per_sample_mat[[i]] <- count_mat_nmf
    #names(per_sample_mat)[i] <- sample_ls[[i]]
    rm(m, m_CPM, m_loged, m_proc, count_mat_cent)
    #print(count_mat_nmf[1:4,1:4])
    ###
    #res <- NMF::nmf(x = count_mat_nmf, rank = rank_lb:rank_ub, nrun = 5, method = "snmf/r", .opt = list(debug=F, parallel=F, shared.memory=F, verbose=T))
    res <- NMF::nmf(x = count_mat_nmf, rank = rank_lb:rank_ub, nrun = 5, method = "snmf/r", .opt = 'vP6', .pbackend=NULL)
    saveRDS(res, file = paste0(out.path, sample.name, ".nmf.rds"))
    1

#}, mc.cores=5)

##

