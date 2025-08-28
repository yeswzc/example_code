library(Seurat)
library(scDblFinder)
options(future.globals.maxSize = 59*1000 * 1024^2)

##############################################################################
##############################################################################
all.samples <- list.files("data/")
all.samples
#all.samples <- grep("png|VDJ", all.samples, invert = T, value = T)
all.samples <- grep("vdj", all.samples, invert = T, value = T)
all.samples

###merge single cell RNA
sc.data.list = lapply(all.samples, function(id){
  cat("Reading", id, "\n")
  sc = Read10X(data.dir = file.path("data", id, "filtered_feature_bc_matrix/"))
  sc <- CreateSeuratObject(counts = sc, min.cells = 3, min.features = 500, project = id); 
  sc <- RenameCells(sc, paste0(id, "_", colnames(sc)))
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  #sc <- add_clonotype(paste0(id, "-VDJ/"), sc, "T")
  sc <- subset(sc, subset = nFeature_RNA > 500 & percent.mt < 20)
  sce <- scDblFinder(sc[['RNA']]$counts)
      #sce <- findDoubletClusters(sce)
  sce.dbl <- data.frame(sce@colData); rm(sce)
  sc <- AddMetaData(sc, sce.dbl)
  #sc <- subset(sc, scDblFinder.class == "singlet")
  sc$orig.ident <- id
  return(sc);
})
##
save(sc.data.list, file  = "sc.rna.filterDoublet.rda")
