library(spacexr)
library(Matrix)
library(data.table)
library(Seurat)

#glioma_snRNA_pair_with_label_transfer_metadata.RDS x@meta.data$nFeature_RNA
if(0){
    cnt <- readMM("../cell2location/GBmap_rawX.mtx")
    #data = readRDS("Abdelfattah.NatureComm.2022.Glioma.RDS")
    #data = subset(data, subset = Assignment != "Other")
    #data = subset(data, nFeature_RNA > 500)
    nUMI <- rowSums(cnt)
    #nUMI <- data@meta.data$nCount_RNA; names(nUMI) <- rownames(data@meta.data)
    #cell_types <- as.character(data@meta.data$Assignment); names(cell_types) <- rownames(data@meta.data)
    cell_types <- read.csv("../cell2location/GBmap.obs.csv", head =T, stringsAsFactors=F)
    cell_types <- cell_types$annotation_level_4
    cell_types <- gsub("\\/", ".", cell_types)
    #need a minimum of 25 cells for each cell type in the reference
    exclude.celltype <- names(which(table(cell_types) < 25))
    exclude.idx <- which(cell_types %in% exclude.celltype)

    cell_types <- factor(cell_types)

    cell_names <- read.csv("../cell2location/GBmap_rawX_cell_name.txt", head =F)[,1]
    gene_names <- read.csv("../cell2location/GBmap_rawX_gene_name.txt", head =F)[,1] 
    colnames(cnt) <- gene_names
    rownames(cnt) <- cell_names
    names(nUMI) <- cell_names
    names(cell_types) <- cell_names

    #nUMI <- nUMI[-exclude.idx]
    #cell_types <- factor(cell_types[-exclude.idx])
    #cnt <- t(cnt[-exclude.idx,])
    cnt <- t(cnt)
    print(cnt[1:4,1:4])

    #cell_types[cell_types == "central nervous system macrophage"] = "macrophage"
    #cell_types[cell_types == "ependymal cell"] = "ependymal"
    #cell_types[cell_types == "endothelial cell"] = "endothelial"
    #cell_types[cell_types == "vascular associated smooth muscle cell"] = "VSMC"
    #cell_types[cell_types == "Bergmann glial cell"] =   "BG"
    #cell_types[cell_types == "choroid plexus epithelial cell"] = "ChP"
    #cell_types[cell_types == "oligodendrocyte precursor cell"] = "OPC"
    #cell_types = as.factor(cell_types)
    ####
    # names(cell_types) matches colnames(counts) and names(nUMI)
    reference <- Reference(cnt, cell_types, nUMI)
    ## Examine reference object (optional)
    print(dim(reference@counts)) #observe Digital Gene Expression matrix
    table(reference@cell_types) #number of occurences for each cell type

    saveRDS(reference, file = "GBMapSCRef.rds")
}else{
    
  reference = readRDS("GBMapSCRef.rds")
}
####
make_puck <- function(sample.name, data.home){
  # set the different file paths of the filtered matrix
  # load the matrix.mtx.gz
  counts <- readMM(file = paste0(data.home, "/outs/filtered_feature_bc_matrix/matrix.mtx.gz"))
  # load the feature.tsv.gz
  feature.names = read.delim(paste0(data.home, "/outs/filtered_feature_bc_matrix/features.tsv.gz"), header=F, stringsAsFactors=F)
  # load the barcodes.tsv.gz
  barcode.names = read.delim(paste0(data.home, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"), header=F, stringsAsFactors=F)
  # set the matrix column and row names                         
  colnames(counts) = barcode.names$V1

  rownames(counts) = feature.names$V1 #V2 is gene sumbol, V1 is ENS...
  dup.gene.idx = which(duplicated(feature.names$V1))
  if(length(dup.gene.idx) >0) counts = counts[-dup.gene.idx,]

  position_file1 <- paste0(data.home, "/outs/spatial/tissue_positions.csv")
  position_file2 <- paste0(data.home, "/outs/spatial/tissue_positions_list.csv")
  position_file3 <- paste0(data.home, "/outs/spatial/tissue_positions_list.csv.gz")
  if(file.exists(position_file1)){
      coords = read.csv(position_file1, head =T, stringsAsFactors = F, row.names = 1)
  }else if(file.exists(position_file2)){
      coords = read.csv(position_file2, head =F, stringsAsFactors = F, row.names = 1)
  }else{
      coords = read.csv(position_file3, head =F, stringsAsFactors = F, row.names = 1)
  }
  coords = coords[,4:5]; colnames(coords) = c("xcoord", "ycoord")
  ####
  nUMI <- colSums(counts)
  puck <- SpatialRNA(coords, counts, nUMI)
  puck
}

norm_cell_frac <- function(object){
#as_AssayObject <- function(object) {
  if (is(object, "RCTD")) {
    if (!requireNamespace("spacexr", quietly = TRUE)) {
      stop("Install spacexr.")
    }
  } else {
    stop("Only RCTD objects supported")
  }
  r <- object@results
  if (length(r) > 1) {
    if (!is.null(r[[1]]$sub_weights)) {
      sw <- data.table::rbindlist(lapply(seq_along(r), function(i)
        data.table(
          barcode = colnames(object@spatialRNA@counts)[i],
          cell_type = names(r[[i]]$sub_weights),
          weight = r[[i]]$sub_weights
        )), fill = TRUE)
      sw$cell_type[is.na(sw$cell_type)] <- "unassigned"
      swd <- data.table::dcast(sw, barcode ~ cell_type, value.var = "weight", fill = 0)
      print(head(swd))
      swm <- as.matrix(swd[, -1])
      rownames(swm) =  swd$barcode
      message("N row swm = ", nrow(swm))
      swm <- spacexr::normalize_weights(swm)
      rownames(swm) =  swd$barcode
      message("N row swm = ", nrow(swm), " after normalize")
      #swm <- t(spacexr::normalize_weights(swm))
      #swm <- rbind(swm, max = apply(swm[!rownames(swm) %in% "unassigned", ], 2, max))
      #swm <- as(swm, "sparseMatrix")
      #return(CreateAssayObject(data = swm))
      return(swm)
    }
  } else if (length(r) == 1) {
    m <- spacexr::normalize_weights(as.matrix(r$weights))
    #m <- t(spacexr::normalize_weights(as.matrix(r$weights)))
    #m <- rbind(m, max = apply(m, 2, max))
    #return(CreateAssayObject(data = m))
    return(m)
  }
}



###

info <- read.csv("/data/wuz6/project/25.pedHGG_spatial/01.preproess/01.cluster2025HE/00.data.folder.csv", head =T, stringsAsFactors = F)
#sample.list = read.table("../../01.preproess/00.sample.list.txt", head =F, stringsAsFactors = F)[,1]

run = lapply(1:nrow(info), function(k){
    sample.name <- info$Sample[k]
    data.home <- info$folder[k]

    if(file.exists( paste0(sample.name, ".GBMap.myRCTD.rds"))){
        message( paste0(sample.name, ".GBMap.myRCTD.rds"), " existed. Next...")
        return(1)
    }
    message(sample.name, "\n")
  #sample.name = "Y963_B8ma"
    puck = make_puck(sample.name, data.home)

    myRCTD <- create.RCTD(puck, reference, max_cores = 10, CELL_MIN_INSTANCE = 20)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
    saveRDS(myRCTD, file = paste0(sample.name, ".GBMap.myRCTD.rds"))
    message("myRCTD finished...")
    ###
    #results <- myRCTD@results
    # normalize the cell type proportions to sum to 1.
    #norm_weights = normalize_weights(results$weights) 
    norm_weights = norm_cell_frac(myRCTD)
    #cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
    write.csv(norm_weights, file = paste0(sample.name, ".GBMap.norm_weights.csv"), quote = F)
    1;
})



