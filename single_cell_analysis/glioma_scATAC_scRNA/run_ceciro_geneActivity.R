library(SnapATAC)
library(cicero)

input.data = "GBM.inhouse.PeakCluster.xsp_final.rda"
output = "test"
load(input.data)


cat("Remove unwantted chromosomes.\n");
#chr.exclude = seqlevels(x.sp@peak)[grep("random|chrEBV|chrX|chrY|chrM|chrUn", seqlevels(x.sp@peak))];
#idy = grep(paste(chr.exclude, collapse="|"), x.sp@peak);
idy = grep("random|chrEBV|chrX|chrY|chrM|chrUn", x.sp@peak);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="pmat"]; cat("removed..\n")}
cat(paste0(unique(seqlevels(x.sp@feature))[1:30], collapse = ","),"... ", length(idy)," bins removed.\n");
x.sp;

#x.sp@barcode
#x.sp@pmat
#x.sp@peak
#identical(row.names(x.sp@metaData), x.sp@barcode)
#cellinfo <- x.sp@metaData
#names(cellinfo)[1] <- "cells"
cellinfo = data.frame( paste0(x.sp@sample, "#", x.sp@barcode))
rownames(cellinfo) = cellinfo[,1]
names(cellinfo)[1] <- "cells"



peakinfo = data.frame(x.sp@peak)[,1:3]
colnames(peakinfo) = c("chr", "bp1", "bp2")
peakinfo[,1] = gsub("[b']", "", peakinfo[,1])
peakinfo$bp1 = peakinfo$bp1 - 1
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name


indata = t(x.sp@pmat)
row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)
colnames(indata) = paste0(x.sp@sample, "#", x.sp@barcode)

# make CDS
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                                              phenoData = pd,
                                              featureData = fd,
                                              expressionFamily=VGAM::binomialff(),
                                              lowerDetectionLimit=0))

input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

## Constructing cis-regulatory networks
### Running Cicero
#### Create a Cicero CDS #required
set.seed(2020)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

# *** if you are using Monocle 3 alpha, you need to run the following line as well!
# NOTE: Cicero does not yet support the Monocle 3 beta release (monocle3 package). We hope
# to update soon!
#input_cds <- preprocessCDS(input_cds, norm_method = "none")
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim= 50,
                             reduction_method = 'tSNE', norm_method = "none") #a few minutes

#Next, we access the tSNE coordinates from the input CDS object where they 
#are stored by Monocle and run make_cicero_cds:
tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)
#?or?
#cicero_cds <- make_cicero_cds(input_cds, k = 50, reduced_coordinates = x.sp@tsne)
#### Run Cicero #Required
#data("human.hg38.genome")
hg38.size = read.table("/data/wuz6/reference/refdata-cellranger-atac-GRCh38-1.2.0.genome.fa.sizes", stringsAsFactors=F)
#sample_genome <- subset(hg38.size, V1 == "chr18")
conns <- run_cicero(cicero_cds, hg38.size) # Takes a few minutes to run
#head(conns)


#Cicero gene activity scores
#### Add a column for the fData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
#devtools::install_github("cran/refGenome")

gtf.file = "/fdb/cellranger-atac/refdata-cellranger-atac-GRCh38-1.2.0/genes/genes.gtf"
gene_annotation = data.frame(rtracklayer::import(gtf.file))
gene_annotation = gene_annotation[gene_annotation$type == "exon", ];
gene_annotation = subset(gene_annotation, 
           select = c("seqnames", "start", "end", "strand", "gene_type", "gene_id", "transcript_id", "gene_name"))
colnames(gene_annotation) = c("chromosome", "start", "end", "strand", "feature", "gene", "transcript", "symbol")
save(gene_annotation, file = "hg38.cellranger1.2.0.gene_annotation.rda")
#load("hg38.cellranger1.2.0.gene_annotation.rda")

pos <- subset(gene_annotation, strand == "+")
pos <- pos[order(pos$start),] 
pos <- pos[!duplicated(pos$transcript),] # remove all but the first exons per transcript
pos$end <- pos$start + 1 # make a 1 base pair marker of the TSS

neg <- subset(gene_annotation, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
neg <- neg[!duplicated(neg$transcript),] # remove all but the first exons per transcript
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,c(1:3, 8)]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"
gene_annotation_sub = gene_annotation_sub[gene_annotation_sub[,1] !="chrX" & gene_annotation_sub[,1] != "chrY" & gene_annotation_sub[,1] != "chrM",]
input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

#head(fData(input_cds))

#### Generate gene activity scores ####
# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)

# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# normalize #dgCMatrix
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
save(input_cds, file = paste0(output, ".input_cds.rda"))
saveRDS(cicero_gene_activities, file = paste0(output, ".ceciro_gene_activity.normalized.rds"))
saveRDS(unnorm_ga, file = paste0(output, ".ceciro_gene_activity.raw.rds"))
saveRDS(conns, file = paste0(output, ".ceciro_connections.rds"))

# if you had two datasets to normalize, you would pass both:
# num_genes should then include all cells from both sets
#unnorm_ga2 <- unnorm_ga
#cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), num_genes)




