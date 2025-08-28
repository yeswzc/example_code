args = commandArgs(trailingOnly = T)
if(length(args) != 2){stop("Usage: [seurat.obj.rda] [output.png]\n\n"); }

library(Seurat)

input = args[1];
output = args[2] 

load(input);
DefaultAssay(scRNA) = "RNA";
scRNA = NormalizeData(scRNA);
mk.list = "/data/wuz6/data/pipeline/celltypeMarker/GliomaSingleCellMarkers.csv";
mk = read.csv(mk.list, head =T, skip =1, stringsAsFactors = F)
#Glioma single cell marker,,
#Gene,Celltype,Paper
#HES1,quiescence,"Anoop et al. Science, 2014"
#TSC22D1,quiescence,"Anoop et al. Science, 2014"
all.marker = unique(mk$Gene[order(mk$Celltype)])

png(output, width = 26, height = 12, unit = "in", res = 100)
FeaturePlot(scRNA, pt.size = 0.02, reduction = "umap", features = all.marker, coord.fixed = T, ncol = 11, slot = "data")
dev.off()
