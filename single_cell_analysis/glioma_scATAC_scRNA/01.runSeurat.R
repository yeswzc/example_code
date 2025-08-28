library(Seurat)
#library(harmony)
library(viridis)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()

#library(clustree)

#library(NMF)
#source("/data/wuz6/project/01.single_cell_ATAC/01.scRNA/GBM/runNMF.R")
#
output = "01.GBM"
data.home = "/data/wuz6/data/aldape_lab/scRNA/cellranger/"
patient.list = c("ABS", "CCH", "MBF", "MBT", "RCA", "ROB", "STB", "STR", "GP27", "GP28")
dim.use = 20;

### Marker genes
mk.list = "/data/wuz6/data/pipeline/celltypeMarker/GliomaSingleCellMarkers.csv";
mk = read.csv(mk.list, head =T, skip =1, stringsAsFactors = F)
all.marker = unique(mk$Gene[order(mk$Celltype)])

if(1){
    #ABS  CCH  MBF  MBT  RCA  ROB  STB  STR  TRCnormal
    #merge data from different samples
    data.list = lapply(patient.list, function(id){
            cat("Reading", id, "\n")
            sc.data = Read10X(data.dir = file.path(data.home, id, "/outs/filtered_feature_bc_matrix/"))
            scRNA <- CreateSeuratObject(counts = sc.data, min.cells = 3, min.features = 200, project = id); rm(sc.data);
            scRNA;
            scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
    #Plot features
            cat("Quantiles of mitochondrial content\n")
            #aggregate(scRNA[['percent.mt']], scRNA[['orig.ident']], median)
            print(quantile(scRNA[['percent.mt']][,1]))
            cat("Quantiles of nCount_RNA/UMI\n")
            print(quantile(scRNA[[]]$nCount_RNA));

            p1 = ggplot(data.frame(scRNA[[]]), aes(x= nCount_RNA, y = percent.mt, col = orig.ident)) +
                  geom_point(size = 0.2, alpha = 0.2) + theme_bw()+ #scale_color_manual(values = patient.col)+
                  guides(color = guide_legend(override.aes = list(size=2, alpha = 1)))
            p2 = ggplot(data.frame(scRNA[[]]), aes(x= nFeature_RNA, y = percent.mt, col = orig.ident)) +
                  geom_point(size = 0.2, alpha = 0.2) + theme_bw()+ #scale_color_manual(values = patient.col)+
                  guides(color = guide_legend(override.aes = list(size=2, alpha = 1)))
            p3 = ggplot(data.frame(scRNA[[]]), aes(x= nCount_RNA, y = nFeature_RNA, col = orig.ident)) +
                  geom_point(size = 0.2, alpha = 0.2) + theme_bw()+ #scale_color_manual(values = patient.col)+
                  guides(color = guide_legend(override.aes = list(size=2, alpha = 1)))
            png(paste0(output,"-",id ,".featureScatter.png"), width = 15, height = 5, units = "in", res = 300)
            print(cowplot::plot_grid(p1, p2, p3, nrow = 1))
            dev.off()
            
            #scRNA <- subset(scRNA, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & nCount_RNA >1800 & percent.mt < 20)
            scRNA <- subset(scRNA, subset = nFeature_RNA > 500 & percent.mt < 20)
            cat("Non sample-wise feature selection using SCT...\n");
            scRNA <- SCTransform(scRNA, verbose = T,vars.to.regress = "percent.mt", variable.features.n = 3000); #variable.features.n default is 3000
            scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA));
            cat("TSNE and cluster...");
            scRNA <- RunTSNE(scRNA, dims = 1:dim.use, verbose = FALSE) #, reduction = "harmony")
            scRNA <- RunUMAP(scRNA, dims = 1:dim.use, verbose = FALSE) #, reduction = "harmony")
            scRNA <- FindNeighbors(scRNA, dims = 1:dim.use, verbose = FALSE) #, reduction = "harmony")
            #scRNA <- FindClusters(scRNA, verbose = FALSE, resolution = seq(0.5, 1.2, 0.1));
            scRNA <- FindClusters(scRNA, verbose = FALSE, resolution = 0.6);
            # SingleR annotation 
            pred <- SingleR(test =scRNA@assays$RNA@counts, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main);
            scRNA = AddMetaData(scRNA, pred$pruned.labels, col.name= "SingleR")

            print(head(scRNA[[]]))

            res = cbind(scRNA[[]], scRNA@reductions$tsne@cell.embeddings, scRNA@reductions$umap@cell.embeddings,
                        scRNA@reductions$pca@cell.embeddings[,1:20])
            res$cluster = res$SCT_snn_res.0.6; 
            #plot PCA
            p.pca = lapply(seq(1, 20, 2), function(i){
                p = ggplot(res, aes_string(x = paste0("PC_", i), y = paste0("PC_", i+1), col = "orig.ident"))+ #seurat PCA
                    geom_point(size = 0.2) + theme_bw() + theme(aspect.ratio = 1) + #scale_color_manual(values = patient.col)+
                      guides(color = guide_legend(override.aes = list(size=2))) + theme(legend.position = "none");
                p;
            })
            png(paste0(output,"-", id, ".PCA.png"), width = 15, height = 8, unit = "in", res = 300)
            print(cowplot::plot_grid(plotlist = p.pca, ncol = 5))
            dev.off()
            ### plot cluster
            #tsne.center = aggregate(subset(res, select = c("tSNE_1", "tSNE_2")), list(res$cluster), median)
            #colnames(tsne.center)[1] = "cluster"
            umap.center = aggregate(subset(res, select = c("UMAP_1", "UMAP_2")), list(res$cluster), median)
            colnames(umap.center)[1] = "cluster"

            p1 = ggplot(res, aes(x=UMAP_1, y = UMAP_2, col = cluster)) +
                  geom_point(size = 0.05)+theme_bw()+
                  guides(color = guide_legend(override.aes = list(size=1), nrow =5))+
                  geom_label(data = umap.center, aes(label = cluster), fill = "grey",
                             alpha = 0.3, colour = "black", fontface = "bold", label.size = NA)+
                  theme(aspect.ratio = 1, legend.position = "none", panel.grid = element_blank(), 
                             axis.text = element_blank(), axis.ticks = element_blank())+
                  labs(x = "", y = "", color = "")
            p2 = ggplot(res, aes(x=UMAP_1, y = UMAP_2, col = SingleR)) +
                  geom_point(size = 0.05)+theme_bw()+
                  #scale_color_viridis(discrete = TRUE, option = "D")+
                  guides(color = guide_legend(override.aes = list(size=1), nrow =5))+
                  geom_label(data = umap.center, aes(label = cluster), fill = "grey",
                             alpha = 0.3, colour = "black", fontface = "bold", label.size = NA)+
                  theme(aspect.ratio = 1, legend.position = "right", panel.grid = element_blank(), 
                             axis.text = element_blank(), axis.ticks = element_blank())+
                  labs(x = "", y = "", color = "");
            p2.legend = legend <- cowplot::get_legend(p2 + theme(legend.box.margin = margin(0, 0, 0, 12)))
            p2 = p2 + theme(legend.position = "none")
            png(paste0(output,"-", id, ".umap.png"), width = 15, height = 5, unit = "in", res = 150)
            print(cowplot::plot_grid(p1, p2, p2.legend, nrow = 1))
            dev.off()
            #plot marker genes
            cat("Plot marker genes...\n")
            scRNA = NormalizeData(scRNA)
            png( paste0(output, "-", id, ".markerGenes.png"), width = 26, height = 12, unit = "in", res = 100)
            p = FeaturePlot(scRNA, pt.size = 0.02, reduction = "umap", features = all.marker, coord.fixed = T, ncol = 11, slot = "data")
            print(p);
            dev.off()
            ### Find marker genes
            sc.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod", verbose = F)  %>% group_by(cluster)
            write.table(sc.markers, file = paste0(output,"-", id, ".Markers.tsv"), sep = "\t", quote=F, row.names = F)
            #save(res, file = paste0(output, "-",id, ".cluster.rda"))
            save(scRNA, file = paste0(output,"-", id, ".SeuratObjCluster.rda"));
            id;
    })
}

sessionInfo()

