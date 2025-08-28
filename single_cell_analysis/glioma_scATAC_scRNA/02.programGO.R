#meta program GO enrichment
args = commandArgs(trailingOnly = T)

if(length(args) != 1) stop("Usage: [dir/consensus.txt]")

in.file = args[1]
base.dir = dirname(in.file)

library(clusterProfiler);
library(dplyr)
library(ggplot2)

spectra = read.table(in.file, head =T, sep = "\t", check.names= F)


db = msigdbr::msigdbr(species = "human", category = "C5")
db = db %>% dplyr::distinct(gene_symbol, gs_name) %>% as.data.frame()


k.gene  = 30;
go = lapply(1:nrow(spectra),function(k){
    G = colnames(spectra)[order(spectra[k,], decreasing = T, na.last = T)[1:k.gene]];
    res = data.frame(enricher(gene = G, TERM2GENE = db), check.names = F);
    if(nrow(res) < 1){
        message("K = ", k, "has no enrichment...")
        return(1);
    }
    write.table(res, file = paste0(base.dir, "/GO.", k, ".xls"));
    if(nrow(res) > 20) res = res[1:20,]
    res$ID = gsub("GO\\w\\w_", "", res$ID)
    p= ggplot(res, aes(x = reorder(ID, -p.adjust), y = -log10(p.adjust))) +
      geom_bar(stat = "identity") + theme_bw() +
      scale_x_discrete(label = function(x) stringr::str_trunc(x, 45))+
      coord_flip() +
      labs(x= "", y = "")

    png(file = paste0(base.dir, "/GO.", k, ".png"), width = 7, height = 5, units = "in", res = 150)
    print(p)
    dev.off()

    res;
});




