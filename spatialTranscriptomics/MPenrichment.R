signature.path = "/Users/wuz6/Documents/Project/data/signatures/"

all.signature.list = readRDS("/Users/wuz6/Documents/Project/25.pHGG_spatial/01.MP/all.signature.list.rds")

#background.genes = readRDS("background.genes.rds")

###
GDF= MP
GDF <- mp1

#GDF = new_metaprograms.df
#GDF = read.csv("01.MP_final.csv.csv", head =T)
GDF[1:4,1:4]

#GDF = hc.splitK1.k5.list.program.df
#GDF = data.frame(readxl::read_xlsx("/Users/wuz6/Documents/Project/data/signatures/C.Greenwald2024Cell.xlsx"))
###
#xx = nmf_programs_sig[, !colnames(nmf_programs_sig) %in% mp.pp[,2]]
#dim(xx)
#GDF = nmf_programs_sig

res = lapply(all.signature.list, function(LL){
    #LL = all.signature.list[[1]]
    #phyper()
    LL[,2] = toupper(LL[,2])
    LL.source = LL[1,3]
    message(LL.source)
    tmp.res = lapply(unique(LL$GS), function(gs){
      #gs = unique(LL$GS)[1]
      gs.genes = LL[LL$GS == gs, 2]
      n.gs.in.background = sum(gs.genes %in% background.genes)
      n.not.in.background.genes = length(background.genes) - n.gs.in.background
      ##
      sig.test = apply(GDF, 2, function(x){
        #x = GDF[,1]
        x = x[!is.na(x)]
        n.overlap = sum(x %in% gs.genes)
        phyper(n.overlap -1 , n.gs.in.background,  n.not.in.background.genes, length(x), lower.tail = F)
      })
      y = cbind(gs, sig.test, names(sig.test), LL.source)
      colnames(y) = c("GeneSet", "P.value", "MP", "GS.source")
      rownames(y) = NULL
      y
    })
    tmp.res = data.frame(do.call(rbind, tmp.res))
    #tmp.res[,2] <- as.numeric(tmp.res[,2])
    #tmp.res <- tmp.res[tmp.res[,2]< 0.01,]
    tmp.res
})


##adjust p.value
res = data.frame(do.call(rbind, res))
head(res)
res$P.value = as.numeric(as.character(res$P.value))

res$p.adjust = NA
for(i in unique(res$MP)){
  idx = which(res$MP == i)
  res$p.adjust[idx] = p.adjust(res$P.value[idx], method = "BH")
}
sum(is.na(res$p.adjust))

if(0){
   res = res[res$p.adjust<0.05,]
   res = res[order(res$MP),]
   write.csv(res, file = "pHGGGBM_all_BanksyMulti_NMFrobust_programs.enrich.csv", row.names = F)
   
   table(res$GeneSet)
   res1 = lapply( unique(res$MP), function(k){
      x = res[res$MP ==k,]
      x = x[order(x$p.adjust, decreasing = F)[1],]
      x
   })
   res1 = do.call(rbind, res1)
   table(res1$GeneSet)
    res2 = unique(res1[grepl("NPC|OPC", res1$GeneSet),]$GeneSet)

    ggplot(res[res$GeneSet %in% res2,], aes(x = GeneSet, y = MP, fill = -log10(p.adjust)))+
  geom_tile()+
  theme_bw()+
  theme(panel.grid = element_blank(), aspect.ratio = 1)

    x = unique(res1[res1$GeneSet %in% res2,]$MP)
    hc = hclust(as.dist(1-all.cor.mean[x,x]), method = "average")
    
    res1[match(hc$labels[hc$order], res1$MP),]
    
    Heatmap(all.cor.mean[x,x], cluster_rows = hc, cluster_columns = hc)
}

    ###
sapply(colnames(GDF), function(mp.name){
    #mp.name = "MP_1"
    #mp.name = colnames(MP)[1]
    message(mp.name)
    p.data = res[res$MP == mp.name & res$p.adjust <0.005, ]
    if(nrow(p.data) < 2){return(NULL)}
    if(nrow(p.data) > 20){p.data = p.data[order(p.data$p.adjust, decreasing = F)[1:20],] }
    p.data = p.data[order(p.data$P.value, decreasing = T),] #has to sort, can't use reorder in ggplot otherwise cannot add second axis
    p.data$x = 1:nrow(p.data)
    
    p = ggplot(p.data, aes(x = .data[["x"]], y = -log10(.data[["p.adjust"]]), fill = .data[["GS.source"]]) )+
        geom_bar(stat = "identity", position=position_dodge())+
        #geom_point()+
        theme_bw()+
        theme(panel.grid = element_blank(), legend.position = "none")+
        coord_flip()+
        scale_x_continuous(breaks = p.data$x, labels = p.data$GeneSet,
                           sec.axis = sec_axis(~., breaks = p.data$x, labels = p.data$GS.source))+
        labs(x = "", y = "", title = mp.name)
        
    png(paste0("01.", mp.name, ".enrich.png"), width = 10, height = 15, units = "in", res = 300)
    print(p)
    dev.off()
    NULL
})



####
if(0){
h_gene_sets = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
C2_gene_sets = msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
C5_gene_sets = msigdbr::msigdbr(species = "Homo sapiens", category = "C5")
C8_gene_sets = msigdbr::msigdbr(species = "Homo sapiens", category = "C8")
msigdb.df = data.frame(rbind(h_gene_sets, C2_gene_sets, C5_gene_sets, C8_gene_sets))

msigdb.use = lapply(unique(msigdb.df$gs_name), function(k){
    idx = msigdb.df$gs_name ==k
   data.frame("GS" = k, Symbol = msigdb.df$gene_symbol[idx], Source = msigdb.df$gs_cat[idx])
})
names(msigdb.use) = unique(msigdb.df$gs_name)
saveRDS(msigdb.use, file="Msigdb.HC2C5C8.rds")
rm(msigdb.df)
msigdb.use = readRDS("../01.MP_pHGG_GBM/Msigdb.HC2C5C8.rds")
}

GDF = MP

res2 = lapply(msigdb.use, function(LL){
  #LL = msigdb.use[[22004]] 
  #phyper()
  LL[,2] = toupper(LL[,2])
  LL.source = LL[1,3]
  #message(LL.source)
  tmp.res = lapply(unique(LL$GS), function(gs){
    #gs = unique(LL$GS)[1]
    gs.genes = LL[LL$GS == gs, 2]
    n.gs.in.background = sum(gs.genes %in% background.genes)
    n.not.in.background.genes = length(background.genes) - n.gs.in.background
    ##
    sig.test = apply(GDF, 2, function(x){
      #x = GDF[,1]
      x = x[!is.na(x)]
      n.overlap = sum(x %in% gs.genes)
      phyper(n.overlap -1 , n.gs.in.background,  n.not.in.background.genes, length(x), lower.tail = F)
    })
    y = cbind(gs, sig.test, names(sig.test), LL.source)
    colnames(y) = c("GeneSet", "P.value", "MP", "GS.source")
    rownames(y) = NULL
    y
  })
  tmp.res = data.frame(do.call(rbind, tmp.res))
  tmp.res[,2] <- as.numeric(tmp.res[,2])
  tmp.res <- tmp.res[tmp.res[,2]< 0.01,]
  tmp.res
})
res2 = data.frame(do.call(rbind, res2))
head(res2)
res2$P.value = as.numeric(as.character(res2$P.value))
res2$p.adjust = NA
for(i in unique(res2$MP)){
  idx = which(res2$MP == i)
  res2$p.adjust[idx] = p.adjust(res2$P.value[idx], method = "BH")
}

###get top5 go for each profram
res3 = lapply(unique(res2$MP), function(k){
    x = res2[res2$MP ==k,]
    x = x[x$P.value < 0.01,]
    x = x$GeneSet[order(x$P.value, decreasing = F)]
    if(length(x) > 5) x = x[1:5]
    x
})
unique(unlist(res3))


sapply(colnames(GDF), function(mp.name){
  #mp.name = "MP_1"
  #mp.name = colnames(MP)[1]
  message(mp.name)
  p.data = res2[res2$MP == mp.name & res2$p.adjust <0.005, ]
  if(nrow(p.data) < 2){return(NULL)}
  if(nrow(p.data) > 20){p.data = p.data[order(p.data$p.adjust, decreasing = F)[1:20],] }
  p.data = p.data[order(p.data$P.value, decreasing = T),] #has to sort, can't use reorder in ggplot otherwise cannot add second axis
  p.data$x = 1:nrow(p.data)
  
  p = ggplot(p.data, aes(x = .data[["x"]], y = -log10(.data[["p.adjust"]]), fill = .data[["GS.source"]]) )+
    geom_bar(stat = "identity", position=position_dodge())+
    #geom_point()+
    theme_bw()+
    theme(panel.grid = element_blank(), legend.position = "none")+
    coord_flip()+
    scale_x_continuous(breaks = p.data$x, labels = p.data$GeneSet,
                       sec.axis = sec_axis(~., breaks = p.data$x, labels = p.data$GS.source))+
    labs(x = "", y = "", title = mp.name)
  
  png(paste0("01.msigdb-", mp.name, ".enrich.png"), width = 10, height = 15, units = "in", res = 300)
  print(p)
  dev.off()
  NULL
})




#
#







###keep go
#MP_1: Oligo --C.Greenwald_2024_GBM
#MP_2: Vasc --C.Greenwald_2024_GBM
#MP_3: Neu-NRGN
head(res)
unique(res$MP)
###Select 3-5 representative entichment for each MP
selc.rows = lapply(unique(res$MP), function(mm){
    #mm = unique(res$MP)[1]
    x = res[res$MP == mm,]
    x = x[order(x$p.adjust, decreasing = F),]
    rownames(x)[1:3]
})

res$GeneSet[as.numeric(selc.rows[[6]])]
selc.rows = unique(unlist(selc.rows))

res$GeneSet2 = paste(res$GeneSet, res$GS.source, sep = "|")
res1 = res[selc.rows,]
res1
nrow(res1)

res2 = res[res$GeneSet2 %in% res1$GeneSet2 ,]
##
res2$y1 = as.numeric(factor(res2$GeneSet, levels = unique(res1$GeneSet)))
yy = unique(res2$y1)
yy1 = res2$GeneSet[match(yy, res2$y1)]
yy2 = res2$GS.source[match(yy, res2$y1)]

max(-log10(res2$p.adjust))
min(-log10(res2$p.adjust))
max(res2$p.adjust)
res2$p.adjust[res2$p.adjust<1e-20] = 1e-20

sort(unique(res2$GeneSet2))
res3 = res2[!res2$GeneSet2 %in% c("low.quality|C.Greenwald_2024IDH", "U1|Nowakowski_2017", "Mac|C.Greenwald_2024GBM","macrophage|C.Greenwald_2024IDH",
                                  "MES..malig..|C.Greenwald_2024GBM","AC..malig..|C.Greenwald_2024GBM",
                                  "Oligo|C.Greenwald_2024GBM","oligo|C.Greenwald_2024IDH", "OPC..malig..|C.Greenwald_2024GBM", 
                                  "U1|Nowakowski_2017","vascular|C.Greenwald_2024IDH", "MP41.Unassigned|Gavish_2023", "MGE-IPC2|Nowakowski_2017"),]
sort(unique(res3$GeneSet2))
res3$y1 = as.numeric(factor(res3$GeneSet, levels = sort(unique(res3$GeneSet))))
yy = unique(res3$y1)
yy1 = res3$GeneSet[match(yy, res3$y1)]
yy2 = res3$GS.source[match(yy, res3$y1)]

p = ggplot(res3, aes(x = MP, y = y1, fill = -log10(p.adjust)))+
  geom_tile()+
  geom_vline(xintercept = 1:17, linetype = "dashed", linewidth = 0.1, color = c(rep(c("blue", "gray"), 8), "blue"))+
  geom_hline(yintercept = sort(yy), linetype = "dashed", linewidth = 0.1, color = rep(c("blue", "gray"), length(yy)/2))+
  #scale_fill_gradient2(low= custom_magma[1], high= custom_magma[333], mid = custom_magma[166], midpoint = 10)+
  scale_fill_gradient2(limits = c(2,20), low= custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 15,na.value = "white")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))+
  scale_y_continuous(breaks = yy, labels = yy1,
                     sec.axis = sec_axis(~., breaks = yy, labels = yy2))+
  labs(x = "", y = "")+
  coord_equal(expand = F)
  
p

pdf("01.pHGGGBM_0-25-12_15-15-5/02.MP.enrichment.pdf", width = 9, height = 8, useDingbats = F)
p
dev.off()
 
