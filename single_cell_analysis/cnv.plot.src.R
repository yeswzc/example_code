library(ComplexHeatmap)
red_blue_50 = c("#0a3b70","#10457e","#15508d","#1b5a9c","#2065ab","#276eb0","#2e77b5",
                "#3480b9","#3b88be","#4291c2","#4f9bc7","#5fa5cd","#6eaed2","#7eb8d7",
                "#8dc2dc","#9bc9e0","#a7d0e4","#b3d6e8","#c0dceb","#cce2ef","#d5e7f1",
                "#ddebf2","#e4eef4","#ecf2f5","#f3f5f6","#f8f4f2","#f9efe9","#fae9df",
                "#fbe4d6","#fcdecd","#fcd7c2","#fbccb4","#f9c2a7","#f7b799","#f5ac8b",
                "#f2a17f","#ec9374","#e6866a","#e17860","#db6b55","#d55d4c","#ce4f45",
                "#c6413e","#bf3338","#b82531","#b1182b","#a21328","#930e26","#840924",
                "#760521");

findsplitK = function(cell.cluster, N = 200, k.init = 2, k.max = 500, n.cut = 2){
  for(kk in k.init:k.max){
    t = sort(table(cutree(cell.cluster, k = kk)), decreasing = T)
    if(t[n.cut] > N){ cat("\nTry K = ",kk,"\n"); break; }else{cat("K", kk,":", t[n.cut],"...")}
  }
  return(NULL)
}

findsplitKrev = function(cell.cluster, N = 200, k.min = 2, k.max = 200, n.cut = 2){
  for(kk in c(k.max:k.min)){
    t = sort(table(cutree(cell.cluster, k = kk)), decreasing = T)
    cat(t[1:4])
    if(t[n.cut] > N){ cat("\nTry K = ",kk,"\n"); break; }else{cat("K", kk,":", t[n.cut],"...")}
  }
  return(NULL)
}

subtractControl = function(cnv, cell.cluster.cut, k.control, choice = "mean"){
  #choice = "median"
  #k.control = 2
  if(choice == "median"){
    cnv.control = apply(cnv[, which(cell.cluster.cut == k.control)], 1, median)
  }else{
    cnv.control = rowMeans(cnv[,which(cell.cluster.cut== k.control)])
  }
  cnv2 = apply(cnv, 2, function(x) x - cnv.control)
  return(cnv2)
}

make_raw_cnv_plot = function(cnv, pos, cell.cluster, png.file){
  #cnv.cor = cor(cnv, method = "pearson")
  #d = 1 - as.dist(cnv.cor)
  #cell.cluster = hclust(d, method = "average")
  
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))
  row.names(cnv) = paste0(pos$seqnames,":", pos$start, "-", pos$end)
  p = Heatmap(t(cnv), cluster_columns = F,
              row_dend_gp = gpar(lwd = 0.2),
              show_row_names = F,show_column_names = F,
              cluster_rows = cell.cluster, show_row_dend = T, 
              column_split = pos$seqnames, cluster_column_slices = F,
              column_title_rot = 90, column_title_side = "bottom",
              border = TRUE, column_gap = unit(0.2, "mm"),
              name = "cnv", col = col_fun,
              use_raster = TRUE, raster_quality = 2,
              heatmap_legend_param = list(direction = "horizontal"))
  
  #EGFR:chr7:55019017-55211628 
  #PDGFRA: chr4:54229293-54298245
  n.chr = length(unique(pos$seqnames))/2;
  n.chr = n.chr + (1+n.chr)/10;
  x1 = which(pos$seqnames==4 & pos$end >= 54229293)[1] #PDGFRA
  x2 = which(pos$seqnames==7 & pos$end >= 55019017)[1] #EGFR
  
  cat(x1/nrow(pos), " ",x2/nrow(pos)," ", n.chr,"\n")
  
  png(file = png.file, width = 10, height = 5, units = "in", res = 100)
  draw(p, heatmap_legend_side = "bottom", merge_legend = FALSE, annotation_legend_side ="bottom")
  decorate_heatmap_body("cnv", {
    grid.lines(x= n.chr* c(x1, x1)/nrow(pos), y = c(0, 1), gp = gpar(lty = 3, lwd = 1))
    grid.lines(x= n.chr* c(x2, x2)/nrow(pos), y = c(0, 1), gp = gpar(lty = 3, lwd = 1))
    #grid.lines(x= c(0, n.chr), y = c(0, 1), gp = gpar(lty = 3, lwd = 3))
    
    #grid.lines(x= c(15, 15)/nrow(pos), y = c(0, 1), gp = gpar(lty = 3, lwd = 1), default.units = "native")
  })
  dev.off()
  return(p)
}


###
make_raw_cnv_plotv2 = function(cnv, pos, cell.cluster, png.file, res.cluster){
  #cnv.cor = cor(cnv, method = "pearson")
  #d = 1 - as.dist(cnv.cor)
  #cell.cluster = hclust(d, method = "average")
  cell.cluster.col = RColorBrewer::brewer.pal(12,"Paired")[1:length(unique(res.cluster))]
  #print(cell.cluster.col)
  names(cell.cluster.col) = unique(res.cluster)
  names(cell.cluster.col)[is.na(names(cell.cluster.col))] = "NA"
  #print(cell.cluster.col)
  
  left_annotation = rowAnnotation(cluster = factor(res.cluster, exclude = NULL), 
                                  col = list(cluster = cell.cluster.col),
                                  annotation_legend_param = list(
                                    cluster = list(direction = "horizontal", nrow = 1)))
  
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))
  p = Heatmap(t(cnv), cluster_columns = F,
              row_dend_gp = gpar(lwd = 0.2),
              show_row_names = F,
              left_annotation = left_annotation,
              cluster_rows = cell.cluster, show_row_dend = T, 
              column_split = pos$seqnames, cluster_column_slices = F,
              column_title_rot = 90, column_title_side = "bottom",
              border = TRUE, column_gap = unit(0.2, "mm"),
              name = " ", col = col_fun,
              use_raster = TRUE, raster_quality = 2,
              heatmap_legend_param = list(direction = "horizontal"))
  png(file = png.file, width = 10, height = 5, units = "in", res = 100)
  draw(p, heatmap_legend_side = "bottom", merge_legend = FALSE, annotation_legend_side ="bottom")
  dev.off()
  return(p)
}


make_raw_cnv_cluster_plot = function(cnv, pos, png.file, distance.method = "euclidean", cluster.method = "complete"){
  cell.hclust = hclust(dist(t(cnv), method = distance.method), method = cluster.method)
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))
  p = Heatmap(t(cnv), cluster_columns = F,
              row_dend_gp = gpar(lwd = 0.2),
              show_row_names = F,
              cluster_rows = cell.hclust,
              show_row_dend = T, 
              column_split = pos$seqnames, cluster_column_slices = F,
              column_title_rot = 90, column_title_side = "bottom",
              border = TRUE, column_gap = unit(0.2, "mm"),
              name = " ", col = col_fun,
              use_raster = TRUE, raster_quality = 2,
              heatmap_legend_param = list(direction = "horizontal"))
  png(file = png.file, width = 10, height = 5, units = "in", res = 100)
  draw(p, heatmap_legend_side = "bottom", merge_legend = FALSE, annotation_legend_side ="bottom")
  dev.off()
  return(cell.hclust)
}



make_rawOrdered_cnv_plot = function(cnv, pos, cell.cluster, png.file){
  #cnv.cor = cor(cnv, method = "pearson")
  #d = 1 - as.dist(cnv.cor)
  #cell.cluster = hclust(d, method = "average")
  cnv = cnv[,cell.cluster$order]
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))
  p = Heatmap(t(cnv), cluster_columns = F,
              row_dend_gp = gpar(lwd = 0.2),
              show_row_names = F,
              cluster_rows = F, show_row_dend = T, 
              column_split = pos$seqnames, cluster_column_slices = F,
              column_title_rot = 90, column_title_side = "bottom",
              border = TRUE, column_gap = unit(0.2, "mm"),
              name = " ", col = col_fun,
              use_raster = TRUE, raster_quality = 2,
              heatmap_legend_param = list(direction = "horizontal"))
  png(file = png.file, width = 10, height = 5, units = "in", res = 100)
  draw(p, heatmap_legend_side = "bottom", merge_legend = FALSE, annotation_legend_side ="bottom")
  dev.off()
  return(p)
}

make_leftAnn_cnvplot = function(cnv, pos, png.file, cell.cluster, cell.cluster.cut, plot.dend = F ){
  #k = length(unique(cell.cluster.cut))
  cell.cluster.col = RColorBrewer::brewer.pal(8,"Paired")[unique(cell.cluster.cut)]
  names(cell.cluster.col) = unique(cell.cluster.cut)
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))
  hm = Heatmap(t(cnv), cluster_columns = F,
               cluster_rows = cell.cluster, show_row_dend = plot.dend, show_row_names = F,
               left_annotation = rowAnnotation(cluster = factor(cell.cluster.cut), 
                                               col = list(cluster = cell.cluster.col),
                                               annotation_legend_param = list(
                                                 cluster = list(direction = "horizontal", nrow = 1))),
               column_split = pos$seqnames, cluster_column_slices = F,
               column_title_rot = 90, column_title_side = "bottom",
               border = TRUE, column_gap = unit(0.2, "mm"),
               name = "cnv", col = col_fun,
               use_raster = TRUE, raster_quality = 2,
               heatmap_legend_param = list(direction = "horizontal"));
  #EGFR:chr7:55019017-55211628 
  #PDGFRA: chr4:54229293-54298245
  n.chr = length(unique(pos$seqnames))/2;
  n.chr = n.chr + (1+n.chr)/10;
  x1 = which(pos$seqnames==4 & pos$end >= 54229293)[1] #PDGFRA
  x2 = which(pos$seqnames==7 & pos$end >= 55019017)[1] #EGFR
  
  cat(x1/nrow(pos), " ",x2/nrow(pos)," ", n.chr,"\n")
  
  png(png.file, width = 10, height = 5, units = "in", res = 100)
  draw(hm, heatmap_legend_side = "bottom", merge_legend = FALSE, annotation_legend_side ="bottom")
  decorate_heatmap_body("cnv", {
    grid.lines(x= n.chr* c(x1, x1)/nrow(pos), y = c(0, 1), gp = gpar(lty = 3, lwd = 1))
    grid.lines(x= n.chr* c(x2, x2)/nrow(pos), y = c(0, 1), gp = gpar(lty = 3, lwd = 1))
    #grid.lines(x= c(0, n.chr), y = c(0, 1), gp = gpar(lty = 3, lwd = 3))
    #grid.lines(x= c(15, 15)/nrow(pos), y = c(0, 1), gp = gpar(lty = 3, lwd = 1), default.units = "native")
  })
  
  dev.off()
  
  return(hm)
}

make_final_cnv_plot = function(cnv2, pos, png.file, cell.cluster.cut, k.malignant){
  #cnv3 = cnv2[,cell.cluster$order[c(which(cell.cluster.cut.reoder ==2), which(cell.cluster.cut.reoder !=2))]] #order by dend
  cnv3 = cnv2[,c(which(cell.cluster.cut == k.malignant), which(cell.cluster.cut != k.malignant))]
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))
  hm = Heatmap(t(cnv3), cluster_columns = F, cluster_rows = F,
               show_row_names = F,
               column_split = pos$seqnames, cluster_column_slices = F,
               column_title_rot = 90, column_title_side = "bottom",
               border = TRUE, column_gap = unit(0.2, "mm"),
               name = " ", col = col_fun,
               use_raster = TRUE, raster_quality = 2,
               heatmap_legend_param = list(direction = "horizontal"))
  png(png.file, width = 10, height = 5, units = "in", res = 100)
  draw(hm, heatmap_legend_side = "bottom", merge_legend = FALSE, annotation_legend_side ="bottom")
  dev.off()
  return(hm)
  return(hm)
}

make_final_cnv_plot_Dendreorder = function(cnv2, pos, png.file, cell.cluster, cell.cluster.cut, k.malignant, cluster.ann = F){
  cell.cluster.cut.reoder = cell.cluster.cut[cell.cluster$order]  #order of cells after reorder
  cnv3 = cnv2[,cell.cluster$order[c(which(cell.cluster.cut.reoder == k.malignant), 
                                    which(cell.cluster.cut.reoder !=k.malignant))]] #order by dend
  cell.cluster.cut = cell.cluster.cut[cell.cluster$order[c(which(cell.cluster.cut.reoder == k.malignant), 
                                                           which(cell.cluster.cut.reoder !=k.malignant))]]
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))
  if(cluster.ann){
    cell.cluster.col = RColorBrewer::brewer.pal(8,"Paired")[unique(cell.cluster.cut)]
    names(cell.cluster.col) = unique(cell.cluster.cut)
    hm = Heatmap(t(cnv3), cluster_columns = F, cluster_rows = F,
                 show_row_names = F,
                 left_annotation = rowAnnotation(cluster = factor(cell.cluster.cut), 
                                                 col = list(cluster = cell.cluster.col),
                                                 annotation_legend_param = list(
                                                   cluster = list(direction = "horizontal", nrow = 1))),
                 column_split = pos$seqnames, cluster_column_slices = F,
                 column_title_rot = 90, column_title_side = "bottom",
                 border = TRUE, column_gap = unit(0.2, "mm"),
                 name = " ", col = col_fun,
                 use_raster = TRUE, raster_quality = 2,
                 heatmap_legend_param = list(direction = "horizontal"))
  }else{
    hm = Heatmap(t(cnv3), cluster_columns = F, cluster_rows = F,
                 show_row_names = F,
                 column_split = pos$seqnames, cluster_column_slices = F,
                 column_title_rot = 90, column_title_side = "bottom",
                 border = TRUE, column_gap = unit(0.2, "mm"),
                 name = " ", col = col_fun,
                 use_raster = TRUE, raster_quality = 2,
                 heatmap_legend_param = list(direction = "horizontal"))
  }
  png(png.file, width = 10, height = 5, units = "in", res = 100)
  draw(hm, heatmap_legend_side = "bottom", merge_legend = FALSE, annotation_legend_side ="bottom")
  dev.off()
  return(hm)
}


make_final_cnv_plot_Dendreorderv2 = function(cnv2, pos, png.file, cell.cluster, cell.cluster.cut, k.malignant, res.cluster){
  cell.cluster.cut.reoder = cell.cluster.cut[cell.cluster$order]  #order of cells after reorder
  cnv3 = cnv2[,cell.cluster$order[c(which(cell.cluster.cut.reoder == k.malignant), 
                                    which(cell.cluster.cut.reoder !=k.malignant))]] #order by dend
  cell.cluster.cut = cell.cluster.cut[cell.cluster$order[c(which(cell.cluster.cut.reoder == k.malignant), 
                                                           which(cell.cluster.cut.reoder !=k.malignant))]]
  res.cluster = res.cluster[cell.cluster$order[c(which(cell.cluster.cut.reoder == k.malignant), 
                                                 which(cell.cluster.cut.reoder !=k.malignant))]]
  cnv.cell.cluster.col = RColorBrewer::brewer.pal(8,"Paired")[unique(cell.cluster.cut)]
  names(cnv.cell.cluster.col) = unique(cell.cluster.cut)
  
  
  cell.cluster.col = RColorBrewer::brewer.pal(12,"Paired")[1:length(unique(res.cluster))]
  #print(cell.cluster.col)
  names(cell.cluster.col) = unique(res.cluster)
  names(cell.cluster.col)[is.na(names(cell.cluster.col))] = "NA"
  #print(cell.cluster.col)
  
  left_annotation = rowAnnotation(cluster = factor(res.cluster, exclude = NULL), 
                                  cnv.cluster = factor(cell.cluster.cut),
                                  col = list(cluster = cell.cluster.col, cnv.cluster = cnv.cell.cluster.col),
                                  annotation_legend_param = list(
                                    cluster = list(direction = "horizontal", nrow = 1)))
  
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#0a3b70", "white", "#760521"))
  hm = Heatmap(t(cnv3), cluster_columns = F, cluster_rows = F,
               show_row_names = F,
               left_annotation = left_annotation,
               column_split = pos$seqnames, cluster_column_slices = F,
               column_title_rot = 90, column_title_side = "bottom",
               border = TRUE, column_gap = unit(0.2, "mm"),
               name = " ", col = col_fun,
               use_raster = TRUE, raster_quality = 2,
               heatmap_legend_param = list(direction = "horizontal"))
  
  png(png.file, width = 10, height = 5, units = "in", res = 100)
  draw(hm, heatmap_legend_side = "bottom", merge_legend = FALSE, annotation_legend_side ="bottom")
  dev.off()
  return(hm)
}



#correlation heatmap
make_correlation_heatmap = function(cnv, png.file, barcodes, cluster.method = "average", correlation.method = "pearson", plot.dend = T){
  #correlation.method = "pearson"
  #cluster.method = "average"
  cnv.cor = cor(cnv, method = correlation.method)
  cnv.cluster = hclust(as.dist((1 - cnv.cor)/2), method = cluster.method)
  #cnv.cluster = hclust(as.dist((1 - cnv.cor) method = cluster.method)
  #return(cnv.cluster)
  umi = log10(barcodes$passed_filters)
  mito.ratio = 100* barcodes$mitochondrial/barcodes$passed_filters
  
  col_fun = circlize::colorRamp2(c(min(umi), median(umi), max(umi)), c("#0a3b70", "white", "#760521"))
  promoter_ratio =  (barcodes$promoter_region_fragments +1) / (barcodes$passed_filters +1)
  col_fun.promoter = circlize::colorRamp2(c(min(promoter_ratio), median(promoter_ratio), 
                                            max(promoter_ratio)), c("#0a3b70", "white", "#760521"))
  
  col_fun.mito = circlize::colorRamp2(c(min(mito.ratio), median(mito.ratio), 
                                        max(mito.ratio)), c("#0a3b70", "white", "#760521"))
  
  row.mean = as.vector(colMeans(abs(cnv)))
  row.mean.mean = mean(row.mean)
  #third quartile + 1.5·IQR 
  fun.max = quantile(row.mean)[4] + 1.5*IQR(row.mean.mean)
  fun.max = ifelse(fun.max < max(row.mean), fun.max, max(row.mean))
  col_fun2 = circlize::colorRamp2(c(min(row.mean), row.mean.mean , fun.max), c("blue", "white", "red"))
  
  #add a ruller ...
  #k = c(1:ncol(cnv))[match(cnv.cluster$order, 1:ncol(cnv))]
  
  row.ann = rowAnnotation(foo = anno_mark(cnv.cluster$order[round(quantile(1:ncol(cnv), probs = seq(0.05,1,0.05)))], 
                                          labels = seq(0.5,10, by = 0.5), labels_gp = gpar(fontsize = 5), 
                                          link_width = unit(0, "mm"), link_height = unit(0, "mm")), 
                          umi = umi, pr = promoter_ratio, cn = row.mean, 
                          col = list(umi = col_fun, pr = col_fun.promoter, cn = col_fun2))
  
  hm = Heatmap(cnv.cor, cluster_rows = cnv.cluster, cluster_columns = cnv.cluster, 
               right_annotation = row.ann,
               row_dend_gp = gpar(lwd = 0.2),
               show_row_dend = plot.dend,
               show_column_dend = F,
               row_dend_width = unit(12, "cm"), 
               show_row_names = F, show_column_names = F, col = red_blue_50, name = "pcc",
               use_raster = T)
  w.png = 10
  if(! plot.dend){ w.png = 5}
  png(file = png.file, width = w.png, height = 4, units = "in", res = 100)
  draw(hm, annotation_legend_side = "right")
  dev.off()
  return(cnv.cluster)
}

