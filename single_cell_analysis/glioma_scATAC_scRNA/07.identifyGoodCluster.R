load("06.ABS.rda")
library(parallel)
all.clusters.n.markers.names = unlist(lapply(all.clusters.n.markers, function(k){ k[[1]]}))
names(all.clusters.n.markers) = all.clusters.n.markers.names;

sig.deg.cnt = lapply(all.clusters.n.markers, function(ll){
  k = ll[[1]];
  deg.list = ll[[3]]
  nn = lapply(deg.list, function(deg){
    n1 = sum(deg$p_val_adj < 0.05 & deg$avg_log2FC > log2(3));
    n2 = sum(deg$p_val_adj < 0.005 & deg$avg_log2FC > log2(3));
    return(c(n1,n2));
  })
  nn = do.call(rbind, nn);
  rownames(nn) = paste0(k, "-", names(deg.list));
  nn;
})

sig.deg.cnt = data.frame(do.call(rbind, sig.deg.cnt))
colnames(sig.deg.cnt) = c("Nsig1", "Nsig2"); #0.05 & 0.005

keep.k = sig.deg.cnt[sig.deg.cnt[,1] >20 & sig.deg.cnt[,2] >15,]
nrow(keep.k)


### Extract the DEGs for the keep.k clusters
keep.k.deg = lapply(rownames(keep.k), function(k.split){
  k = unlist(strsplit(k.split, "-"))[1];
  sp = unlist(strsplit(k.split, "-"))[2];
  n = which(names(all.clusters.n.markers) == k);
  j = which(names(all.clusters.n.markers[[n]][[3]]) == sp);
  deg <- all.clusters.n.markers[[n]][[3]][[j]];
  deg <- deg[deg$p_val_adj <0.05 & deg$avg_log2FC > log2(3), ];
  deg;
})
names(keep.k.deg) = rownames(keep.k);


###calculate the jaccard index between clusters
jaccard.index = lapply(1:(-1+length(keep.k.deg)), function(i){
  res = lapply(c((i+1): length(keep.k.deg)), function(j){
      deg1 = keep.k.deg[[i]];
      deg2 = keep.k.deg[[j]];
      jac = length(intersect(rownames(deg1), rownames(deg2)))/length(union(rownames(deg1), rownames(deg2)));
      c(names(keep.k.deg)[i], names(keep.k.deg)[j], jac)
  })
  res = do.call(rbind, res);
  res;
})

jaccard.index = data.frame(do.call(rbind, jaccard.index), stringsAsFactors = F)
colnames(jaccard.index) = c("c1", "c2", "jac")
jaccard.index[,3] = as.numeric(jaccard.index[,3]);

jaccard.index$c1.Nsig1 = keep.k$Nsig1[match(jaccard.index$c1, rownames(keep.k))];
jaccard.index$c2.Nsig1 = keep.k$Nsig1[match(jaccard.index$c2, rownames(keep.k))];
jaccard.index$c1.Nsig2 = keep.k$Nsig2[match(jaccard.index$c1, rownames(keep.k))];
jaccard.index$c2.Nsig2 = keep.k$Nsig2[match(jaccard.index$c2, rownames(keep.k))];

### (3) For each pair of clusters with Jaccard index above 75%, we excluded the cluster with lower Nsig
exclude.cluster = mclapply(1:nrow(jaccard.index), function(k){
    if(jaccard.index$jac[k] > 0.75){
        if(jaccard.index$c1.Nsig1 > jaccard.index$c2.Nsig1){ #filter c2
            return(jaccard.index$c2[k])
        }else if(jaccard.index$c1.Nsig1 < jaccard.index$c2.Nsig1){
            return(jaccard.index$c2[k])
        }else if(jaccard.index$c1.Nsig1 == jaccard.index$c2.Nsig1 ){ #filter the large cluster?
            #the large cluster is always in the first row(c1);
            if(jaccard.index$c1.Nsig2 == jaccard.index$c2.Nsig2) message("Unluckly find a pair with same Nsig1 and Nsig2 ", rownames(jaccard.index)[k])
            ifelse(jaccard.index$c1.Nsig2>jaccard.index$c2.Nsig2, return(jaccard.index$c2[k]), return(jaccard.index$c2[1]))
        }else{ #jaccard index
            stop("Error: unkown jaccard index: ", k, " ", paste0(jaccard.index[k,], collapse = " "));
        }
    }else{
      return(NULL);
    }
}, mc.cores = 2);

final.exclude.cluster = unique(unlist(exclude.cluster));
final.keep.cluster = keep.k[!rownames(keep.k) %in% final.keep.cluster,]
nrow(final.keep.cluster);


### Extract the clusters
final.clusters = lapply(rownames(final.keep.cluster), function(k.split){
  k = unlist(strsplit(k.split, "-"))[1];
  sp = unlist(strsplit(k.split, "-"))[2];
  n = which(names(all.clusters.n.markers) == k);
  ind = names(which(all.clusters.n.markers[[1]][[2]] == 9));
  ind;
})
names(final.clusters) = rownames(final.keep.cluster)

final.degs = keep.k.deg[rownames(final.keep.cluster)]

save(final.exclude.cluster, final.keep.cluster, final.clusters,final.degs, "test.rda")
q(save = "no")


