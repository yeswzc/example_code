source("score.fun.R")

extractProgram = function(file, k.gene = 100){
    spectra = read.table(file, head =T, sep = "\t", check.names= F, row.names = 1)
    pg = lapply(1:nrow(spectra),function(k){
        G = colnames(spectra)[order(spectra[k,], decreasing = T, na.last = T)[1:k.gene]];
    });
    #names(pg) = paste0(name , ":", 1:length(pg));
    pg;
};

score.cells = function(data, programs, n.control = 100){
    res = lapply(1:ncol(programs), function(k){
            G = programs[,k]
            G = G[G %in% rownames(data)];
            message("Expression program ", k, " has ", length(G), " genes...")
            #score = colMeans(data[G,]);
            score = gene_score(score.G = G, E = data, control.unit = n.control);
            score;
        });
    res = do.call(cbind, res);
    res;
}

####
library(Seurat)
all.spectra = list.files("./", "cNMF.spectra.*consensus.txt", recursive=T);
all.programs = lapply(all.spectra, function(file.name){
    x1 = unlist(strsplit(file.name, "/"))
    x2 = unlist(strsplit(x1[1], "\\-"));
    sample.ID = x2[2]
    message(file.name, ":", sample.ID)
    p = extractProgram(file.name);
    p = do.call(cbind, p);
    colnames(p) = paste0(sample.ID, "-", 1:ncol(p));
    p;
})
all.programs = do.call(cbind, all.programs);

write.table(all.programs, file = "03.cNMF.allPrograms.txt", row.names =F, quote =F, sep = "\t");


all.program.scores = lapply(all.spectra, function(file.name){
    x1 = unlist(strsplit(file.name, "/"))
    x2 = unlist(strsplit(x1[1], "\\-"));
    sample.ID = x2[2]
    count.file = paste0("02.macrophage.cNMF-", sample.ID, ".tsv.gz")
    #02.macrophage.cNMF-CCH.tsv.gz 
    data = t(read.table(count.file, head =T, check.names =F, row.names = 1));
    colnames(data) = paste0(sample.ID, "-", colnames(data));
    data = CreateSeuratObject(data, min.cells = 0, min.features = 0)
    data <- NormalizeData(data)
    #data <- ScaleData(data, features = rownames(data));
    data = data[["RNA"]]@data
    message("Program ", sample.ID, "-")
    scores = score.cells(data, all.programs);
    colnames(scores) = colnames(all.programs);
    scores;
})
all.program.scores = do.call(rbind, all.program.scores);
print(head(all.program.scores[,1:4]));


saveRDS(all.program.scores, file = "03.cNMF.allPrograms.scores.rds");




