library(ComplexHeatmap)
cor.m = readRDS("04.sc.cor.rds")

p = Heatmap(cor.m, 
            clustering_distance_rows = as.dist(1-cor.m), 
            clustering_distance_columns = as.dist(1-cor.m),
            clustering_method_rows = "average",
            clustering_method_columns = "average",
            show_row_names = F,
            show_column_names = F,
            split = 7)

png("04.sc.cor.png", width = 10, height = 10, units = "in", res = 150)
draw(p)
dev.off()
