ls | grep "5$" |xargs -i echo "velocyto run10x -@ 2 -m mm10_rmsk.gtf {}/ /fdb/cellranger/refdata-cellranger-mm10-3.0.0/genes/genes.gtf" |xargs -i quick_sbatch -m velocyto -g 10g -h 8h '{}'
