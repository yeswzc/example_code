cellranger vdj --id 11thhighavid-VDJ --fastqs /home/lij36/lij36/02.single_cell/00_fastq/ --sample 11thhighavid-VDJ --reference /fdb/cellranger/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0/ --jobmode local --localcores 10 --localmem 64
cellranger vdj --id 11thlowavid-VDJ --fastqs /home/lij36/lij36/02.single_cell/00_fastq/ --sample 11thlowavid-VDJ --reference /fdb/cellranger/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0/ --jobmode local --localcores 10 --localmem 64
cellranger vdj --id 6thhighavid-VDJ --fastqs /home/lij36/lij36/02.single_cell/00_fastq/ --sample 6thhighavid-VDJ --reference /fdb/cellranger/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0/ --jobmode local --localcores 10 --localmem 64
cellranger vdj --id 6thlowavid-VDJ --fastqs /home/lij36/lij36/02.single_cell/00_fastq/ --sample 6thlowavid-VDJ --reference /fdb/cellranger/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0/ --jobmode local --localcores 10 --localmem 64
cellranger aggr --id VDJagge --csv vdj.pb.csv
