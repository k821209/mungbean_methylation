#/home/k821209/programs/eismark_v0.13.1/bismark_genome_preparation --bowtie2 --path_to_bowtie /usr/local/bin/ /home/k821209/mungbean_methylation/gw/
#/home/k821209/programs/bismark_v0.13.1/bismark --gzip --bowtie2 --fastq --bam -p 4 -N 1 --samtools_path /usr/bin/ --path_to_bowtie /usr/local/bin/ /home/k821209/mungbean_methylation/gw/ -1 S_1_4_5_1.fastq.gz -2 S_1_4_5_2.fastq.gz
/home/k821209/programs/bismark_v0.13.1/bismark_methylation_extractor --multicore 3 --genome_folder /home/k821209/mungbean_methylation/gw/ --scaffolds --cytosine_report --CX --paired-end --no_overlap --comprehensive --samtools_path /usr/bin/ --gzip --bedGraph --counts --buffer_size 10G ${1}
