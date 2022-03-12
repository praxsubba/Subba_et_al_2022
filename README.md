# Subba_et_al_2022
# RRBS-Workflow
#### Mapping RRBS Reads to the Reference Genome
bsmap -a $input_fastq -d $ref_genome_fasta -o $output_bam -D C-CGG -D T-CGA -w 100 -v 0.08 -r 0 -p 4 -n 0 -s 12 -S 0 -f 5 -q 0 -u -V 2. 

Reference genome: GCF_003957565.2_bTaeGut1.4.pri_genomic.fna

