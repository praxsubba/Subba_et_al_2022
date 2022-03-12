# Subba_et_al_2022
# RRBS-Workflow
#### Mapping RRBS Reads to the Reference Genome
bsmap -a $input_fastq -d $ref_genome_fasta -o $output_bam -D C-CGG -D T-CGA -w 100 -v 0.08 -r 0 -p 4 -n 0 -s 12 -S 0 -f 5 -q 0 -u -V 2. 

#### Reference genome: GCF_003957565.2_bTaeGut1.4.pri_genomic.fna

### Minimap2 

#### This is used to map cDNA sequences from Dong et al., 2009 to the zebrafinch transcriptome GCF_003957565.2_bTaeGut1.4.pri_rna.fna

minimap2 -I 13G $/GCF_003957565.2_bTaeGut1.4.pri_rna.fna $/sb_array_seq.FASTA

#### This is used to map cDNA sequences from Dong et al., 2009 to the zebrafinch genome GCF_003957565.2_bTaeGut1.4.pri_genomic.fna

minimap2 -I 13G -a --splice --sr --junc-bed $/GCF_003957565.2_bTaeGut1.4.pri_genomic.bed $/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna $/sb_array_seq.FASTA 


```{r}
library("rtracklayer")
library("tidyverse")
readGFF("GCF_003957565.2_bTaeGut1.4.pri_genomic.gff")%>%head()
#I changed this code to select "gene" instead of "mRNA", "lncRNA", "transcript" Now there should be only one line for each gene (instead of the multiple transcripts)
```
```{r}
my_tags <- c("Name", "Dbxref","gene")
my_columns <- c("seqid", "start", "end", "strand", "type")
my_filter<-list(type="gene")
dat<-readGFF("GCF_003957565.2_bTaeGut1.4.pri_genomic.gff",tags=my_tags,columns=my_columns,filter=my_filter)
head(dat)
#I made a new column for interval_start and interval_stop, so I could increase the upstream interval to capture potential regulatory sites
```
```{r}
dat$interval_start<-dat$start
dat$interval_stop<-dat$end
dat
#Add 2kb upstream of the start site gene intervals always have start and stop listed smallest to largest (like in a bed file) for genes on the positive strand, start = start-2000 for genes on the negative strand, stop = stop + 2000
```


```{r}
as.data.frame(dat)%>%mutate(interval_start=ifelse(strand=="+",start-2000,start))%>%  #if positive strand, subtract 2000 from the start
  mutate(interval_start=ifelse(interval_start<1, 1,interval_start))%>%  #if start is within 2kb of start of chr, start interval at 1
  mutate(interval_stop=ifelse(strand=="-",end+2000,end))%>%
  saveRDS("PS_bTaeGut_v2p_gene_intervals_plus_10K_v2.RDS")
RDS_file <- readRDS("PS_bTaeGut_v2p_gene_intervals_plus_10K_v2.RDS") #This reads the RDS file into a data frame test
gene_list <- read.table(file = "unique_genes_AASA.txt") #This reads the txt file in a data frame gene_list from Habtuated vs novel, habituated vs silence or novel vs silence
colnames(gene_list) <- "gene" #This step changes the column name/header to "gene" for the gene_list data frame
RRBS_genes <- dplyr::semi_join(RDS_file, gene_list, by="gene") #This step matches filters the RDS for only the unique or significant genes in gene_list
saveRDS(RRBS_genes, file="Condition_RRBS_genes.RDS") #saves the RDS file
```
