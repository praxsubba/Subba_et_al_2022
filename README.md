# Subba_et_al_2022
# RRBS-Workflow
### Mapping RRBS Reads to the Reference Genome using BSMAP
bsmap -a $input_fastq -d $ref_genome_fasta -o $output_bam -D C-CGG -D T-CGA -w 100 -v 0.08 -r 0 -p 4 -n 0 -s 12 -S 0 -f 5 -q 0 -u -V 2. 

#### Reference genome: GCF_003957565.2_bTaeGut1.4.pri_genomic.fna

### Minimap2 

#### This is used to map cDNA sequences from Dong et al., 2009 to the zebrafinch transcriptome GCF_003957565.2_bTaeGut1.4.pri_rna.fna (Available at: https://www.ncbi.nlm.nih.gov/assembly/GCF_003957565.2/)

minimap2 -I 13G ~/GCF_003957565.2_bTaeGut1.4.pri_rna.fna ~/sb_array_seq.FASTA

#### This is used to map cDNA sequences from Dong et al., 2009 to the zebrafinch genome GCF_003957565.2_bTaeGut1.4.pri_genomic.fna

minimap2 -I 13G -a --splice --sr --junc-bed ~/GCF_003957565.2_bTaeGut1.4.pri_genomic.bed ~/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna ~/sb_array_seq.FASTA 

## Rtracklayer to capture genomic internals, promoter regions and transcription start sites.

```{r}
hab_vs_nov_sig_genes <- read.csv(file="AASA_Manuscript_plot.csv") 
hab_vs_nov_sig_genes_final <- hab_vs_nov_sig_genes %>% dplyr::select("gene","clone_id","AASA_fold_clone","AASA_fdr_clone","log2FoldChange","fdr")  %>% filter(AASA_fdr_clone < 0.05) 
```

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
#AASA
#write to file
as.data.frame(dat)%>%mutate(interval_start=ifelse(strand=="+",start-2000,start))%>%  #if positive strand, subtract 2000 from the start
  mutate(interval_start=ifelse(interval_start<1, 1,interval_start))%>%  #if start is within 2kb of start of chr, start interval at 1
  mutate(interval_stop=ifelse(strand=="-",end+2000,end))%>%
  saveRDS("PS_bTaeGut_v2p_gene_intervals_plus_10K_v2.RDS")
test <- readRDS("PS_bTaeGut_v2p_gene_intervals_plus_10K_v2.RDS") #This reads the RDS file into a data frame test
unique_genes_AASA <- hab_vs_nov_sig_genes_final #This reads the txt file in a data frame unique_genes

TEST_2 <- dplyr::semi_join(test, unique_genes_AASA, by="gene") #This step matches filters the RDS for only the unique or significant genes in unique_genes
saveRDS(TEST_2, file="August_2022_sig_Habituated_vs_novel_Rtracklayer.RDS") #saves the RDS file
write.table(TEST_2,file="August_2022_sig_Habituated_vs_novel_Rtracklayer.txt",sep = "\t",
            row.names = FALSE) #this txt file can be used as a BED file in the following steps
```

### Extracting alignments from BSMAP output files using samtools view

samtools view -b -L August_2022_sig_Habituated_vs_novel_Rtracklayer.bed ~/$_output.bam > $_interval.bam

### Extracting Methylated CpG sites using BSMAPz
python ~/methratio.py -o $_methratio.txt -d ~/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna -z -x CG ~/$_interval.bam

## Methylkit for differential methylation analysis, read coverage filtration, and annotation
#### This is an example that shows Habituated vs Silence
#### FA = Familiar, SI = Silence, NO = Novel

---
title: "August 2022 AASS Methylkit"
output:
  html_document:
    df_print: paged
---

```{r}
setwd("/Users/subbaprakrit/Library/CloudStorage/Box-Box/Prakrit Subba research/RRBS_July/February_Manusript/Methratio_files/")
setwd("~/Box/Prakrit Subba research/RRBS_July/February_Manusript/Methratio_files/") #mac pc
```

```{r}
#Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylKit")
library("methylKit")
file.list = list("G1FA160_Aug_2022_new_workflow_methratio.txt", "G1FA211_Aug_2022_new_workflow_methratio.txt", "G1FA237_Aug_2022_new_workflow_methratio.txt","G2FA149_Aug_2022_new_workflow_methratio.txt","G2FA209_Aug_2022_new_workflow_methratio.txt","G2FA222_Aug_2022_new_workflow_methratio.txt","G1SI152_Aug_2022_new_workflow_methratio.txt","G1SI199_Aug_2022_new_workflow_methratio.txt","G1SI246_Aug_2022_new_workflow_methratio.txt","G2SI146_Aug_2022_new_workflow_methratio.txt","G2SI188_Aug_2022_new_workflow_methratio.txt","G2SI218_Aug_2022_new_workflow_methratio.txt")
myobj=methRead( file.list,pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2, coverage.col=6,strand.col=3,freqC.col=5 ),
                sample.id=list("G1FA160", "G1FA211", "G1FA237", "G2FA149", "G2FA209", "G2FA222", "G1SI152", "G1SI199", "G1SI246", "G2SI146", "G2SI188", "G2SI218"),assembly="taeGut1",treatment=c(1,1,1,1,1,1,0,0,0,0,0,0))
```
#get methylation statistics
```{r}
for (i in 1:12) {
  getMethylationStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
}
```

```{r}
#filter by coverage and get coverage stats

filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=NULL)

for (i in 1:12){
  getCoverageStats(filtered.myobj[[i]],plot=TRUE,both.strands=FALSE)
}

for (i in 1:12) {
  getMethylationStats(filtered.myobj[[i]],plot=TRUE,both.strands=FALSE)
}
```

```{r}
meth=unite(filtered.myobj, destrand=FALSE)
meth
meth_destrand=unite(filtered.myobj, destrand=TRUE) 
##see that destrand=TRUE increases coverage
meth_destrand 

```
```{r}
#meth_unfiltered=unite(myobj, destrand=FALSE)
#meth_unfiltered
#meth_unfiltered_destrand=unite(myobj, destrand=TRUE) 
##see that destrand=TRUE increases coverage
#meth_unfiltered_destrand 

#myDiff_unfiltered_d=calculateDiffMeth(meth_unfiltered_destrand)
#getData(myDiff_unfiltered_d) %>% dplyr::arrange(qvalue)
```

```{r}
#see how samples cluster

clusterSamples(meth, dist="correlation", method="ward", plot=TRUE,sd.threshold = .90)
clusterSamples(meth_destrand, dist="correlation", method="ward", plot=TRUE,sd.threshold = .90)
```
```{r}
#plot Principal Components

PCASamples(meth,adj.lim = c(.4, 1), sd.threshold = .90) #not sure about how to set this value?
PCASamples(meth_destrand,adj.lim = c(.5, 1), sd.threshold = .90) #not sure about how to set this value?

```
```{r}
#get a methylDiff object containing the differential methylation statistics and locations for regions or bases
myDiff=calculateDiffMeth(meth)
#write.csv(getData(myDiff),"summary.csv")
getData(myDiff)
myDiff_d=calculateDiffMeth(meth_destrand)
getData(myDiff_d) %>% dplyr::arrange(qvalue)%>% write.csv(file="Aug_2022_new_workflow_meth_diff.csv")
```

```{r}
# get all differentially methylated bases, difference greater than 25
myDiff10p=getMethylDiff(myDiff,qvalue=.01)
getData(myDiff10p)


myDiff10p_destranded=getMethylDiff(myDiff_d,qvalue=.01,difference=25)
diff_25_august <- getData(myDiff10p_destranded) %>% dplyr::arrange(qvalue)
write.csv(diff_25_august,"Aug_2022_new_workflow_meth_diff25.csv")

```
#Annotation Plots
```{r}
library(genomation)
```

```{r}
# read the gene BED file
bed=readTableFast("1_test_genePred.bed.txt",header=FALSE,skip="auto")
head(bed)
gene.obj=readTranscriptFeatures("1_test_genePred.bed.txt",remove.unusual=FALSE)
```

```{r}
#
# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data
#
annotateWithGeneParts(as(myDiff10p_destranded,"GRanges"),gene.obj)
```

```{r}
promoters=regionCounts(filtered.myobj,gene.obj$promoters)

head(promoters[[1]])

```

```{r}
diffAnn=annotateWithGeneParts(as(myDiff10p_destranded,"GRanges"),gene.obj)

# target.row is the row number in myDiff10p_destranded
head(getAssociationWithTSS(diffAnn))
```
```{r}
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
```
```{r}
plotTargetAnnotation(diffAnn,precedence=TRUE,
    main="Familiar vs Silence Differential Methylation Annotation")
```




```{r}
sessionInfo()
```

