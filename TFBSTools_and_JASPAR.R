#================================================================================
#
#       Using 'TFBSTools' and "JASPAR2018' packages to detect TFBS
#
#================================================================================

## Load Packages
library(TFBSTools)
library(JASPAR2018)
library(Biostrings)

## Set working directory
setwd('/Users/mijiarui/biosoft/HOMER/results/significant_delta_vs_alphaBeta_lincRNA_2')

## Prepare two file
### For TFBS detection, there are two files needs to be prepared. The first is a PWM (position weighted matrix)
### which has been stored in JASPAR2018, the other is fasta file of the regions you are interested in
### For fasta file, you can use bedtools getfasta function to get the sequence from the coordinate (you
### also need to prepare the genome fasta sequence file as one of the input for bedtools getfasta).
### The bedtools code is like this: bedtools getfasta -fo black.fa -name -fi 
### /Users/mijiarui/RNA-seq/danio_genome/Danio_rerio.GRCz10.fa -bed promoter.txt

### Make an 'opts' object to help you retrieve the PWM from JASPAR2028
### First let's check which category you want to test. Here I chose 'tax_group' and vertebrates, but of course, you can
### chose other options, especially 'species' and the taxonomy_id of the species (Check it in NCBI taxonomy track)
?getMatrixSet
opts = list()
opts[['species']] <- 7227
opts[['matrixtype']] <- 'PWM'
PWMatrixList <- getMatrixSet(JASPAR2018, opts)

### Make you sequence object
dnaSet <- readDNAStringSet('promoter_significant_delta_vs_alphaBeta_lincRNA.fa', format = 'fasta')
dnaTab <- as.data.frame(dnaSet)

### For loop to get the results. The key function is 'searchSeq()' which contains PWM and sequence 
### information. 'min.score' setting for restriction. The higher the number, the stricter it will be
### strand = "*" means check the positive and negative strand
j <- 0
for (i in rownames(dnaTab)){
  j <- j+1
  seqName <- i
  seqFas <- as.character(dnaTab[seqName,])
  subject <- DNAString(seqFas)
  siteSetList <- searchSeq(PWMatrixList, subject, seqname = seqName, min.score = '95%', strand = '*')
  if (j == 1){
    write.table(writeGFF3(siteSetList, scoreType = 'relative'), file = 'result_vertebrates_2018_delta_DM_95.gff3', quote = F,
                sep = '\t', row.names = F)
  } else {
    write.table(writeGFF3(siteSetList, scoreType = 'relative'), file = 'result_vertebrates_2018_delta_DM_95.gff3', quote = F,
                sep = '\t', append = T, col.names = F, row.names = F)
  }
}


#================================================================================
#
#       Bubble plot for TFBS motif prediction (from HOMER)
#
#================================================================================
setwd('/Users/mijiarui/biosoft/HOMER/results')
bubble <- read.csv('bubble_TF.csv', header = T, quote = "", check.names = F, stringsAsFactors = T)
head(bubble)
library('ggplot2')
library('RColorBrewer')
bubble$Log_normalized_counts <- as.factor(round(bubble$log_normalized_counts))
ggplot(data = bubble, mapping = aes(x = Cell_type, y = Transcription_factor, size = Neg_logPvalue, col = log_normalized_counts)) + 
  geom_point() + theme_classic() + scale_color_gradient(low = 'blue', high = 'red')