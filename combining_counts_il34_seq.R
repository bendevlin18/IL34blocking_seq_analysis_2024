



rm(list = ls())
library('biomaRt')
library('EDASeq')
library('tidyverse')


## post that outlined this workflow
## https://stackoverflow.com/questions/28543517/how-can-i-convert-ensembl-id-to-gene-symbol-in-r


## reading in the counts file
###

## list all the files in given directory
count_direc <- 'C:\\Users\\Ben\\Desktop\\il34_RNASeq_fastq\\counts'
files <- list.files(count_direc)

## load in first file for a template df
df <- read.table(paste(count_direc, files[1], sep = '\\'), fill = TRUE, header = 1)
df <- select(df, c('Length', 'Geneid'))


## read in all the files in the directory and adding them to the template df
for (file in files) {
  if (grepl('.txt', file, fixed=TRUE)) {
    df[file] <- read.table(paste(count_direc, file, sep = '\\'), fill = TRUE, header = 1)[7]
    print(paste('Reading..', file))
  }
}
  
## creating gene list with hgnc symbols from BiomaRt
## grabbing mus musculus gene ensembl ID matches on BiomaRt
mouse <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","gene_biotype", 'external_gene_name'), values = df$Geneid, mart = mouse)

## merging gene names with count table using ensembl IDs
merged <- merge(df, gene_list, by.x="Geneid", by.y="ensembl_gene_id")


write.csv(merged, 'C:\\Users\\Ben\Desktop\\il34_RNASeq_fastq\\combined_counts_RAW.csv')

protein_coding_genes <- merged[merged$gene_biotype == 'protein_coding', ]

write.csv(merged, 'C:\\Users\\Ben\\Desktop\\il34_RNASeq_fastq\\combined_counts_protein_coding.csv')

