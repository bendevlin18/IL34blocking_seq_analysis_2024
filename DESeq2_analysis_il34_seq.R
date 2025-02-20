
rm(list = ls())
library(DESeq2)
library(EnhancedVolcano)
library(fdrtool)
library(tidyverse)

count_df <- read.csv('C:\\Users\\Ben\\Box\\duke_bilbo\\MemberFolders\\Ben\\il34_project\\il34_bulk_RNASeq\\data_analysis\\combined_counts_protein_coding.csv', fileEncoding="UTF-8-BOM")
mgla_metadata <- read.csv('C:\\Users\\Ben\\Box\\duke_bilbo\\MemberFolders\\Ben\\il34_project\\il34_bulk_RNASeq\\data_analysis\\il34_seq_full_sample_metadata_nos4_20_29a.csv', fileEncoding="UTF-8-BOM")
setwd('C:\\Users\\Ben\\Box\\duke_bilbo\\MemberFolders\\Ben\\il34_project\\il34_bulk_RNASeq\\data_analysis')
## cleaning the count df to only include protein coding genes
count_df <- count_df[count_df$gene_biotype == 'protein_coding',]
count_df <- count_df[!duplicated(count_df$external_gene_name), ]
row.names(count_df) <- count_df$external_gene_name
count_df <- count_df[-which(grepl("^mt-",row.names(count_df))), ]


## transposing count df so that we can merge it with metadata
t_count_df <- data.frame(t(count_df))
t_count_df$id <- row.names(t_count_df)

####################################################################################
####################################################################################
####################################################################################
####################################################################################
################ Creating TPM with whole matrix ###############
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')
metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$SampleID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)
raw_counts <- counts(dds)

z <- data.frame(t(raw_counts))
z$sample_ID <- row.names(z)
merged_raw_counts <- merge(z, metadata, by.x = 'sample_ID', by.y = 'SampleID')
t_merged_raw_counts <- setNames(data.frame(t(merged_raw_counts[,-1])), merged_raw_counts[,1])
t_merged_raw_counts$gene_name <- row.names(t_merged_raw_counts)
final_merge_raw_counts <- merge(t_merged_raw_counts, count_df[c('external_gene_name', 'Length')], by.x = 'gene_name', by.y = 'external_gene_name')


#### making TPM ####

## batching through the column names to calculate TPM for each sample
for (cols in colnames(final_merge_raw_counts)) {
  
  if (grepl('X', cols, fixed=TRUE)) {
    print(cols)
    #Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
    lengths <- final_merge_raw_counts$Length
    rpk <- as.numeric(unlist(final_merge_raw_counts[cols])) / lengths
    #Count up all the RPK values in a sample and divide this number by 1,000,000. This is your per million scaling factor.
    per_mil <- sum(rpk, na.rm=T) / 1000000
    #Divide the RPK values by the per million scaling factor. This gives you TPM.
    tpm <- rpk/per_mil
    final_merge_raw_counts[paste(cols, 'tpm', sep = '_')] <- tpm
  }
}

write.csv(final_merge_raw_counts, 'tpm.csv')


################ RUN FROM HERE ###############################
################ FOR KO EXPERIMENT ###############################



## extracting sampleID metadata to only include CD11b+ samples
merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')


######## FILTERS TO CHANGE ########
#merged <- merged[merged$sex == 'f', ]
merged <- merged[merged$TissueType == 'Mgla', ]
merged <- merged[merged$Expt == 'KO', ]
######## FILTERS TO CHANGE ########


metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$sample_ID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)

## keep genes that are expressed by at least 4 samples at 1 count or more
keep <- rowSums(counts(dds) >= 10) >= 4

dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Genotype","IL34 Het","IL34 KO"))
summary(res)
hist(res$pvalue)

sig_genes <- data.frame(subset(res, padj < 0.001))
sig_genes$gene_name <- rownames(sig_genes)

down_genes <- sig_genes[sig_genes$log2FoldChange < -1,]
up_genes <- sig_genes[sig_genes$log2FoldChange > 1,]

write.csv(rbind(down_genes, up_genes), 'KO_mgla_sig_genes.csv')

down_genes <- down_genes %>%
  arrange(down_genes$padj) %>%
  slice(1:20)
up_genes <- up_genes %>%
  arrange(up_genes$padj) %>%
  slice(1:20) 

EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = c('Ifitm3', 'Irf7', 'Ifit3', 'Axl', 'Aif1'),
                #selectLab = c(up_genes$gene_name, down_genes$gene_name),
                x = 'log2FoldChange',
                y = 'padj',
                maxoverlapsConnectors = 150,
                pCutoff = .001,
                FCcutoff = 1,
                labSize = 5,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                colConnectors = 'black',
                legendPosition = 'None',
                title = 'Up in KO vs. Up in Het (Ctrl)'
)




################ RUN FROM HERE ###############################
################ FOR BLOCKING AB EXPT ###############################



## extracting sampleID metadata to only include CD11b+ samples
merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')


######## FILTERS TO CHANGE ########
#merged <- merged[merged$sex == 'f', ]
merged <- merged[merged$TissueType == 'Mgla', ]
merged <- merged[merged$Expt == 'ab', ]
######## FILTERS TO CHANGE ########


metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$sample_ID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)

## keep genes that are expressed by at least 4 samples at 1 count or more
keep <- rowSums(counts(dds) >= 10) >= 4

dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Genotype", "a-gp120", "a-IL34"))
summary(res)
hist(res$pvalue)

sig_genes <- data.frame(subset(res, padj < 0.001))
sig_genes$gene_name <- rownames(sig_genes)

down_genes <- sig_genes[sig_genes$log2FoldChange < -1,]
up_genes <- sig_genes[sig_genes$log2FoldChange > 1,]

write.csv(rbind(down_genes, up_genes), 'blocking_ab_mgla_sig_genes.csv')


down_genes <- down_genes %>%
  arrange(down_genes$padj) %>%
  dplyr::slice(1:10)


up_genes <- up_genes %>%
  arrange(up_genes$padj) %>%
  dplyr::slice(1:10) 


EnhancedVolcano(res,
                lab = rownames(res),
                #selectLab = c('Ifitm3', 'Tmem119', 'P2ry12', 'Axl', 'Mx1', 'Ifit2', 'Sirpa', 'Siglech', 'Cx3cr1', 'Olfml3', 'Hexb', 'Snap25', 'Aqp4', 'Olig2', 'Sox9', 'Aif1'),
                #selectLab = c('Ifitm3', 'Aif1', 'Hexb', 'Cx3cr1', 'Siglech', 'Tmem173', 'Csf1r', 'C3ar1', 'Slc2a5', 'Il6ra', 'Cd44', 'Serping1', 'Emp1', 'Mog', 'Gad1', 'Mbp', 'P2ry6', 'Mx1', 'Irf7'),
                selectLab = c(up_genes$gene_name, down_genes$gene_name, 'Tmem119', 'Mx1'),
                x = 'log2FoldChange',
                y = 'padj',
                maxoverlapsConnectors = 150,
                pCutoff = .001,
                FCcutoff = 1,
                labSize = 5,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                colConnectors = 'black',
                legendPosition = 'None',
                title = 'Up in a-IL34 vs. Up in Ctrl'
)


ggsave('blocking_ab_isolated_microglia_volcano.svg', EnhancedVolcano(res,
                                                                     lab = rownames(res),
                                                                     #selectLab = c('Ifitm3', 'Tmem119', 'P2ry12', 'Axl', 'Mx1', 'Ifit2', 'Sirpa', 'Siglech', 'Cx3cr1', 'Olfml3', 'Hexb', 'Snap25', 'Aqp4', 'Olig2', 'Sox9', 'Aif1'),
                                                                     #selectLab = c('Ifitm3', 'Aif1', 'Hexb', 'Cx3cr1', 'Siglech', 'Tmem173', 'Csf1r', 'C3ar1', 'Slc2a5', 'Il6ra', 'Cd44', 'Serping1', 'Emp1', 'Mog', 'Gad1', 'Mbp', 'P2ry6', 'Mx1', 'Irf7'),
                                                                     selectLab = c(up_genes$gene_name, down_genes$gene_name, 'Tmem119', 'Mx1'),
                                                                     x = 'log2FoldChange',
                                                                     y = 'padj',
                                                                     maxoverlapsConnectors = 150,
                                                                     pCutoff = .001,
                                                                     FCcutoff = 1,
                                                                     labSize = 5,
                                                                     col=c('black', 'black', 'black', 'red3'),
                                                                     colAlpha = 1,,
                                                                     labCol = 'black',
                                                                     labFace = 'bold',
                                                                     drawConnectors = TRUE,
                                                                     colConnectors = 'black',
                                                                     legendPosition = 'None',
                                                                     title = 'Up in a-IL34 vs. Up in Ctrl'
))



################ RUN FROM HERE ###############################
################ FOR KO EXPERIMENT ###############################



## extracting sampleID metadata to only include CD11b+ samples
merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')


######## FILTERS TO CHANGE ########
#merged <- merged[merged$sex == 'f', ]
merged <- merged[merged$TissueType == 'Brain', ]
merged <- merged[merged$Expt == 'KO', ]
######## FILTERS TO CHANGE ########


metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$sample_ID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)

## keep genes that are expressed by at least 4 samples at 1 count or more
keep <- rowSums(counts(dds) >= 10) >= 4

dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Genotype","IL34 Het","IL34 KO"))
summary(res)
hist(res$pvalue)

sig_genes <- data.frame(subset(res, padj < 0.001))
sig_genes$gene_name <- rownames(sig_genes)

down_genes <- sig_genes[sig_genes$log2FoldChange < -1,]
down_genes <- down_genes %>%
  arrange(down_genes$padj) %>%
  slice(1:20)

up_genes <- sig_genes[sig_genes$log2FoldChange > 1,]
up_genes <- up_genes %>%
  arrange(up_genes$padj) %>%
  slice(1:20) 


write.csv(sig_genes, "brain_KO_genes.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                #selectLab = c( 'Tmem119', 'P2ry12', 'Il34', 'Aif1', 'Cx3cr1', "Csf1"),
                selectLab = c(up_genes$gene_name, down_genes$gene_name),
                x = 'log2FoldChange',
                y = 'padj',
                maxoverlapsConnectors = 150,
                pCutoff = .001,
                FCcutoff = 1,
                labSize = 5,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                colConnectors = 'black',
                legendPosition = 'None',
                title = 'Up in KO vs. Up in Het (Ctrl)'
)


ggsave('KO_brain_volcano.svg',EnhancedVolcano(res,
                                              lab = rownames(res),
                                              #selectLab = c( 'Tmem119', 'P2ry12', 'Il34', 'Aif1', 'Cx3cr1', "Csf1"),
                                              selectLab = c(up_genes$gene_name, down_genes$gene_name),
                                              x = 'log2FoldChange',
                                              y = 'padj',
                                              maxoverlapsConnectors = 150,
                                              pCutoff = .001,
                                              FCcutoff = 1,
                                              labSize = 5,
                                              col=c('black', 'black', 'black', 'red3'),
                                              colAlpha = 1,
                                              labCol = 'black',
                                              labFace = 'bold',
                                              drawConnectors = TRUE,
                                              colConnectors = 'black',
                                              legendPosition = 'None',
                                              title = 'Up in KO vs. Up in Het (Ctrl)'
)
)


################ RUN FROM HERE ###############################
################ FOR BLOCKING AB EXPT ###############################



## extracting sampleID metadata to only include CD11b+ samples
merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')


######## FILTERS TO CHANGE ########
#merged <- merged[merged$sex == 'f', ]
merged <- merged[merged$TissueType == 'Brain', ]
merged <- merged[merged$Expt == 'ab', ]
######## FILTERS TO CHANGE ########


metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$sample_ID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)

## keep genes that are expressed by at least 4 samples at 1 count or more
keep <- rowSums(counts(dds) >= 10) >= 4

dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Genotype", "a-gp120", "a-IL34"))
summary(res)
hist(res$pvalue)

sig_genes <- data.frame(subset(res, padj < 0.01))
sig_genes$gene_name <- rownames(sig_genes)

down_genes <- sig_genes[sig_genes$log2FoldChange < -1,]
up_genes <- sig_genes[sig_genes$log2FoldChange > 1,]

write.csv(rbind(down_genes, up_genes), 'blocking_ab_brain_sig_genes.csv')

down_genes <- down_genes %>%
  arrange(down_genes$padj) %>%
  slice(1:10)


up_genes <- up_genes %>%
  arrange(up_genes$padj) %>%
  slice(1:10) 


EnhancedVolcano(res,
                lab = rownames(res),
                #selectLab = c( 'Aif1', 'Hexb'),
                selectLab = c(up_genes$gene_name, down_genes$gene_name),
                x = 'log2FoldChange',
                y = 'padj',
                maxoverlapsConnectors = 150,
                pCutoff = .001,
                FCcutoff = 1,
                labSize = 5,
                labCol = 'black',
                labFace = 'bold',
                colCustom = keyvals,
                colAlpha = 1,
                drawConnectors = TRUE,
                colConnectors = 'black',
                legendPosition = 'None',
                title = 'Up in a-IL34 vs. Up in Ctrl'
)



ggsave('blocking_ab_brain_volcano.svg', EnhancedVolcano(res,
                                                        lab = rownames(res),
                                                        #selectLab = c( 'Aif1', 'Hexb'),
                                                        selectLab = c(up_genes$gene_name, down_genes$gene_name),
                                                        x = 'log2FoldChange',
                                                        y = 'padj',
                                                        maxoverlapsConnectors = 150,
                                                        pCutoff = .001,
                                                        FCcutoff = 1,
                                                        labSize = 5,
                                                        labCol = 'black',
                                                        labFace = 'bold',
                                                        colCustom = keyvals,
                                                        colAlpha = 1,,
                                                        drawConnectors = TRUE,
                                                        colConnectors = 'black',
                                                        legendPosition = 'None',
                                                        title = 'Up in a-IL34 vs. Up in Ctrl'
))




################ RUN FROM HERE ###############################
################ COMPARING MICROGLIA ###############################



## extracting sampleID metadata to only include CD11b+ samples
merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')


######## FILTERS TO CHANGE ########
#merged <- merged[merged$sex == 'f', ]
merged <- merged[merged$TissueType == 'Mgla', ]
merged <- merged[merged$Genotype != 'IL34 Het', ]
merged <- merged[merged$Genotype != 'a-gp120', ]
######## FILTERS TO CHANGE ########


metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$sample_ID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)

## keep genes that are expressed by at least 4 samples at 1 count or more
keep <- rowSums(counts(dds) >= 10) >= 4

dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Genotype", "IL34 KO", "a-IL34"))
summary(res)
hist(res$pvalue)

sig_genes <- data.frame(subset(res, padj < 0.001))
sig_genes$gene_name <- rownames(sig_genes)


down_genes <- sig_genes[sig_genes$log2FoldChange < -1,]
down_genes <- down_genes %>%
  arrange(down_genes$padj) %>%
  slice(1:20)

up_genes <- sig_genes[sig_genes$log2FoldChange > 1,]
up_genes <- up_genes %>%
  arrange(up_genes$padj) %>%
  slice(1:20) 


EnhancedVolcano(res,
                lab = rownames(res),
                #selectLab = c( 'Tmem119', 'P2ry12', 'Il34', 'Aif1', 'Cx3cr1', 'Csf1', 'Ifitm3', 'Oasl2'),
                selectLab = c(up_genes$gene_name, down_genes$gene_name),
                x = 'log2FoldChange',
                y = 'padj',
                maxoverlapsConnectors = 150,
                pCutoff = .001,
                FCcutoff = 1,
                labSize = 5,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                colConnectors = 'black',
                legendPosition = 'None',
                title = 'Up in a-IL34 vs. Up in IL34 KO'
)








################ RUN FROM HERE ###############################
################ COMPARING MICROGLIA ###############################



## extracting sampleID metadata to only include CD11b+ samples
merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')


######## FILTERS TO CHANGE ########
#merged <- merged[merged$sex == 'f', ]
merged <- merged[merged$TissueType == 'Mgla', ]
merged <- merged[merged$Genotype != 'IL34 KO', ]
merged <- merged[merged$Genotype != 'a-IL34', ]
######## FILTERS TO CHANGE ########


metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$sample_ID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)

## keep genes that are expressed by at least 4 samples at 1 count or more
keep <- rowSums(counts(dds) >= 100) >= 4

dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Genotype", "IL34 Het", "a-gp120"))
summary(res)
hist(res$pvalue)

sig_genes <- data.frame(subset(res, padj < 0.01))
sig_genes$gene_name <- rownames(sig_genes)

write.csv(rbind(down_genes, up_genes), 'ctrl_mgla_showdown_sig_genes.csv')

down_genes <- sig_genes[sig_genes$log2FoldChange < -1,]
down_genes <- down_genes %>%
  arrange(down_genes$padj) %>%
  slice(1:20)

up_genes <- sig_genes[sig_genes$log2FoldChange > 1,]
up_genes <- up_genes %>%
  arrange(up_genes$padj) %>%
  slice(1:20) 



EnhancedVolcano(res,
                lab = rownames(res),
                #selectLab = c( 'Tmem119', 'P2ry12', 'Il34', 'Aif1', 'Cx3cr1', 'Csf1', 'Ifitm3', 'Oasl2'),
                selectLab = c(up_genes$gene_name, down_genes$gene_name),
                x = 'log2FoldChange',
                y = 'padj',
                maxoverlapsConnectors = 150,
                pCutoff = .01,
                FCcutoff = 1,
                labSize = 5,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                colConnectors = 'black',
                legendPosition = 'None',
                title = 'Up in a-gp120 vs. Up in IL34 Het'
)








################ RUN FROM HERE ###############################
################ COMPARING MICROGLIA to whole brain ###############################



## extracting sampleID metadata to only include CD11b+ samples
merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')


######## FILTERS TO CHANGE ########
#merged <- merged[merged$sex == 'f', ]
#merged <- merged[merged$TissueType == 'Mgla', ]
merged <- merged[merged$Genotype != 'IL34 Het', ]
merged <- merged[merged$Genotype != 'a-gp120', ]
######## FILTERS TO CHANGE ########


metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$sample_ID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ TissueType)

## keep genes that are expressed by at least 4 samples at 1 count or more
keep <- rowSums(counts(dds) >= 1) >= 20

dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("TissueType", "Brain", "Mgla"))
summary(res)
hist(res$pvalue)

sig_genes <- data.frame(subset(res, padj < 0.01))
sig_genes$gene_name <- rownames(sig_genes)


down_genes <- sig_genes[sig_genes$log2FoldChange < -1,]
down_genes <- down_genes %>%
  arrange(down_genes$padj) %>%
  slice(1:20)

up_genes <- sig_genes[sig_genes$log2FoldChange > 1,]
up_genes <- up_genes %>%
  arrange(up_genes$padj) %>%
  slice(1:20) 


EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = c( 'Tmem119', 'P2ry12', 'Il34', 'Aif1', 'Cx3cr1', 'Csf1', 'Ifitm3', 'Oasl2', 'Snap25', 'Olig2', 'Aqp4', 'Cldn11'),
                #selectLab = c(up_genes$gene_name, down_genes$gene_name),
                x = 'log2FoldChange',
                y = 'padj',
                maxoverlapsConnectors = 150,
                pCutoff = .01,
                FCcutoff = 1,
                labSize = 5,
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                colConnectors = 'black',
                legendPosition = 'None',
                title = 'Up in mgla vs. Up in Whole Brain'
)










################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')

#merged <- merged[merged$TissueType != 'Brain', ]
#merged <- merged[merged$Expt != 'ab', ]

metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$SampleID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)

plotPCA(vst(dds), intgroup = c('TissueType'))





################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################
################################### PCA PLOT ####################################



merged <- merge(mgla_metadata, t_count_df, by.x = 'SampleID', by.y = 'id')

merged <- merged[merged$TissueType != 'Brain', ]
merged <- merged[merged$Expt != 'KO', ]

metadata <- data.frame(merged[, 0:6])
t_merged <- t(merged[, -(1:6)])
count_mtx <- matrix(as.numeric(t_merged), ncol = ncol(t_merged))
count_mtx[is.na(count_mtx)] <- 0
row.names(count_mtx) <- row.names(t_merged)
colnames(count_mtx) <- merged$SampleID

dds <- DESeqDataSetFromMatrix(countData = count_mtx, colData = metadata, design = ~ Genotype)
keep <- rowSums(counts(dds) >= 100) >= 8
dds <- dds[keep,]

plotPCA(vst(dds), intgroup = c('Genotype'))+
  coord_fixed(ratio = 1)

ggsave('KO_mgla_PCA.svg', plotPCA(vst(dds), intgroup = c('Genotype')) + coord_fixed(ratio = 4), width = 4, height = 7)
