#######@FPKM count in linux#############
#--------------------------------------#
#######make by Kogoori Masaki###########
#--------------------------------------#
##########r-base version > 4.2.2########
##Usage:Rscript pipeline.r input.csv experiment_rownumber control_rownumber total_number experiment_name control_name title_name
library(tidyverse,quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(ggplot2)
library(pheatmap)
library(FactoMineR)
argv<-commandArgs(TRUE)
a <- argv[1]
b <- argv[2]
c <- argv[3]
d <- argv[4]
name1 <- argv[5]
name2 <- argv[6]
titlename <- argv[7]
file <- read.csv(a,header = TRUE,row.names = 1)
count <- file[1:d]
count <- as.matrix(count)
sample <- as.data.frame(colnames(file))
sample$deal <- c(rep(name1,b),rep(name2,c))
names(sample) <- c("sam","condition")
rownames(sample) <- sample$sam
sample %>% select(-1) -> sample
dds <- DESeqDataSetFromMatrix(countData = count,colData = sample,design = ~ condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', name1, name2))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, "diffselect.txt", col.names = NA, sep = '\t', quote = FALSE)
res1<-res1[complete.cases(res1),]
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'updownselect.txt', sep = '\t', col.names = NA, quote = FALSE)

p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  
  scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down')) + 
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = titlename, color = '') +  
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') + 
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-15, 15) + ylim(0, 500)  
ggsave(p,filename = "vocalno.pdf",width = 12,height = 9)

vsd <- varianceStabilizingTransformation(dds)
data <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
q <- ggplot(data, aes(PC1, PC2, color = condition)) + 
  geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(title = titlename ) + theme_bw()
ggsave(q,filename = "pca.pdf",width = 12,height = 9)

ddCor <- cor(count, method = "spearman")
pheatmap(file="pheatmap.pdf",ddCor, clustering_method = "average", 
         display_numbers = F,legend=TRUE,show_rownames=F,show_colnames=F,border = F,
         treeheight_row=50)
