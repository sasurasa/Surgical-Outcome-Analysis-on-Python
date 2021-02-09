##Load library
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(magrittr)
library(vsn)
library(ggplot2)
library(NMF)
library(grid)
library(gridExtra)
library(ggrepel)
library(dplyr)
library(plotly)
library(grDevices)
library(AnnotationDbi)
library(hgu133plus2.db)
library(hgu95av2.db)
library(biomaRt)
library(httr)
library(randomcoloR)
###############################################
##Clear all data
rm(list = ls())
###############################################
##Load data
list.files <- list.files(path = ".")
##Loop create data frame
i <- 1
for (j in list.files ) {
  if (i == 1) {
    file <- fread(j,header = TRUE)
    ## Add first sample
    read.counts <- data.frame(file$FPKM,row.names = file$gene_id)
    colnames(read.counts) <- "sample1"
    ## Set condition
    sample.info <- data.frame("sample1","normal")
    colnames(sample.info) <- c("name","condition")
  } else {
    ## Add another sample ##
    file <- fread(j, header = TRUE)
    ## Set sample name ##
    sample <- paste("sample",sep = "",i)
    ## Set condition
    if (i >= 2 & i <= 6) {
      condition <- "normal"
    } else {
      condition <- "cancer"
    }
    ## Insert sample to sample.info and read.counts ##
    count <- data.frame(file$FPKM)
    colnames(count) <- sample
    info <- data.frame(sample,condition)
    colnames(info) <- c("name","condition")
    read.counts <- cbind(read.counts,count)
    sample.info <- rbind(sample.info,info)
  }
  i <- i + 1
}
## Delete decimal number
read.counts[,-90] <- round(read.counts[,-90])

## Insert data to DESeq2 ##
DESeq.ds <- DESeqDataSetFromMatrix(read.counts ,colData = sample.info,~ condition)
colData(DESeq.ds)$condition  <- relevel(colData(DESeq.ds)$condition , "normal")

## Include gene that counts more than 0 ## 
DESeq.ds <- DESeq.ds[ rowSums(counts(DESeq.ds)) > 0, ]

## Calculate size factor of normalization counts ##
DESeq.ds <- estimateSizeFactors(DESeq.ds)

##########################################################
## NORMALIZATION ##
## Normalization and Log2 transform by DESeq2 ##
DESeq.rlog  <- vst(DESeq.ds , blind = TRUE)

## Summarized object data to matrix ##
rlog.norm.counts  <- assay(DESeq.rlog)

## Save dataframe of rlog.norm.counts with gene id
norm.counts <- data.frame(rlog.norm.counts)

## ggplot2 box plot ##
df_norm.counts <- list(counts = as.numeric(rlog.norm.counts), group = rownames(sample.info))
df_norm.counts <- data.frame(df_norm.counts)

# Plot
p <-  ggplot(df_norm.counts, aes(x = group, y = counts)) +
  geom_boxplot(varwidth = FALSE, outlier.colour = "black", outlier.size = 1.5, outlier.stroke = 1,
               outlier.shape = 1, notch = FALSE, fill = randomColor(14, luminosity="dark"),
               color="black", outlier.alpha = 1) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme_bw() + ggtitle("Normalized read counts") + ylab("log 2 read counts") + xlab("Sample")
p

## Mean SD Plot ##
#plot
msd_plot  <- meanSdPlot(rlog.norm.counts ,ranks=FALSE, plot = FALSE)
msd_plot$gg + ggtitle("log 2 read counts") +ylab("standard deviation")

## Crearte distance for each sample by pearson correlation ##
distance.m_rlog  <- as.dist(1 - cor(rlog.norm.counts , method = "pearson" ))

## Make dendogram ##
plot(hclust(distance.m_rlog),labels = colnames(rlog.norm.counts),
     main = "rlog  transformed  read  counts\ndistance: Pearson  correlation")

## Create PCA by ggplot2 ##
df_pca <- prcomp(t(rlog.norm.counts)) #Transform data row->col
df_out <- as.data.frame(df_pca$x) #Transform to data frame
df_out$group <- as.character(sample.info$condition) #Add group
p <- ggplot(df_out,aes(x=PC1,y=PC2, color=group, label=row.names(df_out) ))
p <- p + geom_point(size = 2) + geom_text(size=3,  vjust = -1.2, nudge_y = -0.2)
p
#####################################################
## DIFFERENTIATION ##
## Summarized object data to matrix ##
rlog.norm.counts  <- assay(DESeq.rlog)

## Compactly Display Structure ##
str(colData(DESeq.ds)$condition)

## Set level to differential expression ##
colData(DESeq.ds)$condition  <- relevel(colData(DESeq.ds)$condition , "normal")

## Set analysis ##
DESeq.ds <- DESeq(DESeq.ds) 

## Calculate negative binomial distribution ##
DESeq.ds <- estimateDispersions(DESeq.ds) ##Expression differentiation

## Calculate Wald Test for geralized linea model ##
DESeq.ds <- nbinomWaldTest(DESeq.ds) ##Weight distance of difference from GLM

## Extract result ##
DGE.results  <- results(DESeq.ds , independentFiltering = TRUE , alpha = 0.05)

## Histogram of P-value ##
hist(DGE.results$pvalue ,
     col = "grey", border = "white", xlab = "", ylab = "",
     main = "Frequencies  of P-value")

## MA plot ##
plot.MA <- data.frame(DGE.results)
plot.MA <- plot.MA %>% mutate(significant=padj<0.05) #Create significant column
MA <- ggplot(data = plot.MA, aes(x = baseMean, y = log2FoldChange, col=significant) )+
  geom_point(aes(color=as.factor(significant)), alpha=0.8, size=1.3) + scale_x_log10() + ggtitle("MA Plot of mRNA Expression") +
  theme_bw() 
MA

## Volcano plot ##
plot.volcano <- data.frame(DGE.results)
plot.volcano <- plot.volcano %>% mutate(significant=padj < 0.05)
row.names(plot.volcano) <- row.names(DGE.results)
volcano <- ggplot(plot.volcano) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=significant)) +
  ggtitle("Volcano Plot of mRNA Expression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  xlim(-3, 3) +
  theme_bw()
volcano

## Sort P-value ##
DGE.results.sorted  <- DGE.results[order(DGE.results$padj), ]

## Select first 30 gene ##
DGEgenes  <- rownames(subset(DGE.results.sorted[1:30,] ,))

## Annotate gene name ##
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = DGEgenes,
                  mart = ensembl)

rlog.norm.counts <- data.frame(rlog.norm.counts)

## Prepare data to txt ##
DGE.sort <- data.frame(DGE.results)
## Merge gene to DGE.sort ##
DGE.sort <- DGE.sort[order(DGE.sort$log2FoldChange), ]

## VISULIZATION HEATMAP ##
## Select normalized counts by gene of interest 
hm.mat_DGEgenes  <- rlog.norm.counts[DGEgenes , ]

## Annotate gene name ##
rownames(hm.mat_DGEgenes) <- genemap$hgnc_symbol

## Create heatmap ##
aheatmap(hm.mat_DGEgenes ,
         Rowv = TRUE , Colv = TRUE ,
         distfun = "pearson", hclustfun = "complete", 
         scale = "row")

###########################################################
## VISUALIZATION TOP N GENES##
DGEgenes  <- rownames(subset(DGE.results.sorted[1:10,] ,))
for (i in DGEgenes) {
  singlegene <- subset(genemap, genemap$ensembl_gene_id == i)
  title <- paste("Expression of ",singlegene$hgnc_symbol,"\n(",singlegene$ensembl_gene_id,")", sep="")
  ## Retrive normalized count ##
  plot.gene <- data.frame(t(subset(rlog.norm.counts, row.names(rlog.norm.counts)==i)))
  colnames(plot.gene) <- "count"
  ## Add group ##
  plot.gene$group <- sample.info$condition
  # Build file name
  file.name <- paste(singlegene$hgnc_symbol,"png",sep=".")
  ## Plot by ggplot2##
  p <- ggplot(plot.gene, aes(x = group, y = count)) +
    geom_boxplot(varwidth = T, notch = FALSE, fill=c("lightblue"), 
                 color="black", width = 0.5) +
    geom_jitter(shape=16, position=position_jitter(0.1), size = 2) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_bw() + ggtitle(title) + ylab("log 2 read counts") + xlab("Condition") +
    ggsave(file.name)
}

fwrite(read.counts,"read.count.txt")
fwrite(sample.info,"sample.info.txt")
fwrite(rlog.norm.counts,"rlog.norm.counts.txt")
fwrite(DGE.sort,"DEG.sort.txt")
