#This set of R-code is prepared for 'Genome Data Analysis' course, Biomedical Science Program, Faculty of Medicine, Prince of Songkla University

BiocManager::install('RNAseqData.HNRNPC.bam.chr14')

#Begin from BAM
library(RNAseqData.HNRNPC.bam.chr14)
alnfiles <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
alnfiles

#Defining experiment conditions
exp.des <- factor(rep(c('Control', 'HNRNPC knockdown'), each = 4))
exp.des
names(exp.des) <- RNAseqData.HNRNPC.bam.chr14_RUNNAMES
exp.des

#Alignment files to read-count table
BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)
seqlevels(txdb) <- c('chr14')
seqlevels(txdb)

#Applying Gene, Transcript and Exon 
genelist <- genes(txdb)
txlist <- transcripts(txdb)
exonlist <- exons(txdb)
genelist
txlist
exonlist

BiocManager::install('GenomicAlignments')
library(GenomicAlignments)
counts <- data.frame(row.names=genelist$gene_id)
for (i in 1:length(alnfiles)) {
  print(paste('Reading:',alnfiles[i]))
  aln <- readGAlignments(alnfiles[i])
  reads <- GRanges(seqnames = 'chr14', ranges=ranges(aln), strand='*')
  genecounts <- countOverlaps(genelist, reads)
  a <- data.frame(genecounts)
  names(a) <- RNAseqData.HNRNPC.bam.chr14_RUNNAMES[i]
  counts <- cbind(counts, a)
}
head(counts)
dim(counts)

#Removing low count (<10)
sum(rowSums(counts)==0)
sum(rowSums(counts)<10)
counts <- counts[!rowSums(counts)<10,]
dim(counts)

#Annotating Gene Symbol
BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
sym <- select(org.Hs.eg.db, row.names(counts), keytype = 'ENTREZID', 'SYMBOL')
#See what has been selected
head(sym)
dim(sym)
row.names(sym) <- sym$ENTREZID
head(sym)
sym$SYMBOL[is.na(sym$SYMBOL)] <- sym$ENTREZID[is.na(sym$SYMBOL)]

counts <- cbind(sym$SYMBOL, counts)
row.names(counts) <- as.character(counts$`sym$SYMBOL`)
counts$`sym$SYMBOL` <- NULL
head(counts)

#Creating a csv formatted table containing 'counts'
write.table(counts, file = 'HNRNPC_chr14_gene_read_counts.txt')

#Differential expression calculation
counts <- read.table('HNRNPC_chr14_gene_read_counts.txt')
head(counts)

BiocManager::install('edgeR')
#Single factor analysis with edgeR
library(edgeR)
genexp <- DGEList(counts=counts, group=exp.des)
genexp
## Notice that the libsize is the count for all genes, same number as 'read'
# Calculate normalization factor

genexp <- calcNormFactors(genexp)
genexp

genexp <- estimateCommonDisp(genexp)
genexp


genexp <- estimateTagwiseDisp(genexp)
genexp

#Notice that CommonDisp has only one value when TagwiseDisp has dimension equal to row number
genediff <- exactTest(genexp)
genediff
topTags(genediff)

#Constructing a table showing differential expression with edgeR, trimming out those p-value >= 0.03

genediff$table <- cbind(genediff$table, FDR=p.adjust(genediff$table$PValue, method = 'fdr'))
gsign <- genediff$table[genediff$table$FDR<0.03,]
gsign <- gsign[order(gsign$FDR),]
dim(gsign)
head(gsign)

head(gsign, 10)
##logFC = log(fold change), logCPM = log(count per million of that gene)

#Differential expression with DESeq (Anders and Huber 2010)

BiocManager::install('DESeq')
library(DESeq)
count.set <- newCountDataSet(counts, exp.des)
count.set
##Notice that the row numbers are equal to 'count'

#Normalization
count.set <- estimateSizeFactors(count.set)
sizeFactors(count.set)

##See the correlation between the two normalization factors
plot(genexp$samples$norm.factors, sizeFactors(count.set), xlab='edgeR', ylab='DESeq',
     main = 'Normalization factors', type='n', xlim=c(0.5,1.5), ylim=c(0.5,1.5)) 
abline(0,1,lty=2)
text(genexp$samples$norm.factors, sizeFactors(count.set), labels = names(sizeFactors(count.set)), cex = 0.8) 

#See the difference in normalization between the 2 methods
count.set <- estimateDispersions(count.set)
par(mfrow=c(1,2))
plotDispEsts(count.set, main = 'DESeq: genewide dispersion', cex=1)
plotBCV(genexp, main='edgeR: genewide dispersion normalized counts', cex=1)
par(mfrow=c(1,1))

head(cbind(genexp$tagwise.dispersion, fData(count.set)))
plot(cbind(genexp$tagwise.dispersion, fData(count.set)), ylim=c(0,2), xlim=c(0,2), 
     main='Gene-wide dispersion', xlab='edgeR', ylab= 'DESeq')
abline(0,1,lty=2)

#Differential expression MA plot
genediff.ds <- nbinomTest(count.set, 'Control', 'HNRNPC knockdown')
plot(log2FoldChange ~ baseMean, data = genediff.ds[genediff.ds$padj>0.03,], log='x', 
     xlab='Mean of normalized counts', ylab='log fold change', main ='MA plot', pch=20)
abline(h=0, lwd=2, col='red')
points(genediff.ds[genediff.ds$padj<0.03, c('baseMean', 'log2FoldChange')], col='red') 

gsign.ds <- genediff.ds[genediff.ds$padj < 0.03,]
gsign.ds <- gsign.ds[order(gsign.ds$padj),]
row.names(gsign.ds) <- gsign.ds$id

head(gsign.ds, 10)
head(gsign, 10)

#Comparing edgeR and DESeq
par(mfrow=c(1,3))
plot(gsign[intersect(row.names(gsign), gsign.ds$id), 'logFC'],gsign.ds[intersect(row.names(gsign),gsign.ds$id), 'log2FoldChange'], 
     main='Log fold-changes', xlab='edgeR', ylab='DESeq')
plot(gsign[intersect(row.names(gsign), gsign.ds$id), 'FDR'],gsign.ds[intersect(row.names(gsign),gsign.ds$id), 'padj'], 
     main='Log fold-changes', xlab='edgeR', ylab='DESeq')
plot(gsign[intersect(row.names(gsign), gsign.ds$id), 'FDR'],gsign.ds[intersect(row.names(gsign),gsign.ds$id), 'padj'], 
     main='Log fold-changes', xlab='edgeR', ylab='DESeq', log='xy', xlim=c(1e-20, 1e-1), ylim=c(1e-20, 1e-1))
par(mfrow=c(1,1))

#Using general linear model with edgeR and DESeq
repls <- rep(rep(1:2, each = 2), 2)
repls
design <- model.matrix(~exp.des+repls)
design

#GLM with edgeR

genexp.glm <- DGEList (counts = counts)
genexp.glm
genexp.glm <- estimateGLMCommonDisp(genexp.glm, design)
genexp.glm <- estimateGLMTagwiseDisp(genexp.glm, design)

gene.exp.fit <- glmFit(genexp.glm, design)
gene.exp.fit

#Comparing KO vs Control suing logistic regression
genediff.lrt <- glmLRT(gene.exp.fit, coef=2)
topTags(genediff.lrt)
#Comparing Second design using logistic regression, extended to include repls
genediff.lrt.2 <- glmLRT(gene.exp.fit, coef=3)
topTags(genediff.lrt.2)

genediff.lrt
genediff.lrt$table <- cbind(genediff.lrt$table, FDR=p.adjust(genediff.lrt$table$PValue, method='fdr'))
genediff.lrt$table
gsign.lrt <- genediff.lrt$table[genediff.lrt$table$FDR<0.03,]
gsign.lrt <- gsign.lrt[order(gsign.lrt$FDR),]
dim(gsign.lrt)
dim(gsign)
head(gsign.lrt, 10)
head(gsign, 10)

#GLM with DESeq
exp.des.cmplx <- data.frame(Type=as.character(exp.des), Repl=repls)
row.names(exp.des.cmplx) <- names(exp.des)
exp.des.cmplx

count.set.glm <- newCountDataSet(counts, exp.des.cmplx)
count.set.glm <- estimateSizeFactors(count.set.glm)
count.set.glm <- estimateDispersions(count.set.glm)

#Fitting the model using Chi-square
genediff.ds.mod1 <- fitNbinomGLMs(count.set.glm, count ~ Type+Repl)
genediff.ds.mod2 <- fitNbinomGLMs(count.set.glm, count ~ Type)
genediff.ds.glm.pvals <- nbinomGLMTest(genediff.ds.mod1, genediff.ds.mod2)
genediff.ds.glm.fdr <- p.adjust(genediff.ds.glm.pvals, method='fdr')
genediff.ds.glm <- cbind(genediff.ds.mod1[,2:3], Pval=genediff.ds.glm.pvals, FDR=genediff.ds.glm.fdr)

gsign.ds.glm <- genediff.ds.glm[genediff.ds.glm$FDR<0.03,]
gsign.ds.glm <- genediff.ds.glm[order(gsign.ds.glm),]
head(gsign.ds.glm)



dim(gsign.lrt)
head(gsign.lrt)
dim(gsign.ds.glm)
head(gsign.ds.glm)
common.genes <- intersect(row.names(gsign.lrt), row.names(gsign.ds.glm))
dim(common.genes)
plot(gsign.lrt[common.genes, 'locFC'], gsign.ds.glm[common.genes, 'Type HNRNPC knockdown'], 
     main='Fold change in GLMs', xlab='edgeR', ylab='DESeq')
abline (0,1,lty=2)

#Almost there: Generating a heatmap
library(gplots)



head(gsign)
head(gsign.ds)
diff.exp <- data.frame(edgeR=gsign[intersect(row.names(gsign), gsign.ds$id), 'logFC'],
                       DESeq=gsign.ds[intersect(row.names(gsign), gsign.ds$id), 'log2FoldChange'])
diff.exp
row.names(diff.exp) <- intersect(row.names(gsign), gsign.ds$id)
diff.exp<- diff.exp[!is.infinite(diff.exp$DESeq),]
diff.exp

my.palette <- rev(redgreen(75))
heatmap.2(log2(as.matrix(counts[rownames(gsign),]+1)),
          col=my.palette, margins=c(10,12))
heatmap.2(as.matrix(diff.exp),dendrogram='row', Colv = FALSE,
          col=my.palette, margins=c(12,10))




