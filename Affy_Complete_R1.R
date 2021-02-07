#R program for this session

BiocManager::install("affy")

getwd()

datapath <- ('~/Downloads/GSE29797_RAW')
setwd(datapath)
library(affy)
dat <- ReadAffy()

dat
annotation(dat)
#Chip definition format file
cdfName(dat)

#Putting-in phenotype data

pData(dat)

exp.des <- pData(dat)
#Adding 'Genotype' column with a vector containing 'WT and 'dTsc'.
exp.des$Genotype <- factor(rep(c('WT','dTsc1'),each=8))
#Adding 'Stimulation' column
exp.des$Stimulation <- factor(rep(c('0h','4h'),8))
exp.des
#Return the df 
pData(dat) <- exp.des

image(dat[,1])
image(dat[,1:4])

hist(dat)
hist(dat[,1])

boxplot(exprs(dat))
boxplot(log2(exprs(dat)))

#QC
#RNA Degradation

deg <- AffyRNAdeg(dat)
summaryAffyRNAdeg(deg)
plotAffyRNAdeg(deg)

BiocManager::install('affyQCReport')

library(affyQCReport)
correlationPlot(dat)

#Normalization by Robust multi-array average (RMA)

datrma <- rma(dat)
datrma

#Take expression value out of normalized file
matexp <- exprs(datrma)

par(mfrow=c(1,2))
boxplot(log2(exprs(dat)), main='Before normalization', ylab='log2(intensity)')
boxplot(matexp, main='RMA normalized data')
par(mfrow=c(1,1))

#Plot PCA function


library(affycoretools)

#Find principal component in expression (matexp) that best group 16 samples based on their Genotype or Stimulation
par(mfrow=c(1,2))
plotPCA(matexp, groups = as.numeric(pData(dat)[,2]), groupnames = levels(pData(dat)[,2]))
plotPCA(matexp, groups = as.numeric(pData(dat)[,3]), groupnames = levels(pData(dat)[,3]))
par(mfrow-c(1,1))

#Set experiment matrix
genotype <- factor(pData(dat)[,2], levels = levels(pData(dat)[,2]))
stimuli <- factor(pData(dat)[,3], levels = levels(pData(dat)[,3]))
design <- model.matrix(~genotype)
design.1 <- model.matrix(~stimuli)

BiocManager::install('limma')
library(limma)

#Fit and compare the expression matrix by Genotype, report the top 50 
fit <- lmFit(matexp, design)
fit

fit <- eBayes(fit)
dg.top.50 <- topTable(fit, coef = 2, adjust = 'fdr', n=50)
head(dg.top.50)

#Fit and compare the expression matrix by Stimulation, report the top 50
fit.1 <- lmFit(matexp, design.1)
fit.1 <- eBayes(fit.1)
dg.top.50.stim <- topTable(fit.1, coef = 2, adjust = 'fdr', n=50)
head(dg.top.50.stim)

#Select only fit with p.value < 0.05
#Condition
selected <- fit$p.value[,2] 
selected <- selected < 0.05
selected
#Selection
matexp.gen <- matexp[selected,]
#See new dimension
dim(matexp.gen)
heatmap(matexp.gen)

#Annotate gene names
mouse.annot <-read.delim('~/Downloads/GPL11180.annot', skip = 27, row.names = 1)
head(mouse.annot)
sel.annot <- mouse.annot[row.names(dg.top.50),]
dg.top.50.annotated <- cbind(dg.top.50, sel.annot)
write.table(dg.top.50.annotated[order(dg.top.50.annotated$logFC), c(1,5,7,17:15)],'diff_genes_genotypes.txt')
sel.annot <- mouse.annot[row.names(dg.top.50.stim),]
dg.top.50.stim.annotated <- cbind(dg.top.50.stim, sel.annot)
write.table(dg.top.50.stim.annotated[order(dg.top.50.stim.annotated$logFC), c(1,5,7,17:15)],'diff_genes_stimulation.txt')

dg.top.50.annotated$Gene.symbol
dg.top.50.stim.annotated$Gene.symbol
