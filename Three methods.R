setwd("/Users/huachunyin/desktop/real data")

count<-read.csv("GSE124326.csv",header = T,row.names  = 1)

lable<-read.csv("GSE124326_lable.csv",header = T)

group_list<-factor(lable$lable)

pvalue = 0.05
foldChange = 1
#########################
library("edgeR")

###filtering#PMID: 23975260#####
noint = rownames(count) %in%
  c("no_feature","ambiguous","too_low_aQual",
    "not_aligned","alignment_not_unique")
cpms = cpm(count)
keep = rowSums(cpms>1)>=3 & !noint
count = count[keep,]
class(count)
head(count)
#######normalization

#######
#group_list<-factor(c(rep("disease",4),rep("treatment",4)))


design <- model.matrix(~0+group_list)

rownames(design) = colnames(count)

colnames(design) <- levels(group_list)

DGElist <- DGEList(counts=count, group = group_list)


DGElist <- calcNormFactors(DGElist,method="upperquartile")
DGElist1 <- estimateGLMCommonDisp(DGElist, design)
DGElist1 <- estimateGLMTrendedDisp(DGElist1, design)
DGElist1 <- estimateGLMTagwiseDisp(DGElist1, design)

fit <- glmFit(DGElist1, design)

results <- glmLRT(fit, contrast = c(-1, 1)) 

summary(de<-decideTestsDGE(results))

nrDEG_edgeR <- topTags(results, n = nrow(DGElist1))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
head(nrDEG_edgeR)


edgeRDEGs  <- nrDEG_edgeR[(nrDEG_edgeR$PValue < pvalue & 
                             (nrDEG_edgeR$logFC>foldChange | nrDEG_edgeR$logFC<(-foldChange))),]
edgeRDEGs <- edgeRDEGs[order(edgeRDEGs$logFC),]



###########

library(limma)


dt <- voom(DGElist, design, plot=TRUE)
#VOOM returns normalized data in a log scale
boxplot(dt$E)

design <- model.matrix(~0+factor(group_list))

colnames(design) <-levels(factor(group_list))

rownames(design) <- colnames(dt)

design

fit <- lmFit(dt,design)

fit2 <- eBayes(fit)  

limmaDEGs <- topTable(fit2, coef = 2,
                      p.value = 0.05, 
                      lfc = log2(2), 
                      number = Inf,
                      sort.by="logFC", 
                      adjust = "BH")


##################################

library(DESeq2)

#mycounts<-mycounts[rowMeans(mycounts)>0.01,] 
goup_list<-as.vector(lable$lable)

condition <- factor(goup_list, levels = c("Experimental","Control"))

colData <- data.frame(row.names=colnames(DGElist$counts), condition)
#colData <- data.frame(row.names=colnames(mycounts), condition2)

dds <- DESeqDataSetFromMatrix(countData = DGElist$counts,
                              colData = colData,
                              design = ~condition)
dds$condition<- relevel(dds$condition, ref = "Control") # ??????????????????????????????

dds <- DESeq(dds)

dds
allDEG2 <- as.data.frame(results(dds))
allDEG2<-allDEG2[allDEG2$baseMean>0,] 


diff_signif2 <- allDEG2[(allDEG2$pvalue < pvalue & 
                           (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
DESeq2DEGs <- diff_signif2[order(diff_signif2$log2FoldChange),]



library(VennDiagram)

venn.diagram(list(edgeR = rownames(edgeRDEGs),
                  DEseq2 = rownames(DESeq2DEGs),
                  Limma = rownames(limmaDEGs)),
             filename = "GSE124326FC.tiff", height = 450,width = 450,
             resolution =600,imagetype="tiff",col="transparent",
             fill=c("cornflowerblue","gray","orchid3"),alpha = 0.70,
             cex=0.45,cat.cex=0.45)
