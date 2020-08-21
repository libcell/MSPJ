
################################################################################
#    &&&....&&&    % Project: MSPJ approach for identification of DEGs         #
#  &&&&&&..&&&&&&  % Author: Bo Li, Huachun Yin, Jingxin Tao, Youjin Hao       #
#  &&&&&&&&&&&&&&  % Date: Jun. 1st, 2020                                      #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.6.0;                             #
#       &&&&       % x86_64-w64-mingw32/x64 (64-bit)                           #
#        &         %                                                           #
################################################################################

### ****************************************************************************
### code chunk number 08: Gene selection based on three classic methods.
### ****************************************************************************

# input: eset - a gene expression matrix, genes in lines and samples in columns. 
# output: degs - a genelist, including the gene names, pooled SMDs, and so on. 

### ------------------------------------------------------------------------ ###
### Step-01. LIMMA method for DNA microarray and RNA-seq data. 



### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. edgeR method for RNA-seq data. 

### ------------------------------------------------------------------------ ###
### Step-01. Determing DEGs by using edgeR method.  

library("edgeR")

######################

# For simulated dataset. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

count <- eset

group <- input$sam.lab

cpms  <- cpm(count)

keep <- rowSums(cpms>1) >= 3

count <- count[keep, ]

y <- DGEList(counts = count, group = group)

# Data normalization

norm.y <- calcNormFactors(y, method = "upperquartile")

# preparing the design matrix

design <- model.matrix( ~ group)

# estimating the dispersion

y <- estimateDisp(norm.y, design, robust = TRUE)

y$common.dispersion

plotBCV(y)

# Differential expression by performing the likelihood ratio test. 

fit <- glmFit(y, design)

lrt <- glmLRT(fit, coef = 2)

degTable <- topTags(lrt, n = nrow(count))

degTable <- as.data.frame(degTable)

deg.edgeR <- degTable[abs(degTable$logFC) > 1 & degTable$PValue < 0.05, ]

deg.edgeR <- deg.edgeR[rev(order(abs(deg.edgeR$logFC))), ]

deg.edgeR <- rownames(deg.edgeR)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

######################



### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. LIMMA method for RNA-seq data. 

# eset <- count

library(limma)

if (all(as.integer(eset) == as.numeric(eset))) {
  
  #-- Filtering the genes with low-expression levels. 
  
  cpms <- cpm(eset)
  keep <- rowSums(cpms>1) >= 3
  eset <- eset[keep, ]
  
  #-- Data normalization for gene expression matrix filled by counts. 
  
  DGElist <- DGEList(counts = eset)
  
  DGElist <- calcNormFactors(DGElist, method = "upperquartile")
  
  #. boxplot(log2(DGElist$count))
  
  # plotMDS(DGElist)
  
  tmp <- DGElist$count  
  
  dt <- voom(tmp, design, plot = FALSE)
  
} else {
  
  dt <- eset
  
}


fit <- lmFit(dt, design)

fit2 <- eBayes(fit, trend = FALSE)  

limmaDEGs <- topTable(fit2, coef = 2, number = Inf)

deg.limma <- limmaDEGs[abs(limmaDEGs$logFC) > 1 & limmaDEGs$adj.P.Val < 0.05, ]

deg.limma <- deg.limma[rev(order(abs(deg.limma$logFC))), ]

deg.limma <- rownames(deg.limma)

### End of Step-02.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-03. DESeq2 method for RNA-seq data. 

###############################################################
### for simulated dataset. 
library(DESeq2)

group <- as.vector(group)

group[group == "Control"] <- "untreated"

condition <- as.factor(group)

colData <- data.frame(row.names = colnames(count), condition)

dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = colData,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, 
                         ref = "untreated") # Designated control group

dds <- DESeq(dds)

deg2 <- as.data.frame(results(dds))
deg2 <- deg2[deg2$baseMean > 0, ] 

deg.deseq2 <- deg2[abs(deg2$log2FoldChange) > 1 & deg2$pvalue < 0.05, ]

deg.deseq2 <- deg.deseq2[rev(order(abs(deg.deseq2$log2FoldChange))), ]

deg.deseq2 <- rownames(deg.deseq2)

#################################################################

### End of Step-03.
### ------------------------------------------------------------------------ ###

deg.final <- intersect(intersect(deg.limma, deg.edgeR), deg.deseq2)

intersect(deg.int, deg.final)

intersect(deg.int, deg.limma)

intersect(deg.int, deg.edgeR)

intersect(deg.int, deg.deseq2)

intersect(deg.deseq2, deg.edgeR)



### For us. 

deg.int <- deg.svm[sort(match(deg.int, deg.svm))]

### End of this chunk. 
### ****************************************************************************

