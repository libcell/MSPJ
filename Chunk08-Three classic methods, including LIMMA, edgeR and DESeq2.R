
################################################################################
#    &&&....&&&    % Project: MSPJ approach for identification of DEGs         #
#  &&&&&&..&&&&&&  % Author: Bo Li, Huachun Yin, Jingxin Tao, Youjin Hao       #
#  &&&&&&&&&&&&&&  % Date: Jun. 1st, 2020                                      #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.5.3;                             #
#       &&&&       % Platform: x86_64-pc-linux-gnu (64-bit)                    #
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

count <- get(load("count.RData"))

count[1:6, 1:6]

lable <- get(load("lable.RData"))

group <- as.factor(lable$lable)

# filtering the counts with low values 

cpms  <- cpm(count)

keep <- rowSums(cpms>1) >= 3

count <- count[keep, ]

# Generating the DGEList object

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

dim(deg.edgeR)

### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. LIMMA method for RNA-seq data. 

eset <- count

library(limma)

dt <- voom(norm.y, design, plot = FALSE)









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
  
  eset <- DGElist$count  
  
}










fit <- lmFit(dt, design)

fit2 <- eBayes(fit, trend = FALSE)  

limmaDEGs <- topTable(fit2, coef = 2, number = Inf)

deg.limma <- limmaDEGs[abs(limmaDEGs$logFC) > 0.58 & limmaDEGs$adj.P.Val < 0.05, ]

deg.limma <- rownames(deg.limma)

deg.limma

### End of Step-02.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-03. DESeq2 method for RNA-seq data. 

library(DESeq2)

group <- as.vector(lable$lable)

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

nrow(deg.deseq2)

### End of Step-03.
### ------------------------------------------------------------------------ ###
