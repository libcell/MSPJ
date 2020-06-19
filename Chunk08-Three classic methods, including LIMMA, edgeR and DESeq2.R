
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
### Step-01. LIMMA method for RNA-seq data. 



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

count <- count[keep,]

# Generating the DGEList object

y <- DGEList(counts = count, group = group)

# Data normalization

y <- calcNormFactors(y, method = "upperquartile")

# preparing the design matrix

design <- model.matrix( ~ group)

# estimating the dispersion

y <- estimateDisp(y, design, robust = TRUE)

y$common.dispersion

plotBCV(y)

# Differential expression by performing the likelihood ratio test. 

fit <- glmFit(y, design)

lrt <- glmLRT(fit, coef = 1)

degTable <- topTags(lrt, n = nrow(count))

degTable <- as.data.frame(degTable)

deg.edgeR <- degTable[abs(degTable$logFC) > 0.5849625 & degTable$PValue < 0.05, ]

dim(deg.edgeR)

### End of Step-01.
### ------------------------------------------------------------------------ ###

### End of Step-02.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-03. DESeq2 method for RNA-seq data. 



### End of Step-03.
### ------------------------------------------------------------------------ ###
