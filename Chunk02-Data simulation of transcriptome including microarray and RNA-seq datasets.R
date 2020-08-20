
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
### code chunk number 02: Data simulation of transcriptome.
### ****************************************************************************

### R packages used in this step: madsim, compcodeR, and LIMMA. 
# For windows, work directory. 
# setwd("F:/") 

### ------------------------------------------------------------------------ ###
### Step-01. Define the number of samples and genes, respectively.

num.control <- 15

num.experimental <- 15

num.gene <- 20000

### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. Data simulation of DNA Microarray using madsim package.

# library(madsim)

fparams <- data.frame(m1 = num.control, # the number of samples in control group
                      m2 = num.experimental, # the number of samples in experimental group
                      shape2 = 4, # shape2 refers to the beta distribution shape
                      lb = 4, # lower value
                      ub = 14, # upper value
                      pde = 0.025, # percentage of DEGs
                      sym = 0.5) # partition of down- and up-regulated genes

dparams <- data.frame(lambda1 = 0.13, # Gene Average Level Variation Range
                      lambda2 = 2, # Fold Change Variation Parameters 
                      muminde = 1, # the mu-min-DE (differential genes) 
                      sdde = 0.5) # the sd of DE

sdn <- 0.4 # normal distribution standard deviation for additive noise

rseed <- 50 # computer random number initialization

mcr.data <- madsim.yhc(mdata = NULL, 
                   n = num.gene, 
                   ratio = 0, 
                   fparams, 
                   dparams, 
                   sdn, 
                   rseed)

# DNA microarray dataset simulated in this study.

mcr.matrix <- mcr.data$xdata
dim(mcr.matrix)
rownames(mcr.matrix) <- paste0("g", 1:nrow(mcr.matrix))

# checking the up- and down-regulated genes in mcr.matrix. 

up.id <- which(mcr.data$xid == 1)
down.id <- which(mcr.data$xid == -1)

table(mcr.data$xid)

mcr.file <- paste0("mcr", 
                   "-", 
                   num.control, "C", 
                   "-", 
                   num.experimental, "E", 
                   "-", 
                   num.gene, "g", 
                   "-", 
                   "matrix.RData")

save(mcr.matrix, file = mcr.file)

### End of Step-02.
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-03. Data simulation of RNA-sequencing using compcodeR package.

# Warnning: you must run the module chunk01 before this step. 
# tmpdir <- normalizePath(tempdir(), winslash = "/")

# library(compcodeR)

seqdata.obj <- generateSyntheticData.yhc(dataset = "seq.data", 
                                         n.vars = num.gene, # the number of genes
                                         m1 = num.experimental, # number of samples in experimental group
                                         m2 = num.control, # number of samples in control group
                                         n.diffexp = 500,
                                         fraction.upregulated = 0.5,
                                         output.file = "seqdata.rds")

# Extracting the gene expression matrix. 

seq.matrix <- seqdata.obj@count.matrix

seq.varanno <- seqdata.obj@variable.annotations

seq.samanno <- seqdata.obj@sample.annotations

seq.upvar <- table(seq.varanno$upregulation)

up.gene <- rownames(seq.matrix)[seq.varanno$upregulation == 1]

seq.downvar <- table(seq.varanno$downregulation)

down.gene <- rownames(seq.matrix)[seq.varanno$downregulation == 1]

dim(seq.matrix)

seq.file <- paste0("seq", 
                   "-", 
                   num.control, "C", 
                   "-", 
                   num.experimental, "E", 
                   "-", 
                   num.gene, "g", 
                   "-", 
                   "matrix.RData")

save(seq.matrix, file = seq.file)

### End of Step-03. 
### ------------------------------------------------------------------------ ###

### End of this chunk. 
### ****************************************************************************
