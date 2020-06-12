
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
### code chunk number 02: Data simulation of transcriptome.
### ****************************************************************************

# - R package: madsim, compcodeR, and LIMMA. 

# For windows, work directory. 
# setwd("F:/")

### ------------------------------------------------------------------------ ###
### Step-01. Data simulation of DNA Microarray using madsim package.

fparams <- data.frame(m1 = 15, # the number of samples in control group
                      m2 = 15, # the number of samples in experimental group
                      shape2 = 4, # shape2 refers to the beta distribution shape
                      lb = 4, # 
                      ub = 14, # 
                      pde = 0.025, # percentage of DEGs
                      sym = 0.5) # partition of down- and up-regulated genes

dparams <- data.frame(lambda1 = 0.13, # Gene Average Level Variation Range
                      lambda2 = 2, # Fold Change Variation Parameters 
                      muminde = 1, # the mu-min-DE (differential genes) 
                      sdde = 0.5) # the sd of DE

sdn <- 0.4 # normal distribution standard deviation for additive noise

rseed <- 50 # computer random number initialization

mcr.data <- madsim.yhc(mdata = NULL, 
                   n = 20000, 
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

save(mcr.matrix, file = "mcr.matrix.RData")

### End of Step-01.
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-02. Data simulation of RNA-sequencing using compcodeR package.

# Warnning: you must run the module chunk01 before this step. 
# tmpdir <- normalizePath(tempdir(), winslash = "/")

library(compcodeR)

seqdata.obj <- generateSyntheticData.yhc(dataset = "seq.data", 
                                         n.vars = 20000, # the number of genes
                                         m1 = 15, # number of samples in experimental group
                                         m2 = 10, # number of samples in control group
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

save(seq.matrix, file = "seq.matrix.RData")

### End of Step-02. 
### ------------------------------------------------------------------------ ###
