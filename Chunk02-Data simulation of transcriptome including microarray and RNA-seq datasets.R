
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
### code chunk number 01: Data simulation of transcriptome.
### ****************************************************************************

# - R package: madsim, compcodeR, and LIMMA. 

# For windows, work directory. 
# setwd("F:/")

### Step-01. Data simulation of DNA Microarray using madsim package.

library(madsim)

fparams <- data.frame(m1 = 15, # the number of samples in control group
                      m2 = 15, # the number of samples in case group
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

mcr.data <- madsim(mdata = NULL, 
                   n = 20000, 
                   ratio = 0, 
                   fparams, 
                   dparams, 
                   sdn, 
                   rseed)

# DNA microarray dataset simulated in this study.

mcr.matrix <- mcr.data$xdata
rownames(mcr.matrix) <- paste("gene", 1:nrow(mcr.matrix), sep = "-")

# checking the up- and down-regulated genes in mcr.matrix. 

table(mcr.data$xid)


### Step-02. Data simulation of RNA-sequencing using compcodeR package.

library(compcodeR)

# tmpdir <- normalizePath(tempdir(), winslash = "/")

seqdata.obj <- generateSyntheticData.yhc(dataset = "seq.data", 
                                         n.vars = 2000, # the number of genes
                                         m1 = 15, # 
                                         m2 = 10, # 
                                         n.diffexp = 500,
                                         fraction.upregulated = 0.5,
                                         output.file = "seqdata.rds")

seq.matrix <- seqdata.obj@count.matrix

seq.varanno <- seqdata.obj@variable.annotations
seq.samanno <- seqdata.obj@sample.annotations

seq.upvar <- table(seq.varanno$upregulation)

up.gene <- rownames(seq.matrix)[seq.varanno$upregulation == 1]

seq.downvar <- table(seq.varanno$downregulation)

down.gene <- rownames(seq.matrix)[seq.varanno$downregulation == 1]

dim(seq.matrix)

### End. 

