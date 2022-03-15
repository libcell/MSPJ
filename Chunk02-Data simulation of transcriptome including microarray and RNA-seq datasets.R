################################################################################
#    &&&....&&&    % Project: MSPJ approach for identification of DEGs         #
#  &&&&&&..&&&&&&  % Author: Bo Li, Huachun Yin, Jingxin Tao   
#  &&&&&&&&&&&&&&  % Date: Mar. 1st, 2022                                      #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.5.3;                             #
#       &&&&       % Platform: x86_64-pc-linux-gnu (64-bit)                    #
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
tot.samples<-num.control+num.experimental
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

library(SPsimSeq)

data("zhang.data.sub") 

# filter genes with sufficient expression (important step to avoid bugs) 
zhang.counts <- zhang.data.sub$counts 
MYCN.status  <- zhang.data.sub$MYCN.status #The grouping variable


set.seed(1) #Set seed for reproducibility
# simulate data
sim.data.bulk <- SPsimSeq(n.sim = 1,
                          s.data = zhang.counts,
                          group = MYCN.status,
                          n.genes = num.gene,
                          batch.config = 1,
                          group.config = c(0.5, 0.5),
                          tot.samples = tot.samples, 
                          pDE = 0.1,
                          lfc.thrld = 0.5,
                          result.format = "list",
                          return.details = FALSE)

# Extracting the gene expression matrix. 

seq.list <- sim.data.bulk[[1]]

seq.matrix <- seq.list$counts

seq.samanno <- seq.list$colData

seq.samanno$Group[seq.samanno$Group==1]<-"Experimental"

seq.samanno$Group[seq.samanno$Group==0]<-"Control"

colnames(seq.matrix)<-seq.samanno$Group

seq.rowanno <- seq.list$rowData

seq.DEG<-seq.rowanno[seq.rowanno$DE.ind==TRUE,]

seq.DEG<-rownames(seq.DEG)


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

### End of Chunk-02. 
### ****************************************************************************
