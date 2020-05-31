
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

tmpdir <- normalizePath(tempdir(), winslash = "/")

seqdata.obj <- generateSyntheticData(dataset = "seq.data", 
                                    n.vars = 20000, # the number of genes
                                    samples.per.cond = 15, # 
                                    n.diffexp = 500,
                                    fraction.upregulated = 0.5,
                                    output.file = "seqdata.rds")


library(compareDEtools)

GenerateSyntheticSimulation(working.dir = ".", 
                            #data.types = "KIRC", 
                            # rep = 10, 
                            nsample = c(15, 16), 
                            nvar = 2000, 
                            nDE = 100, 
                            fraction.upregulated = 0.5, 
                            disp.Types = 'same', 
                            modes = c("D"))


### End. 


