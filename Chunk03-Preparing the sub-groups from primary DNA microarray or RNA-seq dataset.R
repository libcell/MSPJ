
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
### code chunk number 03: Generating the sub-groups based on resampling.
### ****************************************************************************

### Loading the gene expression matrix, in .RData format. 
#-- seq.matrix.RData: simulated RNA-seq data; 
#-- mcr.matrix.RData: simulated microarray data. 

### ------------------------------------------------------------------------ ###
### Step-01. Preparing the colors used in this study. 

mypal1 <- terrain.colors(num.control + num.experimental)

mypal2 <- pal_npg("nrc", alpha = 0.7)(10)

### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. Selecting the dataset, for DNA microarray and RNA-seq. 

### @@@@@@@@@@@@@@ Selecting the dataset (microarray/RNA-seq) @@@@@@@@@@@@@@ ###
### ======================================================================== ###
if (!exists("eset")) {
  
  eset <- mcr.matrix
  # eset <- seq.matrix
  
} 

print(eset[1:6, 1:6])
### ======================================================================== ###


# Visualizing the gene expression matrix. 

if (all(as.integer(eset) == as.numeric(eset))) {
  
  boxplot(log2(eset), col = mypal1, main = "RNA-sequencing data")
  
} else {
  
  boxplot(eset, col = mypal1, main = "DNA microarray data")
  
}

### ------------------------------------------------------------------------ ###
### Step-03. Generating multiple sub-groups based resampling for primary study. 

# ord.gene: which gene you focused on. 

sample.sets <- generateSubGroup(eset, 
                                set.n = 40, # the number of sub-groups
                                size.min = 10, # the lower limit of sample size
                                size.max = 20) # the maximum sample size

sample.sets[[1]][1:5, 1:5]

### End of Step-03.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-04. Preprocessing for DNA microarray or RNA-seq data, alternatively. 

# For all sub-datasets.  

for (d in 1:length(sample.sets)) {
  
  data <- sample.sets[[d]]
  
  if (all(as.integer(eset) == as.numeric(eset))) {
    
    #-- Data normalization for gene expression matrix filled by counts. 
    
    DGElist <- DGEList(counts = data)
    
    DGElist <- calcNormFactors(DGElist, method = "upperquartile")
    
    #. boxplot(log2(DGElist$count))
    
    # plotMDS(DGElist)
    
    data <- DGElist$count  
    
  } else {
    
    tmp <- normalize.quantiles(data)
    
    rownames(tmp) <- rownames(data)
    
    colnames(tmp) <- colnames(data)
    
    data <- tmp
  }
  
  sample.sets[[d]] <- data
  
}

sample.sets[[1]][1:5, 1:5]
#. 
### End of Step-04.
### ------------------------------------------------------------------------ ###
