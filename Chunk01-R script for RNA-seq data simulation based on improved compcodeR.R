
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
### code chunk number 01: Improved R script for RNA-seq data simulation.
### ****************************************************************************

# Setting the work directory, under the Windows Operation System. 
setwd("E:/R-work-dir")

### ------------------------------------------------------------------------ ###
### Step-01. Install all R packages required in this step. 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("edgeR")) BiocManager::install("edgeR")
if (!require("compcodeR")) BiocManager::install("compcodeR")
if (!require("preprocessCore")) BiocManager::install("preprocessCore")
if (!require("meta")) install.packages("meta")
if (!require("e1071")) install.packages("e1071")
if (!require("coin")) install.packages("coin")
if (!require("ggsci")) install.packages("ggsci")

library(preprocessCore)
library(madsim)
library(edgeR)
library(compcodeR)
library(meta)
library(e1071)
library(coin)
library(venn)
library(ggsci)

### End of Step-01. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-01. Defined the function madsim.yhc (for microarray data simulation) 

madsim.yhc <- function (mdata = NULL, 
                        n = 10000, 
                        ratio = 0, 
                        fparams = data.frame(m1 = 7, 
                                             m2 = 7, 
                                             shape2 = 4, 
                                             lb = 4, 
                                             ub = 14, 
                                             pde = 0.02, 
                                             sym = 0.5), 
                        dparams = data.frame(lambda1 = 0.13, 
                                             lambda2 = 2, 
                                             muminde = 1, 
                                             sdde = 0.5), 
                        sdn = 0.4, 
                        rseed = 50) {
  set.seed(rseed)
  m1 <- fparams$m1
  m2 <- fparams$m2
  m <- m1 + m2 + 1
  n1 <- length(mdata)
  if (n < 100) 
    n <- 100
  if (n1 > 100) {
    if (n == 0) {
      n <- n1
      x2 <- mdata
    }
    else {
      if (n < n1) {
        x2 <- sample(mdata, n, replace = FALSE)
      }
      else {
        x2 <- sample(mdata, n, replace = TRUE)
      }
    }
  }
  else {
    x <- rbeta(n, 2, fparams$shape2)
    x2 <- fparams$lb + (fparams$ub * x)
  }
  xdat <- matrix(c(rep(0, n * m)), ncol = m)
  xid <- matrix(c(rep(0, n)), ncol = 1)
  for (i in 1:n) {
    alpha <- dparams$lambda1 * exp(-dparams$lambda1 * x2[i])
    xi_val <- runif(m, min = (1 - alpha) * x2[i], max = (1 + 
                                                           alpha) * x2[i])
    if (sample(1:1000, 1) > (1000 * fparams$pde)) {
      xdat[i, ] <- xi_val
    }
    else {
      xi1 <- xi_val[1:(m1 + 1)]
      mude <- dparams$muminde + rexp(1, dparams$lambda2)
      if (sample(1:1000, 1) > (1000 * fparams$sym)) {
        xi2 <- xi_val[(m1 + 2):m] + rnorm(m2, mean = mude, 
                                          sd = dparams$sdde)
        xid[i] <- 1
      }
      else {
        xi2 <- xi_val[(m1 + 2):m] - rnorm(m2, mean = mude, 
                                          sd = dparams$sdde)
        xid[i] <- -1
      }
      xdat[i, ] <- c(xi1, xi2)
    }
  }
  cont <- paste("Control", 1:m1, sep = "-")
  test <- paste("Experimental", 1:m2, sep = "-")
  colnames(xdat) <- c("V1", cont, test)
  xsd <- sd(xdat[, 1])
  if (sdn > 0) {
    ndat <- matrix(c(rnorm(n * m, mean = 0, sd = sdn)), ncol = m)
    xdat <- xdat + ndat
  }
  xdata <- matrix(c(rep(0, n * (m - 1))), ncol = (m - 1))
  if (ratio) {
    xdata <- xdat[, 2:m] - xdat[, 1]
  }
  else {
    xdata <- xdat[, 2:m]
  }
  list(xdata = xdata, xid = xid, xsd = xsd)
}

### End of Function-01. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-02. Defined the function generateSyntheticData.yhc (for RNA-seq data simulation)

library(edgeR)
library(compcodeR)

generateSyntheticData.yhc <- function(dataset, 
                                      n.vars, 
                                      m1, 
                                      m2, 
                                      n.diffexp, 
                                      repl.id = 1, 
                                      seqdepth = 1e+07, 
                                      minfact = 0.7, 
                                      maxfact = 1.4, 
                                      relmeans = "auto", 
                                      dispersions = "auto", 
                                      fraction.upregulated = 1, 
                                      between.group.diffdisp = FALSE, 
                                      filter.threshold.total = 0, 
                                      filter.threshold.mediancpm = 0, 
                                      fraction.non.overdispersed = 0, 
                                      random.outlier.high.prob = 0, 
                                      random.outlier.low.prob = 0, 
                                      single.outlier.high.prob = 0, 
                                      single.outlier.low.prob = 0, 
                                      effect.size = 1.5, 
                                      output.file = NULL) {
  
  
  if (!is.null(output.file)) {
    if (!(substr(output.file, nchar(output.file) - 3, nchar(output.file)) == 
          ".rds")) {
      stop("output.file must be an .rds file.")
    }
  }
  uID <- paste(sample(c(0:9, letters, LETTERS), 10, replace = TRUE), 
               collapse = "")
  condition <- rep(c(1, 2), times = c(m1, m2))
  S1 <- which(condition == 1)
  S2 <- which(condition == 2)
  if (length(effect.size) == 1) {
    n.upregulated <- floor(fraction.upregulated * n.diffexp)
    if (fraction.upregulated != 0 & n.diffexp != 0) {
      genes.upreg <- 1:n.upregulated
    }
    else {
      genes.upreg <- NULL
    }
    if (fraction.upregulated != 1 & n.diffexp != 0) {
      genes.downreg <- (n.upregulated + 1):n.diffexp
    }
    else {
      genes.downreg <- NULL
    }
    genes.nonreg <- setdiff(1:n.vars, union(genes.upreg, 
                                            genes.downreg))
  }
  else {
    if (length(effect.size) != n.vars) {
      stop("The length of the effect.size vector must be the same as the number of simulated genes.")
    }
    else {
      genes.upreg <- which(effect.size > 1)
      genes.downreg <- which(effect.size < 1)
      genes.nonreg <- which(effect.size == 1)
      n.upregulated <- length(genes.upreg)
      n.diffexp <- length(genes.upreg) + length(genes.downreg)
      fraction.upregulated <- n.upregulated/n.diffexp
    }
  }
  differential.expression <- rep(0, n.vars)
  differential.expression[genes.upreg] <- 1
  differential.expression[genes.downreg] <- 1
  upregulation <- rep(0, n.vars)
  upregulation[genes.upreg] <- 1
  downregulation <- rep(0, n.vars)
  downregulation[genes.downreg] <- 1
  if (is.character(relmeans) | is.character(dispersions)) {
    mu.phi.estimates <- system.file("extdata", "Pickrell.Cheung.Mu.Phi.Estimates.rds", 
                                    package = "compcodeR")
    mu.phi.estimates <- readRDS(mu.phi.estimates)
    mu.estimates <- mu.phi.estimates$pickrell.cheung.mu
    phi.estimates <- mu.phi.estimates$pickrell.cheung.phi
    to.include <- sample(1:length(mu.estimates), n.vars, 
                         replace = ifelse(n.vars > length(mu.estimates), TRUE, 
                                          FALSE))
    truedispersions.S1 <- phi.estimates[to.include]
    truemeans.S1 <- mu.estimates[to.include]
  }
  if (!is.character(relmeans)) {
    if (length(relmeans) != n.vars) 
      stop("The length of the relmeans vector must be the same as the number of simulated genes.")
    truemeans.S1 <- c(relmeans)
  }
  if (!is.character(dispersions)) {
    if (nrow(cbind(dispersions)) != n.vars) 
      stop("The number of provided dispersions must be the same as the number of simulated genes.")
    truedispersions.S1 <- cbind(dispersions)[, 1]
    if (ncol(cbind(dispersions)) > 1) {
      truedispersions.S2 <- cbind(dispersions)[, 2]
    }
    else {
      truedispersions.S2 <- truedispersions.S1
    }
  }
  nfacts <- runif(m1+m2, min = minfact, max = maxfact)
  seq.depths <- nfacts * seqdepth
  overdispersed <- rep(1, n.vars)
  if (fraction.non.overdispersed > 0) {
    overdispersed[genes.upreg[1:round(fraction.non.overdispersed * 
                                        length(genes.upreg))]] <- 0
    overdispersed[genes.downreg[1:round(fraction.non.overdispersed * 
                                          length(genes.downreg))]] <- 0
    overdispersed[genes.nonreg[1:round(fraction.non.overdispersed * 
                                         length(genes.nonreg))]] <- 0
  }
  prob.S1 <- truemeans.S1
  prob.S2 <- rep(0, length(prob.S1))
  if (length(effect.size) == 1) {
    for (i in 1:n.vars) {
      if (i %in% genes.upreg) {
        prob.S2[i] <- (effect.size + rexp(1, rate = 1)) * 
          prob.S1[i]
      }
      else {
        if (i %in% genes.downreg) {
          prob.S2[i] <- 1/(effect.size + rexp(1, rate = 1)) * 
            prob.S1[i]
        }
        else {
          prob.S2[i] <- prob.S1[i]
        }
      }
    }
  }
  else {
    prob.S2 <- c(effect.size) * prob.S1
  }
  true.log2foldchange <- log2(prob.S2/prob.S1)
  sum.S1 <- sum(prob.S1)
  sum.S2 <- sum(prob.S2)
  if (is.character(dispersions)) {
    truedispersions.S2 <- truedispersions.S1
    if (between.group.diffdisp == TRUE) {
      for (i in 1:length(truedispersions.S2)) {
        sample.base <- phi.estimates[abs(log10(mu.estimates) - 
                                           log10(prob.S2[i])) < 0.05]
        if (length(sample.base) < 50) {
          sample.base <- phi.estimates[order(abs(log10(mu.estimates) - 
                                                   log10(prob.S2[i])))][1:500]
        }
        truedispersions.S2[i] <- sample(sample.base, 
                                        1)
      }
    }
  }
  truedispersions.S1 <- truedispersions.S1 * overdispersed
  truedispersions.S2 <- truedispersions.S2 * overdispersed
  Z <- matrix(0, n.vars, length(S1) + length(S2))
  for (i in 1:n.vars) {
    for (j in 1:ncol(Z)) {
      if (j %in% S1) {
        if (overdispersed[i] == 1) {
          Z[i, j] <- rnbinom(n = 1, mu = prob.S1[i]/sum.S1 * 
                               seq.depths[j], size = 1/truedispersions.S1[i])
        }
        else {
          Z[i, j] <- rpois(n = 1, lambda = prob.S1[i]/sum.S1 * 
                             seq.depths[j])
        }
      }
      else {
        if (overdispersed[i] == 1) {
          Z[i, j] <- rnbinom(n = 1, mu = prob.S2[i]/sum.S2 * 
                               seq.depths[j], size = 1/truedispersions.S2[i])
        }
        else {
          Z[i, j] <- rpois(n = 1, lambda = prob.S2[i]/sum.S2 * 
                             seq.depths[j])
        }
      }
    }
  }
  random.outliers <- matrix(0, nrow(Z), ncol(Z))
  random.outliers.factor <- matrix(1, nrow(Z), ncol(Z))
  if (random.outlier.high.prob != 0 | random.outlier.low.prob != 
      0) {
    for (i in 1:nrow(Z)) {
      for (j in 1:ncol(Z)) {
        tmp <- runif(1)
        if (tmp < random.outlier.high.prob) {
          random.outliers[i, j] <- 1
          random.outliers.factor[i, j] <- runif(1, min = 5, 
                                                max = 10)
        }
        else if (tmp < random.outlier.low.prob + random.outlier.high.prob) {
          random.outliers[i, j] <- (-1)
          random.outliers.factor[i, j] <- 1/runif(1, 
                                                  min = 5, max = 10)
        }
      }
    }
    Z <- round(random.outliers.factor * Z)
  }
  has.single.outlier <- rep(0, n.vars)
  single.outliers <- matrix(0, nrow(Z), ncol(Z))
  single.outliers.factor <- matrix(1, nrow(Z), ncol(Z))
  if (single.outlier.high.prob != 0 | single.outlier.low.prob != 
      0) {
    has.single.outlier[genes.upreg[1:floor((single.outlier.high.prob + 
                                              single.outlier.low.prob) * length(genes.upreg))]] <- 1
    has.single.outlier[genes.downreg[1:floor((single.outlier.high.prob + 
                                                single.outlier.low.prob) * length(genes.downreg))]] <- 1
    has.single.outlier[genes.nonreg[1:floor((single.outlier.high.prob + 
                                               single.outlier.low.prob) * length(genes.nonreg))]] <- 1
    for (i in 1:nrow(Z)) {
      if (has.single.outlier[i] == 1) {
        the.sample <- sample(1:(ncol(Z)), 1)
        if (runif(1) < (single.outlier.high.prob/(single.outlier.high.prob + 
                                                  single.outlier.low.prob))) {
          single.outliers[i, the.sample] <- 1
          single.outliers.factor[i, the.sample] <- runif(1, 
                                                         min = 5, max = 10)
        }
        else {
          single.outliers[i, the.sample] <- (-1)
          single.outliers.factor[i, the.sample] <- 1/runif(1, 
                                                           min = 5, max = 10)
        }
      }
    }
    Z <- round(single.outliers.factor * Z)
  }
  rownames(Z) <- 1:n.vars
  n.random.outliers.up.S1 <- apply(random.outliers[, S1] > 0, 1, sum)
  n.random.outliers.up.S2 <- apply(random.outliers[, S2] > 0, 1, sum)
  n.random.outliers.down.S1 <- apply(random.outliers[, S1] < 0, 1, sum)
  n.random.outliers.down.S2 <- apply(random.outliers[, S2] < 0, 1, sum)
  n.single.outliers.up.S1 <- apply(single.outliers[, S1] > 0, 1, sum)
  n.single.outliers.up.S2 <- apply(single.outliers[, S2] > 0, 1, sum)
  n.single.outliers.down.S1 <- apply(single.outliers[, S1] < 0, 1, sum)
  n.single.outliers.down.S2 <- apply(single.outliers[, S2] < 0, 1, sum)
  nf <- calcNormFactors(Z)
  norm.factors <- nf * colSums(Z)
  common.libsize <- exp(mean(log(colSums(Z))))
  pseudocounts <- sweep(Z + 0.5, 2, norm.factors, "/") * common.libsize
  log2.pseudocounts <- log2(pseudocounts)
  M.value <- apply(log2.pseudocounts[, S2], 1, mean) - apply(log2.pseudocounts[, S1], 1, mean)
  A.value <- 0.5 * (apply(log2.pseudocounts[, S2], 1, mean) + 
                      apply(log2.pseudocounts[, S1], 1, mean))
  variable.annotations <- data.frame(truedispersions.S1 = truedispersions.S1, 
                                     truedispersions.S2 = truedispersions.S2, truemeans.S1 = prob.S1, 
                                     truemeans.S2 = prob.S2, n.random.outliers.up.S1 = n.random.outliers.up.S1, 
                                     n.random.outliers.up.S2 = n.random.outliers.up.S2, n.random.outliers.down.S1 = n.random.outliers.down.S1, 
                                     n.random.outliers.down.S2 = n.random.outliers.down.S2, 
                                     n.single.outliers.up.S1 = n.single.outliers.up.S1, n.single.outliers.up.S2 = n.single.outliers.up.S2, 
                                     n.single.outliers.down.S1 = n.single.outliers.down.S1, 
                                     n.single.outliers.down.S2 = n.single.outliers.down.S2, 
                                     M.value = M.value, A.value = A.value, truelog2foldchanges = true.log2foldchange, 
                                     upregulation = upregulation, downregulation = downregulation, 
                                     differential.expression = differential.expression)
  rownames(variable.annotations) <- rownames(Z)
  sample.annotations <- data.frame(condition = condition, depth.factor = nfacts)
  info.parameters <- list(n.diffexp = n.diffexp, fraction.upregulated = fraction.upregulated, 
                          between.group.diffdisp = between.group.diffdisp, filter.threshold.total = filter.threshold.total, 
                          filter.threshold.mediancpm = filter.threshold.mediancpm, 
                          fraction.non.overdispersed = fraction.non.overdispersed, 
                          random.outlier.high.prob = random.outlier.high.prob, 
                          random.outlier.low.prob = random.outlier.low.prob, single.outlier.high.prob = single.outlier.high.prob, 
                          single.outlier.low.prob = single.outlier.low.prob, effect.size = effect.size, 
                          m1=m1, m2=m2, repl.id = repl.id, 
                          dataset = dataset, uID = uID, seqdepth = seqdepth, minfact = minfact, 
                          maxfact = maxfact)
  s <- apply(Z, 1, sum)
  keep.T <- which(s >= filter.threshold.total)
  Z.T <- Z[keep.T, ]
  variable.annotations.T <- variable.annotations[keep.T, ]
  filtering <- paste("total count >=", filter.threshold.total)
  cpm <- sweep(Z.T, 2, apply(Z.T, 2, sum), "/") * 1e+06
  m <- apply(cpm, 1, median)
  keep.C <- which(m >= filter.threshold.mediancpm)
  Z.TC <- Z.T[keep.C, ]
  variable.annotations.TC <- variable.annotations.T[keep.C, 
                                                    ]
  filtering <- paste(filtering, "; ", paste("median cpm >=", 
                                            filter.threshold.mediancpm))
  rownames(Z.TC) <- paste("g", 1:nrow(Z.TC), sep = "")
  cont <- paste("Experimental", 1:m1, sep = "-") # naming the experimental samples
  test <- paste("Control", 1:m2, sep = "-") # naming the control samples
  colnames(Z.TC) <- c(cont, test)
  rownames(sample.annotations) <- colnames(Z.TC)
  rownames(variable.annotations.TC) <- rownames(Z.TC)
  data.object <- compData(count.matrix = Z.TC, variable.annotations = variable.annotations.TC, 
                          sample.annotations = sample.annotations, filtering = filtering, 
                          info.parameters = info.parameters)
  if (!is.null(output.file)) {
    saveRDS(data.object, file = output.file)
  }
  return(invisible(data.object))
  
}

### End of Function-02. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-03. Define the function to generate multiple sub-groups from primary study. 

generateSubGroup <- function(dataset = eset, 
                             set.n = 40, 
                             size.min = 10, 
                             size.max = 20) {
  
  if (!is.matrix(dataset)) {
    
    stop("Please input the propriate dataset!")
    
  }
  
  sample.sets <- list()
  
  for (i in 1:set.n) {
    
    sam.index <- sample(1:ncol(eset), 
                        sample(size.min:size.max, 1, replace = FALSE), 
                        replace = FALSE)
    
    if (length(grep("Experimental", colnames(eset[, sam.index]))) >= 3 & length(grep("Control", colnames(eset[, sam.index]))) >= 3) 
      sample.sets[[i]] <- eset[, sam.index] else {
        
        sample.sets[[i]] <- NA
        
        next
        
      }  
    
  }
  
  real.sample.sets <- sample.sets[!is.na(sample.sets)]
  
  names(real.sample.sets) <- paste("sampling_set", 1:length(real.sample.sets), sep = "-")
  
  return(real.sample.sets)
  
}

# dim(subgroups[[5]])

### End of Function-03. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-04. Define the function to identify the DEGs based on meta-analysis.  

batchMeta <- function(data.list, cutoff = 0.5, g.start = 1, g.end = 100) {
  
  # Parameter 1: determin the threhold of SMDs. 
  
  na.index <- NULL
  up.index <- NULL
  down.index <- NULL
  
  i <- 0
  
  # Display the progress bar!
  pb <- txtProgressBar(min = 0, max = g.end - g.start + 1, style = 3, char = "+")
  
  for (ord.gene in g.start:g.end) {
    
    # Sys.sleep(0.1)
    
    i <- i + 1
    ### ord.gene <- 10
    
    stat.mat <- data.frame(matrix(NA, set.n, 8))
    
    names(stat.mat) <- c("study", "year", 
                         "n.e", "mean.e", "sd.e", 
                         "n.c", "mean.c", "sd.c")
    
    stat.mat$study <- paste("sampling_set", 1:set.n, sep = "-")
    
    stat.mat$year <- sample(2000:2020, set.n, replace = TRUE)
    
    # x <- sample.sets[[1]]
    
    f.ne <- function(x) length(grep("Experimental", colnames(x)))
    
    f.meane <- function(x) mean(x[ord.gene, grep("Experimental", colnames(x))])
    
    f.sde <- function(x) sd(x[ord.gene, grep("Experimental", colnames(x))])
    
    f.nc <- function(x) length(grep("Control", colnames(x)))
    
    f.meanc <- function(x) mean(x[ord.gene, grep("Control", colnames(x))])
    
    f.sdc <- function(x) sd(x[ord.gene, grep("Control", colnames(x))])
    
    stat.mat$n.e <- unlist(lapply(data.list, f.ne))
    stat.mat$n.c <- unlist(lapply(data.list, f.nc))
    
    stat.mat$mean.e <- unlist(lapply(data.list, f.meane))
    stat.mat$mean.c <- unlist(lapply(data.list, f.meanc))
    
    stat.mat$sd.e <- unlist(lapply(data.list, f.sde))
    stat.mat$sd.c <- unlist(lapply(data.list, f.sdc))
    
    # DT::datatable(stat.mat)
    
    ### ------------------------------------------------------------------------ ###
    ### Step-05. Implementation of meta-analysis for a specific gene (ord.gene).
    
    # Forest plot for a given gene (ord.gene). 
    
    res <- metacont(n.e, # Number of observations in experimental group
                    mean.e, # Estimated mean in experimental group
                    sd.e, # Standard deviation in experimental group
                    n.c, # Number of observations in control group
                    mean.c, # Estimated mean in control group
                    sd.c, # Standard deviation in control group
                    studlab = study, # An optional vector with study labels
                    data = stat.mat, # Data frame containing the study information
                    sm = "SMD")  # One of three measures ("MD", "SMD" and "ROM")
    
    # Extracting the detail model parameters: 
    
    # - The SMD statistics in random effect model. 
    
    zTE.r <- res$TE.random # Estimated treatment effect (TE) and standard error of individual studies
    lTE.r <- res$lower.random
    uTE.r <- res$upper.random
    
    # - The SMD statistics in fixed effect model. 
    
    zTE.f <- res$TE.fixed # Estimated treatment effect (TE) and standard error of individual studies
    lTE.f <- res$lower.fixed
    uTE.f <- res$upper.fixed
    
    # - Heterogeneity statistic I2, requiring < 50%. 
    Heter.I2 <- res$I2
    Heter.p <- res$pval.Q
    
    # Update progress bar
    # cat("\n")
    setTxtProgressBar(pb, i)
    
    # Finally, DEGs were identified by above indexes. 
    
    if (!is.na(lTE.r) & !is.na(lTE.f) & !is.na(uTE.r) & !is.na(uTE.f)) {
      
      if (lTE.r > cutoff & lTE.f > cutoff) {
        
        up.index <- c(up.index, ord.gene)
        
        # up.inf <- paste("The gene", ord.gene, "was up-regulated!", sep = " ")
        
        # print(up.inf)
        
      } else 
        
        if (uTE.r < -cutoff & uTE.f < -cutoff) {
          
          down.index <- c(down.index, ord.gene)
          
          # down.inf <- paste("The gene", ord.gene, "was down-regulated!", sep = " ")
          
          # print(down.inf)
          
        } else {
          
          # non.inf <- paste("The gene", ord.gene, "was not changed at the significant level with 0.05!", sep = " ")
          
          # print(non.inf)
          
        }
      
    } else {
      
      na.index <- c(na.index, ord.gene) ############
      
      next
      
    }
    
    flush.console()
  }
  
  MetaResult <- list(UpGene = up.index, DownGene = down.index)
  
  close(pb)
  
  return(MetaResult)
}

### End of Function-04. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-05. Define the function to generate the meta object.  

singleMeta <- function(data.list, ord.gene = 100) {
  
  stat.mat <- data.frame(matrix(NA, set.n, 8))
  
  names(stat.mat) <- c("study", "year", 
                       "n.e", "mean.e", "sd.e", 
                       "n.c", "mean.c", "sd.c")
  
  stat.mat$study <- paste("sampling_set", 1:set.n, sep = "-")
  
  stat.mat$year <- sample(2000:2020, set.n, replace = TRUE)
  
  # x <- sample.sets[[1]]
  
  f.ne <- function(x) length(grep("Experimental", colnames(x)))
  
  f.meane <- function(x) mean(x[ord.gene, grep("Experimental", colnames(x))])
  
  f.sde <- function(x) sd(x[ord.gene, grep("Experimental", colnames(x))])
  
  f.nc <- function(x) length(grep("Control", colnames(x)))
  
  f.meanc <- function(x) mean(x[ord.gene, grep("Control", colnames(x))])
  
  f.sdc <- function(x) sd(x[ord.gene, grep("Control", colnames(x))])
  
  stat.mat$n.e <- unlist(lapply(data.list, f.ne))
  stat.mat$n.c <- unlist(lapply(data.list, f.nc))
  
  stat.mat$mean.e <- unlist(lapply(data.list, f.meane))
  stat.mat$mean.c <- unlist(lapply(data.list, f.meanc))
  
  stat.mat$sd.e <- unlist(lapply(data.list, f.sde))
  stat.mat$sd.c <- unlist(lapply(data.list, f.sdc))
  
  res <- metacont(n.e, # Number of observations in experimental group
                  mean.e, # Estimated mean in experimental group
                  sd.e, # Standard deviation in experimental group
                  n.c, # Number of observations in control group
                  mean.c, # Estimated mean in control group
                  sd.c, # Standard deviation in control group
                  studlab = study, # An optional vector with study labels
                  data = stat.mat, # Data frame containing the study information
                  sm = "SMD")  # One of three measures ("MD", "SMD" and "ROM")
  
  return(res)
}

### End of Function-05. 
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Function-06. A set of SVM-RFE functions for gene selection. 
### Note: derived from johncolby/SVM-RFE/msvmRFE.R in Github. 

#++ Function-6.1 svmRFE.wrap

svmRFE.wrap <- function(test.fold, X, ...) {
  # Wrapper to run svmRFE function while omitting a given test fold
  train.data = X[-test.fold, ]
  test.data  = X[test.fold, ]
  
  # Rank the features
  features.ranked = svmRFE(train.data, ...)
  
  return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
}

#++ Function-6.2 svmRFE

svmRFE <- function(X, k = 1, halve.above = 5000) {
  # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
  n = ncol(X) - 1
  
  # Scale data up front so it doesn't have to be redone each pass
  cat('Scaling data...')
  X[, -1] = scale(X[, -1])
  cat('Done!\n')
  flush.console()
  
  pb = txtProgressBar(1, n, 1, style=3)
  
  i.surviving = 1:n
  i.ranked    = n
  ranked.list = vector(length=n)
  
  # Recurse through all the features
  while(length(i.surviving) > 0) {
    if(k > 1) {
      # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
      folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
      folds = lapply(1:k, function(x) which(folds == x))
      
      # Obtain weights for each training set
      w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
      w = do.call(rbind, w)
      
      # Normalize each weights vector
      w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
      
      # Compute ranking criteria
      v    = w * w
      vbar = apply(v, 2, mean)
      vsd  = apply(v, 2, sd)
      c    = vbar / vsd
    } else {
      # Only do 1 pass (i.e. regular SVM-RFE)
      w = getWeights(NULL, X[, c(1, 1+i.surviving)])
      c = w * w
    }
    
    # Rank the features
    ranking = sort(c, index.return=T)$ix
    if(length(i.surviving) == 1) {
      ranking = 1
    }
    
    if(length(i.surviving) > halve.above) {
      # Cut features in half until less than halve.above
      nfeat = length(i.surviving)
      ncut  = round(nfeat / 2)
      n     = nfeat - ncut
      
      cat('Features halved from', nfeat, 'to', n, '\n')
      flush.console()
      
      pb = txtProgressBar(1, n, 1, style=3)
      
    } else ncut = 1
    
    # Update feature list
    ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
    i.ranked    = i.ranked - ncut
    i.surviving = i.surviving[-ranking[1:ncut]]
    
    setTxtProgressBar(pb, n-length(i.surviving))
    flush.console()
  }
  
  close(pb)
  
  return (ranked.list)
}

#++ Function-6.3 getWeights

getWeights <- function(test.fold, X) {
  # Fit a linear SVM model and obtain feature weights
  train.data = X
  if(!is.null(test.fold)) train.data = X[-test.fold, ]
  
  svmModel = svm(train.data[, -1], train.data[, 1], cost=10, cachesize=500,
                 scale=F, type="C-classification", kernel="linear")
  
  t(svmModel$coefs) %*% svmModel$SV
}

#++ Function-6.4 WriteFeatures

WriteFeatures <- function(results, input, save=T, file='features_ranked.txt') {
  # Compile feature rankings across multiple folds
  featureID = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$ix
  avg.rank  = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$x
  feature.name = colnames(input[, -1])[featureID]
  features.ranked = data.frame(FeatureName=feature.name, FeatureID=featureID, AvgRank=avg.rank)
  if(save==T) {
    write.table(features.ranked, file=file, quote=F, row.names=F)
  } else {
    features.ranked
  }
}

#++ Function-6.5 FeatSweep.wrap

FeatSweep.wrap <- function(i, results, input) {
  # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
  svm.list = lapply(results, function(x) tune(svm,
                                              train.x      = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                              train.y      = input[x$train.data.ids, 1],
                                              validation.x = input[x$test.data.ids, 1+x$feature.ids[1:i]],
                                              validation.y = input[x$test.data.ids, 1],
                                              # Optimize SVM hyperparamters
                                              ranges       = tune(svm,
                                                                  train.x = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                                                  train.y = input[x$train.data.ids, 1],
                                                                  ranges  = list(gamma=2^(-12:0), cost=2^(-6:6)))$best.par,
                                              tunecontrol  = tune.control(sampling='fix'))$perf)
  
  error = mean(sapply(svm.list, function(x) x$error))
  return(list(svm.list=svm.list, error=error))
}

#++ Function-6.6 PlotErrors

PlotErrors <- function(errors, errors2=NULL, no.info=0.5, ylim=range(c(errors, errors2), na.rm=T), xlab='Number of Features',  ylab='10x CV Error') {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col='black') {
    lines(which(!is.na(errors)), na.omit(x), col=col)
    points(which.min(x), min(x, na.rm=T), col='red')
    text(which.min(x), min(x, na.rm=T), paste(which.min(x), '-', format(min(x, na.rm=T), dig=3)), pos=4, col='red', cex=0.75)
  }
  
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=3)
}

### End of Function-06. 
### ------------------------------------------------------------------------ ###

### End of this chunk. 
### ****************************************************************************
