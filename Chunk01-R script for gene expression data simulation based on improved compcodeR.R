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
### code chunk number 01: Improved R script for gene expression data simulation.
### ****************************************************************************

# Setting the work directory, under the Windows Operation System. 

setwd("..")
if (!dir.exists("MSPJ")) {
  dir.create("MSPJ")
  setwd("MSPJ")
} else {
  setwd("MSPJ")
}

# getwd()

### ------------------------------------------------------------------------ ###
### Step-01. Install all R packages required in this step. 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("edgeR")) BiocManager::install("edgeR")
if (!require("preprocessCore")) BiocManager::install("preprocessCore")
if (!require("meta")) install.packages("meta")
if (!require("e1071")) install.packages("e1071")
if (!require("coin")) install.packages("coin")
if (!require("ggsci")) install.packages("ggsci")
if (!require("SPsimSeq")) remotes::install_github("CenterForStatistics-UGent/SPsimSeq")


library(SPsimSeq)
library(preprocessCore)
library(madsim)
library(edgeR)
library(compcodeR)
library(meta)
library(e1071)
library(coin)
library(venn)
library(ggsci)
library(KEGGdzPathwaysGEO)
library(GEOquery)
library(affy)
library(limma)
library(pROC)
library(multtest)
library(preprocessCore)
library(FSA)
library(pROC)
library(limma)
library(samr)
library(simpleaffy)
library(RankProd)
library(bspec)
library(netClass)
library(randomForest)
library(FSinR)
library(multtest)
library(DT)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(Rmisc)

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
        xi2 <- xi_val[(m1 + 2):m] + rnorm(m2, 
                                          mean = mude, 
                                          sd = dparams$sdde)
        xid[i] <- 1
      }
      else {
        xi2 <- xi_val[(m1 + 2):m] - rnorm(m2, 
                                          mean = mude, 
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
    ndat <- matrix(c(rnorm(n * m, mean = 0, sd = sdn)), 
                   ncol = m)
    xdat <- xdat + ndat
  }
  xdata <- matrix(c(rep(0, n * (m - 1))), 
                  ncol = (m - 1))
  if (ratio) {
    xdata <- xdat[, 2:m] - xdat[, 1]
  }
  else {
    xdata <- xdat[, 2:m]
  }
  list(xdata = xdata, 
       xid = xid, 
       xsd = xsd)
}

### End of Function-01. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-02. Define the function to generate multiple sub-groups from primary study. 

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
                        sample(size.min:size.max, 
                               1, 
                               replace = FALSE), 
                        replace = FALSE)
    if (length(grep("Experimental", colnames(eset[, sam.index]))) >= 3 & length(grep("Control", colnames(eset[, sam.index]))) >= 3) 
      sample.sets[[i]] <- eset[, sam.index] else {
        sample.sets[[i]] <- NA
        next
      }  
    
  }
  real.sample.sets <- sample.sets[!is.na(sample.sets)]
  names(real.sample.sets) <- paste("sampling_set", 
                                   1:length(real.sample.sets), 
                                   sep = "-")
  return(real.sample.sets)
}
# dim(subgroups[[5]])

### End of Function-02. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-03. Define the function to identify the DEGs based on meta-analysis.  

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
    stat.mat$study <- paste("sampling_set", 
                            1:set.n, 
                            sep = "-")
    stat.mat$year <- sample(2000:2020, 
                            set.n, 
                            replace = TRUE)
    
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
  MetaResult <- list(UpGene = up.index, 
                     DownGene = down.index)
  close(pb)
  return(MetaResult)
}

### End of Function-03. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-04. Define the function to generate the meta object.  

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

### End of Function-04. 
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Function-05. A set of SVM-RFE functions for gene selection. 
### Note: derived from johncolby/SVM-RFE/msvmRFE.R in Github. 

#++ Function-5.1 svmRFE.wrap

svmRFE.wrap <- function(test.fold, X, ...) {
  # Wrapper to run svmRFE function while omitting a given test fold
  train.data = X[-test.fold, ]
  test.data  = X[test.fold, ]
  
  # Rank the features
  features.ranked = svmRFE(train.data, ...)
  
  return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
}

#++ Function-5.2 svmRFE

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

#++ Function-5.3 getWeights

getWeights <- function(test.fold, X) {
  # Fit a linear SVM model and obtain feature weights
  train.data = X
  if(!is.null(test.fold)) train.data = X[-test.fold, ]
  
  svmModel <- svm(
    train.data[,-1],
    train.data[, 1],
    cost = 10,
    cachesize = 500,
    scale = F,
    type = "C-classification",
    kernel = "linear"
  )
  
  t(svmModel$coefs) %*% svmModel$SV
}

#++ Function-5.4 WriteFeatures

WriteFeatures <-
  function(results,
           input,
           save = T,
           file = 'features_ranked.txt') {
    # Compile feature rankings across multiple folds
    featureID = sort(apply(sapply(results, function(x)
      sort(x$feature, index.return = T)$ix), 1, mean), index = T)$ix
    avg.rank  = sort(apply(sapply(results, function(x)
      sort(x$feature, index.return = T)$ix), 1, mean), index = T)$x
    feature.name = colnames(input[,-1])[featureID]
    features.ranked = data.frame(FeatureName = feature.name,
                                 FeatureID = featureID,
                                 AvgRank = avg.rank)
    if (save == T) {
      write.table(
        features.ranked,
        file = file,
        quote = F,
        row.names = F
      )
    } else {
      features.ranked
    }
  }

#++ Function-5.5 FeatSweep.wrap

FeatSweep.wrap <- function(i, results, input) {
  # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
  svm.list = lapply(results, function(x)
    tune(
      svm,
      train.x      = input[x$train.data.ids, 1 +
                             x$feature.ids[1:i]],
      train.y      = input[x$train.data.ids, 1],
      validation.x = input[x$test.data.ids, 1 +
                             x$feature.ids[1:i]],
      validation.y = input[x$test.data.ids, 1],
      # Optimize SVM hyperparamters
      ranges       = tune(
        svm,
        train.x = input[x$train.data.ids, 1 +
                          x$feature.ids[1:i]],
        train.y = input[x$train.data.ids, 1],
        ranges  = list(gamma =
                         2 ^ (-12:0), cost = 2 ^ (-6:6))
      )$best.par,
      tunecontrol  = tune.control(sampling =
                                    'fix')
    )$perf)
  
  error = mean(sapply(svm.list, function(x)
    x$error))
  return(list(svm.list = svm.list, error = error))
}

#++ Function-5.6 PlotErrors

PlotErrors <-
  function(errors,
           errors2 = NULL,
           no.info = 0.5,
           ylim = range(c(errors, errors2), na.rm = T),
           xlab = 'Number of Features',
           ylab = '10x CV Error') {
    # Makes a plot of average generalization error vs. number of top features
    AddLine <- function(x, col = 'black') {
      lines(which(!is.na(errors)), na.omit(x), col = col)
      points(which.min(x), min(x, na.rm = T), col = 'red')
      text(
        which.min(x),
        min(x, na.rm = T),
        paste(which.min(x), '-', format(min(x, na.rm = T), dig = 3)),
        pos = 4,
        col = 'red',
        cex = 0.75
      )
    }
    
    plot(
      errors,
      type = 'n',
      ylim = ylim,
      xlab = xlab,
      ylab = ylab
    )
    AddLine(errors)
    if (!is.null(errors2))
      AddLine(errors2, 'gray30')
    abline(h = no.info, lty = 3)
  }

geneticAlgorithm.yhc <-
  function (popSize = 20,
            pcrossover = 0.8,
            pmutation = 0.1,
            maxiter = 10,
            run = 10,
            verbose = FALSE)
  {
    geneticAlgorithmSearch <- function(data, class, featureSetEval) {
      if (attr(featureSetEval, "kind") == "Individual measure") {
        stop("Only feature set measures can be used")
      }
      column.names <- names(data)
      class.position <- which(column.names == class)
      features <- column.names[-class.position]
      metricTarget <- attr(featureSetEval, "target")
      if (metricTarget == "maximize") {
        max <- TRUE
      }
      else if (metricTarget == "minimize") {
        max <- FALSE
      }
      else {
        max <- ifelse(is.factor(data[, class]), TRUE, FALSE)
      }
      fitness <- function(ind, data, class) {
        feat <- features[which(ind == 1)]
        if (length(feat) == 0) {
          if (max) {
            value <- 0
          }
          else {
            value <- -Inf
          }
        }
        else {
          if (max) {
            value <- featureSetEval(data, class, feat)
          }
          else {
            value <- -(featureSetEval(data, class, feat))
          }
        }
        return(value)
      }
      if (!verbose) {
        capture.output(
          GA <- GA::ga(
            type = "binary",
            fitness = fitness,
            data = data,
            class = class,
            nBits = length(features),
            popSize = popSize,
            pcrossover = pcrossover,
            pmutation = pmutation,
            maxiter = maxiter,
            run = run
          )
        )
      }
      else {
        if (max) {
          GA <- GA::ga(
            type = "binary",
            fitness = fitness,
            data = data,
            class = class,
            nBits = length(features),
            popSize = popSize,
            pcrossover = pcrossover,
            pmutation = pmutation,
            maxiter = maxiter,
            run = run
          )
        }
        else {
          capture.output(
            GA <- GA::ga(
              type = "binary",
              fitness = fitness,
              data = data,
              class = class,
              nBits = length(features),
              popSize = popSize,
              pcrossover = pcrossover,
              pmutation = pmutation,
              maxiter = maxiter,
              run = run
            )
          )
        }
      }
      res <- list(NULL)
      res[[1]] <- GA@solution
      colnames(res[[1]]) <- features
      res[[2]] <- GA@fitnessValue
      if (!max) {
        res[[2]] <- -res[[2]]
      }
      ncol <- ncol(data)
      res[[3]] <- matrix(cbind(GA@population, GA@fitness),
                         ncol = ncol,
                         dimnames = list(c(), c(paste0("x", 1:((ncol(data) -
                                                                  1)
                         )), c("fitness"))))
      names(res) <- c("bestFeatures", "bestFitness", "population")
      if (!max) {
        res[[3]][, "fitness"] <- -res[[3]][, "fitness"]
      }
      res
    }
    attr(geneticAlgorithmSearch, "shortName") <- "ga"
    attr(geneticAlgorithmSearch, "name") <- "Genetic Algorithm"
    return(geneticAlgorithmSearch)
  }
### End of Function-05. 
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Function-06. SNR. 

snr <- function(x) {
  y1 <- mean(x[e.p]) - mean(x[c.p])
  y2 <- sd(x[e.p]) + sd(x[c.p])
  y <- y1 / y2
  return(y)
}

### End of Function-06. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Function-07. Jaccard index. 
###  

#####Function-7.1 the Jaccard score between A and B method

Jscore <- function(x, y) {
  int <- length(intersect(x, y))
  un <- length(union(x, y))
  z <- int / un
  return(z)
}

#####Function-7.2 the matrix of Jaccard score  between various methods

Jscaore.YHC <- function(x) {
  n <- length(x)
  nl <- gtools::permutations(n, 2, repeats = TRUE)
  m <- length(nl) / 2
  
  NAm <- matrix(nrow = n, ncol = n)
  
  colnames(NAm) <- names(x)
  rownames(NAm) <- names(x)
  for (i in 1:m) {
    
    col_index <- nl[i,1]
    row_index <- nl[i,2]
    
    X1 <- x[[col_index]]
    X2 <- x[[row_index]]
    
    R1 <- Jscore(X1, X2)
    
    NAm[row_index, col_index] <- R1
    
  }
  return(NAm)
}

### End of Function-07. 
### ------------------------------------------------------------------------ ###
