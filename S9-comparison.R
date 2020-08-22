
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
### code chunk number 09: Other methods used for gene selection.
### ****************************************************************************

# input: eset - a gene expression matrix, genes in lines and samples in columns. 
# output: degs - a genelist, including the gene names, pooled SMDs, and so on. 

### ------------------------------------------------------------------------ ###
### Step-01. Preparing the packages used in gene selection. 

library(limma)
library(samr)
library(pROC)
library(simpleaffy)
library(FSA)

library(openxlsx)

data <- read.xlsx("Simulation_of_microarray_data.xlsx", 1)

c.split <- function(x) strsplit(x, ":")[[1]][1]
e.split <- function(x) strsplit(x, ":")[[1]][2]

num1 <- unlist(lapply(data$Sample, c.split))
num2 <- unlist(lapply(data$Sample, e.split))

red.ACC <- matrix(NA, ncol = 16)

colnames(red.ACC) <- c("simulation.num", "MSP.num", "MSP.AUC", "MSP.ACC", 
                       "limma.num", "limma.AUC", "limma.ACC", "SAM.num", 
                       "SAM.AUC", "SAM.ACC", "t-test.num", "t-test.AUC", 
                       "t-test.ACC", "t-test & FC.num", "t-test & FC.AUC", 
                       "t-test & FC.ACC")




for (num in 1:76) {
  
  num.control <- num1[num]
  
  num.experimental <- num2[num]
  
  num.gene <- 20000
  
  group<-c(rep("Control",num1[num]),rep("Experimental",num2[num]))
  
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
  
  rownames(mcr.matrix) <- paste0("g", 1:nrow(mcr.matrix))
  
  # checking the up- and down-regulated genes in mcr.matrix. 
  
  up.id <- which(mcr.data$xid == 1)
  down.id <- which(mcr.data$xid == -1)

  
  mcr.file <- paste0("mcr", 
                     "-", 
                     num.control, "C", 
                     "-", 
                     num.experimental, "E", 
                     "-", 
                     num.gene, "g", 
                     "-", 
                     "matrix.RData")
  
  eset <- mcr.matrix
  
  eset <- eset[which(apply(eset, 1, min) > 0.000001),]
  
  tmp1 <- normalize.quantiles(eset)
  
  rownames(tmp1) <- rownames(eset)
  
  colnames(tmp1) <- colnames(eset)
  
  eset <- tmp1
  
  
  ###########
  
  
  sample.sets <- generateSubGroup(eset, 
                                  set.n = 40, # the number of sub-groups
                                  size.min = 10, # the lower limit of sample size
                                  size.max = 20) # the maximum sample size
  
  set.n <- length(sample.sets)
  
  num.gene <- nrow(eset)
  
  # 1) Identifying the differentially expressed genes by meta-analysis, in bulk. 
  
  deg <- batchMeta(data.list = sample.sets, 
                   cutoff = 0.5, 
                   g.start = 1, 
                   g.end = num.gene)
  
  deg.meta <- c(rownames(eset)[deg$UpGene], 
                rownames(eset)[deg$DownGene])
  
  # 2) Identifying the differentially expressed genes by SVM-RFE, one by one. 
  
  # -- . 
  
  sam.lab <- group
  
  names(sam.lab) <- NULL
  
  eset.mat <- as.data.frame(t(eset))
  
  input <- cbind(sam.lab, eset.mat)
  
  input <- as.data.frame(input)
  
  # ranked.feat <- svmRFE(input, k = 5, halve.above = 1000)
  
  nfold <- 10
  nrows <- nrow(input)
  folds <- rep(1:nfold, len=nrows)[sample(nrows)]
  folds
  folds <- lapply(1:nfold, function(x) which(folds == x))
  
  results <- lapply(folds, svmRFE.wrap, input, k = 5, halve.above = 1000)
  
  
  # Obtain top features across ALL folds
  
  top.features <- WriteFeatures(results, input, save = FALSE)
  
  
  # Extracting top 500 genes. 
  
  deg.svm <- top.features$FeatureName[1:500]
  
  deg.svm <- as.character(deg.svm)
  
  print("----------------------------------------------------------------------")
  
  # 3) Further identifying the DEGs by permutation test, based on 500 genes from deg.svm. 
  
  all.genes <- deg.svm
  
  gene.count <- length(all.genes)
  
  deg.count <- NULL
  
  i <- 0
  
  pb <- txtProgressBar(min = 0, max = gene.count, style = 3, char = "+")
  
  for (g in all.genes) {
    
    i <- i + 1
    
    setTxtProgressBar(pb, i)
    
    # Display the progress bar!
    
    tmp <- Summarize(get(g) ~ sam.lab, data = input, digits = 3)
    
    # print(tmp)
    # boxplot(get(g) ~ sam.lab, data = input)
    
    deg.per <- try(independence_test(get(g) ~ sam.lab, 
                                     data = input), 
                   silent = FALSE)
    
    # deg.Z <- deg.per@statistic@teststatistic
    deg.p <- deg.per@distribution@pvalue(deg.per@statistic@teststatistic)
    
    flush.console()
    
    # This variable, deg.count, stores all the differentially expressed genes.
    
    if (!is.na(deg.p) & deg.p < 0.01) deg.count <- c(deg.count, g) else next
    
  }
  
  close(pb)
  
  deg.per <- deg.count; rm(deg.count)
  
  print("per")
  
  deg.int <- intersect(intersect(deg.meta, 
                                 deg.svm), 
                       deg.per)
  
  
  
  sam.iris <- input[, c("sam.lab", deg.int[1:10])]
  
  X <- sam.iris[, names(sam.iris) != "sam.lab"]
  
  y <- as.character(sam.iris$sam.lab)
  
  # please select the fold number, for k-fold CV. 
  
  folds <- 10
  
  test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning
  
  all.pred.tables <- lapply(1:folds, function(i) {
    
    test.id <- test.fold[[i]]
    
    X.train <- X[-test.id, ]
    
    y.train <- as.factor(y[-test.id])
    
    model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 5) # some tuning may be needed
    
    predict.test <- predict(model, X[test.id, ], prob = TRUE)
    
    prob.benign <- attr(predict.test, "probabilities")[, 2]
    
    data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
    
  })
  
  full.pred.table <- do.call(rbind, all.pred.tables)
  
  res.roc <- roc(full.pred.table$y.test, 
                 full.pred.table$y.pred, 
                 plot = TRUE, 
                 legacy.axes = TRUE)
  
  MSP.auc.value <- auc(res.roc)
  
  print("MSP")
  ################
  
  ##########
  
  
  label<-gsub("Control", 2, group)
  label<-gsub("Experimental", 1, label)
  
  samfit<-SAM(eset,label,resp.type="Two class unpaired")
  
  genes.up<-samfit$siggenes.table$genes.up
  genes.lo<-samfit$siggenes.table$genes.lo
  sam1<-as.data.frame(genes.up)
  sam2<-as.data.frame(genes.lo)
  
  deg.sam <-rbind(sam1,sam2)
  
  FC<-as.vector(deg.sam$`Fold Change`)
  
  qvalue<-as.vector(deg.sam$`q-value(%)`)
  
  deg.sam <- deg.sam[abs(as.numeric(FC)) > 1 & as.numeric(qvalue) < 0.05, ]
  
  deg.sam <- deg.sam[rev(order(abs(as.numeric(as.vector(deg.sam$`Fold Change`))))), ]
  
  deg.sam <- deg.sam$`Gene Name`
  
  deg.sam<-rownames(eset)[deg.sam]
  
  print("sam")
  
  sam.iris <- input[, c("sam.lab", deg.sam)]
  
  X <- sam.iris[, names(sam.iris) != "sam.lab"]
  
  y <- as.character(sam.iris$sam.lab)
  
  # please select the fold number, for k-fold CV. 
  
  folds <- 5
  
  test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning
  
  all.pred.tables <- lapply(1:folds, function(i) {
    
    test.id <- test.fold[[i]]
    
    X.train <- X[-test.id, ]
    
    y.train <- as.factor(y[-test.id])
    
    model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 5) # some tuning may be needed
    
    predict.test <- predict(model, X[test.id, ], prob = TRUE)
    
    prob.benign <- attr(predict.test, "probabilities")[, 2]
    
    data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
    
  })
  
  full.pred.table <- do.call(rbind, all.pred.tables)
  
  res.roc <- roc(full.pred.table$y.test, 
                 full.pred.table$y.pred, 
                 plot = TRUE, 
                 legacy.axes = TRUE)
  sam.auc.value <- auc(res.roc)
  
  #########################3
  colnames(eset)<-group
  
  all.tt<-matrix(1,2,4)
  
  colnames(all.tt)<-c("","FC","Pvalue","fdr")
  
  Pvalue<-c(rep(0,nrow(eset)))
  
  FC<-c(rep(0,nrow(eset)))
  
  for(ttt in 1:nrow(eset)){
    
    if(sd(eset[ttt,which(colnames(eset)=="Control")])==0&&eset[ttt,which(colnames(eset)=="Experimental")]==0){
      
      Pvalue <- "NA"
      
      FC<- "NA"
      
    }else{
      
      y=t.test(as.numeric(eset[ttt,which(colnames(eset)=="Control")]),as.numeric(eset[ttt,which(colnames(eset)=="Experimental")]))
      
      Pvalue<-y$p.value
      
      FC<-mean(as.numeric(eset[ttt,which(colnames(eset)=="Experimental")]))/mean(as.numeric(eset[ttt,which(colnames(eset)=="Control")]))
      
      FC<-log2(FC)
      
      fdr=p.adjust(Pvalue, "BH")
      
      p.tt<-cbind(rownames(eset)[ttt],FC,Pvalue,fdr)
      all.tt<-rbind(all.tt,p.tt)
      all.tt<-as.data.frame(all.tt)
      
    }
    
  }
  
  deg.tt <- all.tt[as.numeric(as.vector(all.tt$fdr))< 0.05, ]
  
  rownames(deg.tt)<-deg.tt[,1]
  deg.tt<-deg.tt[,-1]
  
  head(deg.tt)
  deg.tt<- deg.tt[order(deg.tt$fdr), ]
  
  deg.t <- rownames(deg.tt)
  
  print("t-test")
  
  sam.iris <- input[, c("sam.lab", deg.t[1:10])]
  
  X <- sam.iris[, names(sam.iris) != "sam.lab"]
  
  y <- as.character(sam.iris$sam.lab)
  
  # please select the fold number, for k-fold CV. 
  
  folds <- 5
  
  test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning
  
  all.pred.tables <- lapply(1:folds, function(i) {
    
    test.id <- test.fold[[i]]
    
    X.train <- X[-test.id, ]
    
    y.train <- as.factor(y[-test.id])
    
    model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 5) # some tuning may be needed
    
    predict.test <- predict(model, X[test.id, ], prob = TRUE)
    
    prob.benign <- attr(predict.test, "probabilities")[, 2]
    
    data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
    
  })
  
  full.pred.table <- do.call(rbind, all.pred.tables)
  
  res.roc <- roc(full.pred.table$y.test, 
                 full.pred.table$y.pred, 
                 plot = TRUE, 
                 legacy.axes = TRUE)
  t.auc.value <- auc(res.roc)
  ########################
  
  deg.Ft <- all.tt[abs(as.numeric(as.vector(all.tt$FC))) > log2(1.2) &as.numeric(as.vector(all.tt$fdr))< 0.05, ]
  
  deg.Ft<- deg.Ft[rev(order(abs(as.numeric(as.vector(deg.Ft$FC))))), ]
  
  rownames(deg.Ft) <- deg.Ft[,1]
  
  deg.tf<-rownames(deg.Ft)
  print("tf")
  set.seed(5)
  
  sam.iris <- input[, c("sam.lab", deg.tf[1:10])]
  
  X <- sam.iris[, names(sam.iris) != "sam.lab"]
  
  y <- as.character(sam.iris$sam.lab)
  
  # please select the fold number, for k-fold CV. 
  
  folds <- 5
  
  test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning
  
  all.pred.tables <- lapply(1:folds, function(i) {
    
    test.id <- test.fold[[i]]
    
    X.train <- X[-test.id, ]
    
    y.train <- as.factor(y[-test.id])
    
    model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 5) # some tuning may be needed
    
    predict.test <- predict(model, X[test.id, ], prob = TRUE)
    
    prob.benign <- attr(predict.test, "probabilities")[, 2]
    
    data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
    
  })
  
  full.pred.table <- do.call(rbind, all.pred.tables)
  
  res.roc <- roc(full.pred.table$y.test, 
                 full.pred.table$y.pred, 
                 plot = TRUE, 
                 legacy.axes = TRUE)
  
  
  tf.auc.value<- auc(res.roc)
  
  
  
  #################
  
  
  design <- model.matrix( ~ group)
  
  fit <- lmFit(eset, design)
  
  fit2 <- eBayes(fit, trend = FALSE)  
  
  limmaDEGs <- topTable(fit2, coef = 2, number = Inf)
  
  deg.limma <- limmaDEGs[abs(limmaDEGs$logFC) > 1 & limmaDEGs$adj.P.Val < 0.05, ]
  
  deg.limma <- deg.limma[rev(order(abs(deg.limma$logFC))), ]
  
  deg.limma <- rownames(deg.limma)
  
  print("limma")
  
  sam.iris <- input[, c("sam.lab", deg.limma[1:10])]
  
  X <- sam.iris[, names(sam.iris) != "sam.lab"]
  
  y <- as.character(sam.iris$sam.lab)
  
  # please select the fold number, for k-fold CV. 
  
  folds <- 5
  
  test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning
  
  all.pred.tables <- lapply(1:folds, function(i) {
    
    test.id <- test.fold[[i]]
    
    X.train <- X[-test.id, ]
    
    y.train <- as.factor(y[-test.id])
    
    model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 5) # some tuning may be needed
    
    predict.test <- predict(model, X[test.id, ], prob = TRUE)
    
    prob.benign <- attr(predict.test, "probabilities")[, 2]
    
    data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
    
  })
  
  full.pred.table <- do.call(rbind, all.pred.tables)
  
  res.roc <- roc(full.pred.table$y.test, 
                 full.pred.table$y.pred, 
                 plot = TRUE, 
                 legacy.axes = TRUE)
  limma.auc.value <- auc(res.roc)
  ###############
  
  
  mcr.up<-paste0("g",up.id)
  mcr.down<-paste0("g",down.id)
  
  sam.ACC.up<-intersect(mcr.up,deg.sam)
  sam.ACC.down<-intersect(mcr.down,deg.sam)
  sam.ACC<-c(sam.ACC.up,sam.ACC.down)
  
  t.ACC.up<-intersect(mcr.up,deg.t)
  t.ACC.down<-intersect(mcr.down,deg.t)
  t.ACC<-c(t.ACC.up,t.ACC.down)
  
  tf.ACC.up<-intersect(mcr.up,deg.tf)
  tf.ACC.down<-intersect(mcr.down,deg.tf)
  tf.ACC<-c(tf.ACC.up,tf.ACC.down)
  
  MSP.ACC.up<-intersect(mcr.up,deg.int)
  MSP.ACC.down<-intersect(mcr.down,deg.int)
  MSP.ACC<-c(MSP.ACC.up,MSP.ACC.down)
  
  
  limma.ACC.up<-intersect(mcr.up,deg.limma)
  limma.ACC.down<-intersect(mcr.down,deg.limma)
  limma.ACC<-c(limma.ACC.up,limma.ACC.down)
  
  limma.red<-c(length(deg.limma), limma.auc.value,length(limma.ACC))
  MSP.red<-c(length(deg.int),MSP.auc.value, length(MSP.ACC))
  sam.red<-c(length(deg.sam), sam.auc.value,length(sam.ACC))
  t.red<-c(length(deg.t),t.auc.value, length(t.ACC))
  tf.red<-c(length(deg.tf), tf.auc.value,length(tf.ACC))
  
  simulation.num<-length(up.id)+length(down.id)
  all.red<-c(simulation.num,MSP.red, limma.red, sam.red, t.red, tf.red)
  
  mm<-t(as.data.frame(all.red))
  
  mm<-as.data.frame(mm)
  
  rownames(mm)<-paste0(num1[num],":",num2[num])
  
  colnames(mm)<-c("simulation.num","MSP.num","MSP.AUC","MSP.ACC",
                  "limma.num","limma.AUC","limma.ACC",
                  "SAM.num","SAM.AUC","SAM.ACC",
                  "t-test.num","t-test.AUC","t-test.ACC",
                  "t-test & FC.num","t-test & FC.AUC","t-test & FC.ACC")
  
  red.ACC<-rbind(red.ACC,mm)
}













int <- c("g1", "g224", "g288", "g377", "g52")

limma <- c("g423", "g35", "g129", "g137", "g449")

for (a in 3:200) {
  
  library(pROC)
  
  sam.iris <- input[, c("sam.lab", deg.deseq2[2:a])]
  
  X <- sam.iris[, names(sam.iris) != "sam.lab"]
  
  y <- as.character(sam.iris$sam.lab)
  
  # please select the fold number, for k-fold CV. 
  
  folds <- 5
  
  test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning
  
  all.pred.tables <- lapply(1:folds, function(i) {
    
    test.id <- test.fold[[i]]
    
    X.train <- X[-test.id, ]
    
    y.train <- as.factor(y[-test.id])
    
    model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 5) # some tuning may be needed
    
    predict.test <- predict(model, X[test.id, ], prob = TRUE)
    
    prob.benign <- attr(predict.test, "probabilities")[, 2]
    
    data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
    
  })
  
  full.pred.table <- do.call(rbind, all.pred.tables)
  
  res.roc <- roc(full.pred.table$y.test, 
                 full.pred.table$y.pred, 
                 plot = TRUE, 
                 legacy.axes = TRUE)
  
  roc.df <- data.frame(TPP = res.roc$sensitivities * 100, 
                       FPP = (1 - res.roc$specificities) * 100, 
                       Thresholds = res.roc$thresholds)
  
  # windowsFonts(A = windowsFont("Times New Roman"))
  
  plot.roc(res.roc, 
           col = "#377EB8", 
           legacy.axes = TRUE, 
           percent = TRUE, 
           xlab = "False Positive Percentage", 
           ylab = "TRUE Positive Percentage", 
           lwd = 4, 
           print.auc = TRUE, 
           print.auc.x = 0.45, 
           print.auc.y = 0.45,
           auc.polygon = TRUE, 
           auc.polygon.col = "#377EB822") 
  
  plot.roc(res.roc, 
           col = "#377EB8", 
           legacy.axes = TRUE, 
           percent = TRUE, 
           xlab = "False Positive Percentage", 
           ylab = "TRUE Positive Percentage", 
           lwd = 4, 
           print.auc = TRUE, 
           print.auc.x = 0.85, 
           print.auc.y = 0.85,
           auc.polygon = TRUE, 
           auc.polygon.col = "#377EB822", 
           add = TRUE) 
  
  legend("bottomright", 
         legend = c("SVM-RFE", "LIMMA", "edgeR", "DESeq2"), 
         col = mypal2[1:4], 
         lwd = 4)
  
  par(pty = "m")
  
  # lines(perf, col = "green")
  
  auc.value <- auc(res.roc)
  
  print(auc.value)
  
  Sys.sleep(3)
  
}


f <- function(x) x[2] / (x[1] + x[2])

f(table(gsub("g", "", deg.int) < 500))

f(table(gsub("g", "", deg.limma) < 500))

f(table(gsub("g", "", deg.deseq2) < 500))

f(table(gsub("g", "", deg.edgeR) < 500))




data <- iris[1:100, ]

train.id <- sample(1:100, 50)

train.set <- data[train.id, ]

test.set <- data[-train.id, -5]

glm.out <- glm(Species ~ Sepal.Width + Sepal.Length + Petal.Width + Petal.Length,
               data = train.set,
               family = binomial) # family = binomial required for logistic regression

summary(glm.out)

plant1 <- data.frame(Sepal.Length=6.4, Sepal.Width=2.8, Petal.Length=4.6, Petal.Width=1.8)

predict(glm.out, test.set[, -5], type="response")




