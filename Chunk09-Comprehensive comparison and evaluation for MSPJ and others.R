################################################################################
#    &&&....&&&    % Project: MSPJ approach for identification of DEGs         #
#  &&&&&&..&&&&&&  % Author: Bo Li, Huachun Yin, Jingxin Tao
#  &&&&&&&&&&&&&&  % Date: Mar. 1st, 2022                                      #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.5.3;                             #
#       &&&&       % Platform: x86_64-pc-linux-gnu (64-bit)                    #
#        &         %                                                           #
################################################################################
library(gtools)
library(pheatmap)
library(UpSetR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pROC)

### ------------------------------------------------------------------------ ###
### Step-01. The overlap of 10 methods

method.list <- list(deg.int,
                    deg.GA,
                    deg.limma,
                    deg.mt,
                    deg.RF,
                    deg.RP,
                    deg.sam,
                    deg.snr,
                    deg.t,
                    deg.wt)

names(method.list) <- c(
  "MSPJ",
  "GA",
  "limma",
  "multtest",
  "RF",
  "RankProd",
  "SAM",
  "SNR",
  "T.test",
  "Wilcoxon's test"
)

upset(
  fromList(method.list),
  nsets = length(method.list),
  point.size = 3,
  line.size = 1,
  number.angles = 0,
  order.by = "freq",
  text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1),
  # ytitle, ylabel, xtitle, xlabel, sets, number
  matrix.color = "black",
  main.bar.color = 'black',
  mainbar.y.label = 'Intersection Size',
)
### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. Jaccard score of different methods

b <- Jscaore.YHC(method.list)

#the heatmap of the Jaccard score
p <- pheatmap(
  b,
  show_rownames = T,
  show_colnames = T,
  cluster_cols = T,
  cluster_rows = T,
  fontsize_row = 9,
  fontsize_col = 9,
  border_color = "NA",
  scale = "none",
  angle_col = 90,
  color = colorRampPalette(c("#2C4192", "white", "#E21E22"))(100),
  clustering_distance_rows = 'euclidean',
  clustering_method = 'single',
)

### End of Step-02.
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-03. The Jaccard score of GO terms

GOlist <- list(deg.int,
               deg.GA,
               deg.limma,
               deg.mt,
               deg.RF,
               deg.RP,
               deg.sam,
               deg.snr,
               deg.t,
               deg.wt)

GOname <- list()

for (Ji in 1:length(GOlist)) {
  gene <- unlist(GOlist[Ji])
  
  #geneID <- bitr(gene, fromType = "SYMBOL",
  #               toType = c("ENSEMBL", "ENTREZID"),
  #               OrgDb = org.Hs.eg.db)
  #
  gene.GO <- enrichGO(
    gene = gene,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "none",
    pvalueCutoff = 0.05
  )
  
  GOresult <- gene.GO@result
  
  GOname[[Ji]] <- GOresult$ID
  
}

names(GOname) <- c(
  "MSPJ",
  "GA",
  "limma",
  "multtest",
  "RF",
  "RankProd",
  "SAM",
  "SNR",
  "T.test",
  "Wilcoxon's test"
)

b <- Jscaore.YHC(GOname)


#the heatmap of the Jaccard score
p <- pheatmap(
  b,
  show_rownames = T,
  show_colnames = T,
  cluster_cols = T,
  cluster_rows = T,
  fontsize_row = 9,
  fontsize_col = 9,
  border_color = "NA",
  scale = "none",
  angle_col = 90,
  color = colorRampPalette(c("#2C4192", "white", "#E21E22"))(100),
  clustering_distance_rows = 'euclidean',
  clustering_method = 'single',
)
### End of Step-03.
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-04. Constructing the SVM models and validation with nested k-fold CV.

deg.cer <- deg.int[1:10] ##### top 10 genes as the signatures

sam.iris <- input[, c("sam.lab", deg.cer)]

X <- sam.iris[, names(sam.iris) != "sam.lab"]

y <- as.character(sam.iris$sam.lab)

# please select the fold number, for k-fold CV.

folds <- 10

set.seed(5)

test.fold <- split(sample(1:length(y)), 1:folds) # ignore warning

all.pred.tables <- lapply(1:folds, function(i) {
  test.id <- test.fold[[i]]
  
  X.train <- X[-test.id,]
  
  y.train <- as.factor(y[-test.id])
  
  model <-
    svm(X.train,
        y.train,
        kernel = "radial",
        prob = TRUE,
        cross = 5) # some tuning may be needed
  
  predict.test <- predict(model, X[test.id,], prob = TRUE)
  
  prob.benign <- attr(predict.test, "probabilities")[, 2]
  
  data.frame(y.test = y[test.id], y.pred = prob.benign) # returning this
  
})

full.pred.table <- do.call(rbind, all.pred.tables)

res.roc <- roc(
  full.pred.table$y.test,
  full.pred.table$y.pred,
  plot = TRUE,
  legacy.axes = TRUE
)

roc.df <- data.frame(
  TPP = res.roc$sensitivities * 100,
  FPP = (1 - res.roc$specificities) * 100,
  Thresholds = res.roc$thresholds
)

# windowsFonts(A = windowsFont("Times New Roman"))

plot.roc(
  res.roc,
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
  auc.polygon.col = "#377EB822"
)

plot.roc(
  res.roc,
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
  add = TRUE
)

#######calculated the value of specificity,sensitivity and accuracy

rets <- c("threshold", "specificity", "sensitivity", "accuracy")

feedback <-
  try(pROC::ci.coords(res.roc,
                      x = "best",
                      input = "threshold",
                      ret = rets),
      silent = TRUE)

error_or <-
  "Error in enforce.best.policy(res, best.policy) : \n  More than one \"best\" threshold was found, aborting. Change 'best.policy' to alter this behavior.\n"

if (error_or %in% feedback == TRUE) {
  feedback1 <-
    try(pROC::ci.coords(res.roc,
                        x = "best",
                        input = "threshold",
                        ret = rets),
        silent = TRUE)
  
  if (error_or %in% feedback1 == TRUE) {
    thresholdx <- 0.9
    ciorder <-
      pROC::ci.coords(res.roc,
                      x = thresholdx,
                      input = "threshold",
                      ret = rets)
    
  } else{
    ciorder <- feedback1
  }
  
} else{
  ciorder <- feedback
  
}

order_ci <- c(as.character(a), method[met], pROC::ci(res.roc), ciorder)

order_ci <- t(as.matrix(unlist(order_ci)))

### End of Step-04.
### ------------------------------------------------------------------------ ###
