
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




