

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


f <- function(x) x[1] / x[2]

f(table(gsub("g", "", deg.int) < 500))

f(table(gsub("g", "", deg.limma) < 500))

f(table(gsub("g", "", deg.deseq2) < 500))

f(table(gsub("g", "", deg.edgeR) < 500))



