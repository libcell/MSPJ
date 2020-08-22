



setwd("/Users/touyasushishin/Desktop")


####................sep one .....prepare data...................................

GSE72056_melanoma <- read.table(gzfile("GSE72056_melanoma_single_cell_revised_v2.txt.gz"), 
                                header = TRUE)

##构建黑色素瘤的表达矩阵
expre <- GSE72056_melanoma[-c(1,2,3), ]
colnames(expre) <- colnames(GSE72056_melanoma)

expre <- expre[!duplicated(expre[,1]), ]
dim(expre)
rownames(expre) <- expre$Cell

#expre <- aggregate(.~Cell,expre ,mean)


##记录样本的信息
sample_info <- GSE72056_melanoma[c(1,2,3), ]
rownames(sample_info) <- c("tumor","mailgant","cellType")
sample_info <- sample_info[, -1]

#用字符标注
celltype <- as.character(sample_info[3, ])

sum(table(celltype))

celltype[celltype == "1"] <- "T-cell"
celltype[celltype == "2"] <- "B-cell"
celltype[celltype == "3"] <- "Macrophages"
celltype[celltype == "4"] <- "Endothelial"
celltype[celltype == "5"] <- "CAF"
celltype[celltype == "6"] <- "NK-cell"
table(celltype)
sum(table(celltype))


non_op <- which(celltype == "0")

for (i in non_op) {
  if(sample_info[2, i] == 2) {
    celltype[i] <- "malignant" 
  } else if (sample_info[2, i] == 0) {
    celltype[i] <- "unresolved" 
    } else {celltype[i] <- "unknow-non-malignant-cell" 
  }
}

table(celltype)
sum(table(celltype))

sample_info <- rbind(celltype, sample_info)


###..................end........................................................






###................step two .....sample.........................................

 
Tcell <- colnames(sample_info[, celltype == "T-cell"])
Bcell <- colnames(sample_info[, celltype == "B-cell"])
NKcell <- colnames(sample_info[, celltype == "NK-cell"])
Maccell <- colnames(sample_info[, celltype == "Macrophages"])

Endocell <- colnames(sample_info[, celltype == "Endothelial"])
CAFcell <- colnames(sample_info[, celltype == "CAF"])
UNknowcell <- colnames(sample_info[, celltype == "unknow-non-malignant-cell"])
TUNresolvedcell <- colnames(sample_info[, celltype == "unresolved"])

malcell <- colnames(sample_info[, celltype == "malignant" ])


immunecell <- c(Tcell, Bcell, NKcell, Maccell)

##先随机100的随机数
set.seed(12345) 
f <- rnorm(117, mean = 0.33, sd = 0.3)
##对随机数进行截取尾
f <- f[f <= 0.99 & f >= 0]


##确定抽样个数
num_cancer <- round(500*f)

num_immune <- round(500*(1-f))


cancerlist <-list()

info  <- NULL
for (i in 1:100) {
  #对于第一个样本
  cancername <- sample(malcell, num_cancer[i], replace = FALSE, prob = NULL)
  immunename <- sample(immunecell, num_immune[i], replace = FALSE, prob = NULL)
  samplename <- c(cancername, immunename)
  tmp <- expre[, samplename]
  cancerlist[[i]] <- apply(tmp, 1, mean)
  proportion <- sample_info[1, samplename]
  info <- c(info, proportion)
  
}

simulationcancer <- matrix(unlist(cancerlist),nrow=23684,ncol=100, dimnames=list(expre$Cell,paste0(rep("sample",100),1:100)))
bulk_data <- as.data.frame(simulationcancer)

cellproportion <- matrix(unlist(info),nrow=100,ncol=500, dimnames=list(paste0(rep("sample",100),1:100),paste0(rep("cell",500),5:100)))

save(cellproportion,file = "sampleinformation.Rdata" )
save(bulk_data,file = "bulk_data.Rdata")
###..................end........................................................






###.................deconbolution......................................................



library(immunedeconv)
##.............1.cibersort......................................................
sig_matrix <- read.table("LM22.txt",header=T,sep="\t",row.names=1,
                         check.names=F,skip = 0,comment.char = "")
source("microenvirionment/row-data/Cibersort.R")

res_cib   <- TaoCIBERSORT(sig_matrix,bulk_data, perm=10, QN=F)

##.............2.EPIC...........................................................
res_epic <-  deconvolute(bulk_data, "epic")
##............ 3.cquantiseq.....................................................
res_quantiseq <- deconvolute(bulk_data , "quantiseq", tumor = TRUE)

##............ 4.mcpcounter.....................................................
res_mcp_counter <- deconvolute(bulk_data , "mcp_counter")

##........... .5.timer..........................................................
indi <- rep("sample",time=100)
res_timer <- deconvolute(bulk_data, "timer", indications = indi)

##.............6.xcell..........................................................
res_xcell <- deconvolute(bulk_data , "xcell")

###...............END...........................................................





###...................评估结果..................................................
#########估计比值truefraction and true fraction
##真实抽样B细胞结果
estimatedcell <- apply(cellproportion, 1, table)
tureSingle <- NULL
for (i in 1:100) {
  tmp <-  (estimatedcell[[i]][which(names(estimatedcell[[i]])=="B-cell")])/500
  tureSingle <- c(tureSingle,tmp)
}
tureSingle<- as.numeric(tureSingle)


Bpp <- data.frame(
  tureSingle= as.numeric(tureSingle),
  ##cibersort B细胞结果
  cibBcell =as.numeric(apply(res_cib[,c(1,2)], 1, sum)),
  ##epic B细胞结果
  epicBcell = as.numeric(res_epic[1,-1]),
  ##mcpB细胞结果
  mcpBcell = as.numeric(res_mcp_counter[5,-1]),
  ##quantiseq细胞结果
  quantiseqBcell= as.numeric(res_quantiseq[1,-1]),
  ##xcell细胞结果
  xcellBcell = as.numeric(apply(res_xcell[c(2,24,27,32),-1], 2, sum))
)

rownames(Bpp) <-paste0(rep("sample",100),1:100)

class(Bpp$quantiseqBcell)

#计算相关性
library(Hmisc)#加载包
res2 <- cor(Bpp )
res2


###....................end......................................................

term <- c("genomics", "transcriptomics", "proteomics", "metabolomics")
pm <- pubmed_trend(term, year = 2000:2020)
plot(pm)

library(rentrez)


