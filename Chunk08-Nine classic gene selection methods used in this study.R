
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
### code chunk number 08: Nine gene selection methods used in this study.
### ****************************************************************************

# input: eset - a gene expression matrix, genes in lines and samples in columns. 
# output: degs - a genelist, including the gene names, pooled SMDs, and so on. 

c.p <- which(sam.lab == "Control")
e.p <- which(sam.lab == "Experimental")

pvalue<-0.05

### ------------------------------------------------------------------------ ###
### Step-01. LIMMA method for DNA microarray and RNA-seq data. 

design <- model.matrix( ~ sam.lab)

#y<-voom(eset, design) #RNA-seq 

#fit <- lmFit(y, design) #RNA-seq 

fit <- lmFit(eset, design)

fit2 <- eBayes(fit)  

limmaDEGs <- topTable(fit2, coef = 2, number = Inf)

deg.limma <- limmaDEGs[limmaDEGs$adj.P.Val< pvalue, ]

deg.limma <- deg.limma[order(deg.limma$adj.P.Val), ]

deg.limma <- rownames(deg.limma)

deg.limma<-na.omit(deg.limma)


### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. SAM method for DNA microarray and RNA-seq data. 

### ------------------------------------------------------------------------ ###
label<-gsub("Control",1,sam.lab)

label<-gsub("Experimental",2,label)

genenames <- rownames(eset)  # Get row names = gene names or IDs

data <- list(x = eset, y = label, geneid = genenames, genenames = genenames, logged2 = FALSE)

samr.obj <- samr(data, resp.type = "Two class unpaired", nperms = 100)
#assay.type=c("array","seq")
#samfit<-SAM(tpm1,label,resp.type="Two class unpaired")

#names(samr.obj)  # Get object names from samr.obj

delta.table <- samr.compute.delta.table(samr.obj)  # Compute thresholds for different deltas

delta.table1 <- as.data.frame(delta.table)

alldelta<-delta.table1$`median FDR`<0.1

if(length(which(alldelta==TRUE))==0){
  delta <-0
}else{
  delta <-delta.table[delta.table[, "median FDR"] < 0.1, ][1, ]  # Check delta corresponding to median FDR ~0.1
  
  delta <-delta[1]
  
  names(delta)<-NULL
}

#samr.plot(samr.obj, delta)

siggenes.table <- samr.compute.siggenes.table(samr.obj, delta, data, delta.table)  # Summarize significant genes

names(siggenes.table)  

genes.up<-siggenes.table$genes.up 

genes.lo<-siggenes.table$genes.lo

degsam <-rbind(genes.up,genes.lo)

degsam <-as.data.frame(degsam)

degsam$`q-value(%)`<-as.numeric(as.character(degsam$`q-value(%)`))

deg.sam <- degsam[degsam$`q-value(%)`<pvalue,]

deg.sam <- deg.sam$`Gene Name`

deg.sam <-as.character(deg.sam)

deg.sam<-na.omit(deg.sam)


### End of Step-02.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-03. T-test method

###############################################################

t.deg<-apply(eset, 1, function(x) {t.test(x[c.p], x[e.p])})

degt<-matrix(NA,1)

num.t<-length(names(t.deg))

for (ii in 1:num.t) {
  
  t<-t.deg[[ii]]$p.value
  
  names(t)<-names(t.deg)[ii]
  
  tmp.t<-as.matrix(t)
  
  degt<-rbind(degt,tmp.t)
  
}

colnames(degt)<-c("P.value")

degt<-degt[-1,]

fdr=p.adjust(degt, "BH") #adjust the p value by BH

deg.t<-as.data.frame(fdr)

FC <- apply(etpm, 1, function(x) {mean(x[e.p])/mean(x[c.p])})

FC<-as.data.frame(log2(FC))

deg.FC<-as.data.frame(FC)

degtt<-cbind(deg.FC,degt)

colnames(degtt)<-c("logFC","P.value")

deg.tt <-degtt[degtt$P.value<pvalue,]

deg.tt<- deg.tt[order(deg.tt$P.value), ]

deg.t <- rownames(deg.tt)

deg.t<-na.omit(deg.t)

#################################################################

### End of Step-03.
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-04. Wilcoxon's test

###############################################################

wt.deg<-apply(eset, 1, function(x) {wilcox.test(x[c.p], x[e.p])})

degwt<-matrix(NA,1)

num.wt<-length(names(wt.deg))

for (ii in 1:num.wt) {
  
  wt<-wt.deg[[ii]]$p.value
  
  names(wt)<-names(wt.deg)[ii]
  
  tmp.wt<-as.matrix(wt)
  
  degwt<-rbind(degwt,tmp.wt)
  
}


colnames(degwt)<-c("P.value")

degwt<-degwt[-1,]

fdrw=p.adjust(degwt, "BH") #adjust p value 

degwt<-as.data.frame(fdrw)

degwt<-cbind(deg.FC,degwt)

colnames(degwt)<-c("logFC","P.value")

degwt <-degwt[degwt$P.value<pvalue,]

degwt<- degwt[order(degwt$P.value), ]

deg.wt <-rownames(degwt)

deg.wt<-na.omit(deg.wt)

### End of Step-04.
### ------------------------------------------------------------------------ ###



### ------------------------------------------------------------------------ ###
### Step-05. RankProd method

###############################################################

labels<-gsub("Control","0",sam.lab)
labels<-gsub("Experimental","1",labels)

RPDEGs=RP(eset,labels, num.perm=500,logged=TRUE)

#topGene(RPDEGs,num.gene=10,gene.names=rownames(eset))

RPDEGs1=topGene(RPDEGs,cutoff=pvalue,method="pfp",logged=TRUE,logbase=2,gene.names=rownames(tmp1))
#pfp,pval

upgene<-RPDEGs1$Table1

downgene<-RPDEGs1$Table2

degRP<-rbind(upgene,downgene)

degRP<-as.data.frame(degRP)

degRP <- degRP[order(degRP$`RP/Rsum`), ]

deg.RP <-rownames(eset)[degRP$gene.index]

deg.RP<-na.omit(deg.RP)

### End of Step-05.
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-06. multtest

###############################################################

label<-gsub("2",0,labels)

label<-as.numeric(label)

maxttt<-mt.maxT(eset,label)

degmt <- maxttt[as.numeric(as.vector(maxttt$adjp))< pvalue, ]

degmt<- degmt[order(degmt$adjp), ]

deg.mt <- rownames(degmt)

deg.mt<-na.omit(deg.mt)


### End of Step-06.
### ------------------------------------------------------------------------ ###



### ------------------------------------------------------------------------ ###
### Step-07. SNR method

###############################################################

SNRdeg<-apply(eset,1,FUN = snr)

SNRdeg<-as.data.frame(SNRdeg)

SNRdeg$value<-SNRdeg$SNRdeg

SNRdeg <- SNRdeg[rev(order(abs(SNRdeg$SNRdeg))), ]

deg.snr<-rownames(SNRdeg)[1:floor(nrow(eset)*0.3)]

deg.snr<-na.omit(deg.snr)


### End of Step-07.
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-08. GA method

###############################################################

filter_evaluator <- filterEvaluator('binaryConsistency')
#Generates a filter function to be used as an evaluator in the feature selection proccess. 

search_method <- searchAlgorithm('geneticAlgorithm',
                                 list(popSize =200,pcrossover = 0.8,
                                      pmutation = 0.1,maxiter = 10,
                                      run = 100))
#list()??????modify parameters
res.GA <- featureSelection(input, 'sam.lab', search_method, filter_evaluator)
#data_x,data_y

GA.features<-res.GA$bestFeatures

GADEG<-t(GA.features)

GADEG<-apply(GADEG, 1, sum)

GADEG<-names(GADEG)[rev(order(GADEG))]

deg.GA<-GADEG[1:floor(nrow(eset)*0.3)]

deg.GA<-na.omit(deg.GA)


### End of Step-08.
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-09. RF

###############################################################

input1<-input

colnames(input1)<-paste0("A",colnames(input1))

colnames(input1)[1]<-"sam.lab"

set.seed(222)

RF<-randomForest(sam.lab ~., data=input1, 
                 mtry = 2, ntree=1000, 
                 importance=TRUE, proximity=TRUE)

deg.RF<-as.data.frame(RF$importance)

deg.RF<-deg.RF[rev(order(deg.RF$MeanDecreaseAccuracy)),]

deg.RF<-rownames(deg.RF10)[1:floor(nrow(eset)*0.3)]

deg.RF<-gsub("A","",deg.RF)

deg.RF<-na.omit(deg.RF)

### End of Step-09.
### ------------------------------------------------------------------------ ###
