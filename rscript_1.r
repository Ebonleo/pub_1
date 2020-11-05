#import copied table. clinical information, suvival and RNA-seq, which are manually downloaded manually, are needed.
x=read.delim("clipboard",stringsAsFactors = F,check.names = F)

library(limma)
#grouping matrix
group <- factor(rep(c('tumor', 'normal'), each = x), levels = c('tumor', 'normal'))
#the x should be the real count of patients in each group
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(group)

#voom normalizatoin
norm <- voom(exprSet, design, plot = TRUE)

fit <- lmFit(norm, design, method = 'ls')#fitting, ?lmFit for details
contrast <- makeContrasts('treat-control', levels = design)#contrasting matrix

#Bayes fitting
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

#DEGs
diff_gene <- topTable(fit2, number = Inf, adjust.method = 'fdr')
head(diff_gene, 10)

#univariate Cox regression
library(survival)
library(survminer)
uni_cox=coxph(Surv(time,event)~variable,dataset)
#time=survival time, event=survival status,grouping=grouping factor,dataset=dataset that including time, event and variable

#LASSO
cvfit = cv.glmnet(x, y, family = "cox")
#x is a matrix that contains the variables,such as gene expression, and y is a matrix that contains survival time and status
coef(cvfit,s="lambda.1se")
#to observe included variables

#multivariate Cox regression
multicox=coxph(Surv(time,event)~x+y+z,dataset)#x, y, and z are variale names
multicox

#ROC
library(survivalROC)
my_roc=survivalROC(Stime = time,status = status,marker = marker,predict.time = time2,method = "KM")
#Stime=survival time, status= survival status, marker= predict marker, predict.time= time for prediction
plot(x=my_roc$FP,y=my_roc$TP,type = "l")
abline(0,1,col="gray",lty=2)

#Survival analysis
my_survfit=survfit(Surv(time,event)~variable,dataset)
ggsurvplot(my_survfit)
#?ggsurvplot for more details

#mutation analysis
my_maf = read.maf(maf =maf_file)
#maf_file is the maf file manually downloaded
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#Summary plot of the maf file
oncoplot(maf = var_maf, top = 400, writeMatrix=T,removeNonMutated = F,showTumorSampleBarcodes=T)
#Oncoplot of the maf file
