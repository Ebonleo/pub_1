tb.cln=read.delim("clipboard",stringsAsFactors = F,header = T)#copy clinical table and run this
tb.seq=read.delim("clipboard",stringsAsFactors = F,header = T,check.names = F)#copy RNA-seq table and run this
tb.srv=read.delim("clipboard",stringsAsFactors = F,header = T)#copy survival table and run this
tb.mrg=merge(tb.srv,tb.cln,by.x = "sample",by.y = "sampleID",all = F)
smp_name=tb.mrg$sample
seq_name=colnames(tb.seq)
p_name=intersect(smp_name,seq_name)
tb.msq=tb.seq[,p_name]
rownames(tb.mrg)=tb.mrg$sample
tb.wsh=tb.mrg[p_name,]
rownames(tb.seq)=tb.seq$sample
tb.seq_wsh=tb.seq[gn_name,]
tb.seq_wsh=tb.seq_wsh[,-1]
tb.seq_wsh=t(tb.seq_wsh)
tb.seq_wsh=as.data.frame(tb.seq_wsh)
tb.seq_wsh$name=rownames(tb.seq_wsh)
tb.wsh=merge(tb.wsh,tb.seq_wsh,by.x = "sample",by.y = "name")#the table for analysis
tb.wsh$random_number=sample(nrow(tb.wsh),replace = F)
tb.trn=subset(tb.wsh,random_number>=nrow(tb.wsh)/3)
tb.tst=subset(tb.wsh,random_number<nrow(tb.wsh)/3)
library(survival)
library(survivalROC)
library(survminer)
fml_multicox=as.formula(paste0("Surv(OS.time,OS)~",paste(gn_name,collapse = "+")))
res_cox=coxph(fml_multicox,data = tb.trn)

attach(tb.trn)
tb.trn$risk_score=exp( (0.023211 * BIRC5) + (-0.360271 * FYN) + (-0.059331 * IGF1) + (-0.062212 * MASP1) + (-0.002987 * NR3C2) + (-0.057015 * TGFBR3))
detach(tb.trn)

roc_trn=survivalROC(Stime = tb.trn$OS.time,status = tb.trn$OS,marker = tb.trn$risk_score,predict.time = 365*5,method = "KM")
tb.roc=data.frame(roc_trn$TP,roc_trn$FP)
ggplot(tb.roc,aes(x = roc_trn.FP, y = roc_trn.TP,label = 365*5))+geom_roc(stat = "identity",labels = F)+
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_pubr()+labs(x="FP\nROC of Training Set (AUC=0.7)",y="TP")

which(roc_trn$TP-roc_trn$FP==max(roc_trn$TP-roc_trn$FP))
roc_trn$cut.values[168]

tb.trn$group=tb.trn$risk_score<=0.01222419
tb.trn$OSY=tb.trn$OS.time/365
svf_trn=survfit(Surv(OSY,OS)~group,data = tb.trn)
p_srv_trn=ggsurvplot(svf_trn,pval = T,pval.method = T,surv.median.line = "hv",legend.title="Group",legend.labs=c("High Risk","Low Risk"),conf.int = T)+labs(x="Time (Years)")

attach(tb.tst)
tb.tst$risk_score=exp( (0.023211 * BIRC5) + (-0.360271 * FYN) + (-0.059331 * IGF1) + (-0.062212 * MASP1) + (-0.002987 * NR3C2) + (-0.057015 * TGFBR3))
detach(tb.tst)

roc_tst=survivalROC(Stime = tb.tst$OS.time,status = tb.tst$OS,marker = tb.tst$risk_score,predict.time = 365*5,method = "KM")
tb.roc2=data.frame(roc_tst$TP,roc_tst$FP)
ggplot(tb.roc2,aes(x = roc_tst.FP, y = roc_tst.TP,label = 365*5))+geom_roc(stat = "identity",labels = F)+
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_pubr()+labs(x="FP\nROC of Training Set (AUC=0.7)",y="TP")

which(roc_tst$TP-roc_tst$FP==max(roc_tst$TP-roc_tst$FP))
roc_trn$cut.values[34]
tb.tst$group=tb.tst$risk_score<=0.005782399
svf_tst=survfit(Surv(OSY,OS)~group,data = tb.tst)
p_srv_tst=ggsurvplot(svf_tst,pval = T,pval.method = T,surv.median.line = "hv",legend.title="Group",legend.labs=c("High Risk","Low Risk"),conf.int = T)+labs(x="Time (Years)")

attach(tb.wsh)
tb.wsh$risk_score=exp( (0.023211 * BIRC5) + (-0.360271 * FYN) + (-0.059331 * IGF1) + (-0.062212 * MASP1) + (-0.002987 * NR3C2) + (-0.057015 * TGFBR3))
detach(tb.wsh)
roc_wsh=survivalROC(Stime = tb.wsh$OS.time,status = tb.wsh$OS,marker = tb.wsh$risk_score,predict.time = 365*5,method = "KM")
tb.roc3=data.frame(roc_wsh$TP,roc_wsh$FP)
which(roc_wsh$TP-roc_wsh$FP==max(roc_wsh$TP-roc_wsh$FP))
roc_wsh$cut.values[166]

tb.wsh$group=tb.wsh$risk_score<=0.009170362
tb.wsh$OSY=tb.wsh$OS.time/365
svf_wsh=survfit(Surv(OSY,OS)~group,data = tb.wsh)
p_srv_wsh=ggsurvplot(svf_wsh,pval = T,pval.method = T,surv.median.line = "hv",legend.title="Group",legend.labs=c("High Risk","Low Risk"),conf.int = T)+labs(x="Time (Years)")


tb.trn=tb.trn[order(tb.trn$risk_score),]
tb.trn$order_number=1:nrow(tb.trn)/nrow(tb.trn)
ggplot(tb.trn,aes(x=risk_score,y=order_number,color=group))+geom_point()+theme_pubr()+labs(x="Risk Score",y="",color="Group")+scale_colour_discrete(labels=c("High Risk","Low Risk"))->p_ctr_trn

tb.tst=tb.tst[order(tb.tst$risk_score),]
tb.tst$order_number=1:nrow(tb.tst)/nrow(tb.tst)
ggplot(tb.tst,aes(x=risk_score,y=order_number,color=group))+geom_point()+theme_pubr()+labs(x="Risk Score",y="",color="Group")+scale_colour_discrete(labels=c("High Risk","Low Risk"))->p_ctr_tst

tb.wsh=tb.wsh[order(tb.wsh$risk_score),]
tb.wsh$order_number=1:nrow(tb.wsh)/nrow(tb.wsh)
ggplot(tb.wsh,aes(x=risk_score,y=order_number,color=group))+geom_point()+theme_pubr()+labs(x="Risk Score",y="",color="Group")+scale_colour_discrete(labels=c("High Risk","Low Risk"))->p_ctr_wsh

library(pheatmap)
tb.trn$set="TRS"
tb.tst$set="TES"
tb.ttl=rbind(tb.trn,tb.tst)
ann_row=data.frame(Group=tb.ttl$group,Set=factor(c(rep("TRS",248),rep("TES",123)),levels = c("TRS","TES")))
ann_row$Group[which(ann_row$Group==T)]="L"
ann_row$Group[which(ann_row$Group==F)]="H"
rownames(ann_row)=rownames(tb.ttl)
pheatmap(tb.ttl[,gn_name],cluster_rows = F,cluster_cols = F,show_rownames = F,annotation_row = ann_row,scale = "row")->p_heat_ttl
library(ggplotify)
p_heat_ttl=as.ggplot(p_heat_ttl)
write.csv(tb.ttl,"ttl.csv")

tb.dwn=read.delim("clipboard",check.names = F,stringsAsFactors = F,row.names = 1)
tb.dwn$group=tb.dwn$risk_score<=0.009170362

tb.dwn$log2=scale(tb.dwn$risk_score,center = T,scale = T)
ggplot(tb.dwn[which(tb.dwn$person_neoplasm_cancer_status!=""),],aes(x=person_neoplasm_cancer_status,y=log2,fill=person_neoplasm_cancer_status))+
  geom_boxplot(width=0.2)+stat_compare_means(label = "p.format",method = "wilcox.test")+
  theme_pubr(legend = "none")+labs(fill="",x="",y="Scaled Risk Score")->p_tumor_burden

ggplot(tb.dwn[which(tb.dwn$pathologic_T!=""),],
       aes(x=pathologic_T,y=log2,fill=pathologic_T))+
  geom_boxplot(width=0.2)+stat_compare_means(label = "p.format",method = "wilcox.test")+
  theme_pubr(legend = "none")+labs(fill="",x="",y="Scaled Risk Score")->p_t_stage

ggplot(tb.dwn[which(tb.dwn$pathologic_M!=""&tb.dwn$pathologic_M!="MX"),],
       aes(x=pathologic_M,y=log2,fill=pathologic_M))+
  geom_boxplot(width=0.2)+stat_compare_means(label = "p.format",method = "wilcox.test")+
  theme_pubr(legend = "none")+labs(fill="",x="",y="Scaled Risk Score")->p_m_stage

ggplot(tb.dwn[which(tb.dwn$pathologic_N!=""&tb.dwn$pathologic_N!="NX"),],
       aes(x=pathologic_N,y=log2,fill=pathologic_N))+
  geom_boxplot(width=0.2)+stat_compare_means(label = "p.format",method = "wilcox.test")+
  theme_pubr(legend = "none")+labs(fill="",x="",y="Scaled Risk Score")->p_n_stage

ggplot(tb.dwn[which(tb.dwn$pathologic_stage!=""&tb.dwn$pathologic_stage!="NX"),],
       aes(x=pathologic_stage,y=log2,fill=pathologic_stage))+
  geom_boxplot(width=0.2)+stat_compare_means(label = "p.format",method = "wilcox.test")+
  theme_pubr(legend = "none")+labs(fill="",x="",y="Scaled Risk Score")->p_stage

ggplot(tb.dwn[which(tb.dwn$neoplasm_histologic_grade!=""&tb.dwn$neoplasm_histologic_grade!="NX"),],
       aes(x=neoplasm_histologic_grade,y=log2,fill=neoplasm_histologic_grade))+
  geom_boxplot(width=0.2)+stat_compare_means(label = "p.format",method = "wilcox.test")+
  theme_pubr(legend = "none")+labs(fill="",x="",y="Scaled Risk Score")->p_g_stage

library(cowplot)

plot_grid(p_tumor_burden,p_t_stage,p_m_stage,p_n_stage,p_g_stage,p_stage,ncol = 2)

require(TCGAbiolinks)
require(maftools)

tb.mut=subset(LIHC_mutect2,Hugo_Symbol!="TTN"&Hugo_Symbol!="MUC16")
mut_name_l=tb.dwn$X_GENOMIC_ID_TCGA_LIHC_mutation_ucsc_maf_gene[which(tb.dwn$group==T)]
mut_name_h=tb.dwn$X_GENOMIC_ID_TCGA_LIHC_mutation_ucsc_maf_gene[which(tb.dwn$group==F)]

mut_l=tb.mut[which(tb.mut$Tumor_Sample_Barcode%in%mut_name_l),]
mut.h=tb.mut[which(tb.mut$Tumor_Sample_Barcode%in%mut_name_h),]
mut2_l=read.maf(mut_l)
mut2_h=read.maf(mut.h)

p_mut_l=as.ggplot(~plotmafSummary(maf = mut2_l, rmOutlier = TRUE, addStat = 'median', dashboard = T, titvRaw = FALSE))
p_mut_h=as.ggplot(~plotmafSummary(maf = mut2_h, rmOutlier = TRUE, addStat = 'median', dashboard = T, titvRaw = FALSE))

oncop_h=as.ggplot(~oncoplot(mut2_h,titleText = "High Risk Group",drawRowBar = F))
oncop_l=as.ggplot(~oncoplot(mut2_l,titleText = "Low Risk Group",drawRowBar = F))
