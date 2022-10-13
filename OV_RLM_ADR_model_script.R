#Lipid metabolism pathway selection in the development of drug resistance

library(GEOquery)
library(limma)
library(pheatmap)
library(vcd)

#Incorporate A2780 and W1
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism/GSE73935")
GSE73935_exp<-read.table("GSE73935_exp.txt",header=T)
GSE73935_note<-read.csv("GSE73935_note.csv",sep=",")

head(GSE73935_exp)
GSE73935_note

#check
colnames(GSE73935_note)<-c("class","cell","ID_REF")
all(is.na(GSE73935_exp)==F)
all(is.na(GSE73935_note)==F)

max(GSE73935_exp[,-1])
min(GSE73935_exp[,-1])
#No log processing is required

#Select cell resistant to paclitaxel: A2780
GSE73935_note[intersect(grep("Paclitaxel*",GSE73935_note$class),grep("A2780*",GSE73935_note$cell)),]
A2780_pac<-GSE73935_note[intersect(grep("Paclitaxel*",GSE73935_note$class),grep("A2780*",GSE73935_note$cell)),][,3]
A2780_pac
A2780_pac<-GSE73935_exp[,c("ID_REF",A2780_pac)]
head(A2780_pac)

GSE73935_note[intersect(grep("Control*",GSE73935_note$class),grep("A2780*",GSE73935_note$cell)),]
A2780_control<-GSE73935_note[intersect(grep("Control*",GSE73935_note$class),grep("A2780*",GSE73935_note$cell)),][,3]
A2780_control
A2780_control<-GSE73935_exp[,c("ID_REF",A2780_control)]
head(A2780_control)

#Select cell resistant to paclitaxel: W1
GSE73935_note[intersect(grep("Paclitaxel*",GSE73935_note$class),grep("W1*",GSE73935_note$cell)),]
W1_pac<-GSE73935_note[intersect(grep("Paclitaxel*",GSE73935_note$class),grep("W1*",GSE73935_note$cell)),][,3]
W1_pac
W1_pac<-GSE73935_exp[,c("ID_REF",W1_pac)]
head(W1_pac)

GSE73935_note[intersect(grep("Control*",GSE73935_note$class),grep("W1*",GSE73935_note$cell)),]
W1_control<-GSE73935_note[intersect(grep("Control*",GSE73935_note$class),grep("W1*",GSE73935_note$cell)),][,3]
W1_control
W1_control<-GSE73935_exp[,c("ID_REF",W1_control)]
head(W1_control)

#GSE73935 chip note
GPL13667<-read.csv("GPL13667-15572 .csv",sep=",")
head(GPL13667)

#target gene
#Lipid_metabolism：ACLY,ACC1,ACC2,FASN,SCD,CPT1A,CPT1B,CPT1C,CPT2

#ACLY
GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]

#ACC1->ACACA, ACC2->ACACB
GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]

#FASN
GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]

#SCD
GPL13667[grep("SCD",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
#SCD only
GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]

#CPT1A
GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]

#CPT1B
GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]

#CPT1C
GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]

#CPT2
GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]

#union
GPL13667_lipidmeta_gene_probe<-rbind(
  GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)],
  GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)],
  GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)],
  GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)],
  GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)],
  GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)],
  GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)],
  GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)],
  GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
)
GPL13667_lipidmeta_gene_probe
GPL13667_lipidmeta_gene_probe<-GPL13667_lipidmeta_gene_probe[,c(1,3)]
GPL13667_lipidmeta_gene_probe
colnames(GPL13667_lipidmeta_gene_probe)<-c("probe","gene")

#check data distribution
boxplot(cbind(A2780_pac[,-1],A2780_control[,-1]))
boxplot(cbind(W1_pac[,-1],W1_control[,-1]))

max(GSE73935_exp[,-1])
min(GSE73935_exp[,-1])
#No log processing is required

#A2780 limma
A2780_for_limma<-merge(A2780_pac,A2780_control,by="ID_REF")
head(A2780_for_limma)

rownames(A2780_for_limma)<-A2780_for_limma$ID_REF
A2780_for_limma$ID_REF<-NULL
colnames(A2780_for_limma)<-c("t1","t2","t3","t4","t5","t6","c1","c2","c3")
head(A2780_for_limma)

A2780_group<-c("t","t","t","t","t","t","c","c","c")
A2780_design<-model.matrix(~0+factor(A2780_group))
colnames(A2780_design)<-levels(factor(A2780_group))
rownames(A2780_design)<-colnames(A2780_for_limma)
A2780_design

A2780_contrasts_fit<-makeContrasts(t-c,levels = A2780_design)
A2780_fit<-lmFit(A2780_for_limma,A2780_design)
A2780_fit2<-contrasts.fit(A2780_fit,A2780_contrasts_fit)
A2780_fit2<-eBayes(A2780_fit2)
A2780_limma_result<-topTable(A2780_fit2,coef = 1,n=Inf,adjust="BH")
head(A2780_limma_result)

A2780_pac[A2780_pac$ID_REF=="11749260_a_at",]
A2780_control[A2780_control$ID_REF=="11749260_a_at",]


#select probe
#ACLY
GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1)],]
A2780_ACLY<-mean(A2780_limma_result[GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_ACLY

#ACACA
GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1)],]
A2780_ACACA<-mean(A2780_limma_result[GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_ACACA

#ACACB
GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1)],]
A2780_ACACB<-mean(A2780_limma_result[GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_ACACB

#FASN
GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1)],]
A2780_FASN<-mean(A2780_limma_result[GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_FASN

#SCD
GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1)],]
A2780_SCD<-mean(A2780_limma_result[GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_SCD

#CPT1A
GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1)],]
A2780_CPT1A<-mean(A2780_limma_result[GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_CPT1A

#CPT1B
GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1)],]
A2780_CPT1B<-mean(A2780_limma_result[GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_CPT1B

#CPT1C
GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1)],]
A2780_CPT1C<-mean(A2780_limma_result[GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_CPT1C

#CPT2
GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_limma_result[GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1)],]
A2780_CPT2<-mean(A2780_limma_result[GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1)],][,1])
A2780_CPT2

#union
A2780_lpmeta_result<-as.data.frame(c(A2780_ACLY,A2780_ACACA,A2780_ACACB,A2780_FASN,A2780_SCD,A2780_CPT1A,A2780_CPT1B,A2780_CPT1C,A2780_CPT2))
A2780_lpmeta_result
rownames(A2780_lpmeta_result)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")
colnames(A2780_lpmeta_result)<-"logFC"

GPL13667_lipidmeta_gene_probe
rownames(GPL13667_lipidmeta_gene_probe)<-GPL13667_lipidmeta_gene_probe[,1]
A2780_lpmeta_raw<-A2780_limma_result[GPL13667_lipidmeta_gene_probe$probe,]
A2780_lpmeta_raw<-merge(A2780_lpmeta_raw,GPL13667_lipidmeta_gene_probe,by="row.names")
A2780_lpmeta_raw<-A2780_lpmeta_raw[,c(8,2,9)]
A2780_lpmeta_raw


#W1 limma

W1_for_limma<-merge(W1_pac,W1_control,by="ID_REF")
head(W1_for_limma)
rownames(W1_for_limma)<-W1_for_limma$ID_REF
W1_for_limma$ID_REF<-NULL
colnames(W1_for_limma)<-c("t1","t2","t3","c1","c2","c3")
head(W1_for_limma)

W1_group<-c("t","t","t","c","c","c")
W1_design<-model.matrix(~0+factor(W1_group))
colnames(W1_design)<-levels(factor(W1_group))
rownames(W1_design)<-colnames(W1_for_limma)
W1_design

W1_contrasts_fit<-makeContrasts(t-c,levels = W1_design)
W1_fit<-lmFit(W1_for_limma,W1_design)
W1_fit2<-contrasts.fit(W1_fit,W1_contrasts_fit)
W1_fit2<-eBayes(W1_fit2)
W1_limma_result<-topTable(W1_fit2,coef = 1,n=Inf,adjust="BH")
head(W1_limma_result)

W1_pac[W1_pac$ID_REF=="11759048_at",]
W1_control[W1_control$ID_REF=="11759048_at",]

#select gene

#ACLY
GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1)],]
W1_ACLY<-mean(W1_limma_result[GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_ACLY

#ACACA
GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1)],]
W1_ACACA<-mean(W1_limma_result[GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_ACACA

#ACACB
GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1)],]
W1_ACACB<-mean(W1_limma_result[GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_ACACB

#FASN
GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1)],]
W1_FASN<-mean(W1_limma_result[GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_FASN

#SCD
GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1)],]
W1_SCD<-mean(W1_limma_result[GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_SCD

#CPT1A
GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1)],]
W1_CPT1A<-mean(W1_limma_result[GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_CPT1A

#CPT1B
GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1)],]
W1_CPT1B<-mean(W1_limma_result[GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_CPT1B

#CPT1C
GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1)],]
W1_CPT1C<-mean(W1_limma_result[GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_CPT1C

#CPT2
GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_limma_result[GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1)],]
W1_CPT2<-mean(W1_limma_result[GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1)],][,1])
W1_CPT2

#union
W1_lpmeta_result<-as.data.frame(c(W1_ACLY,W1_ACACA,W1_ACACB,W1_FASN,W1_SCD,W1_CPT1A,W1_CPT1B,W1_CPT1C,W1_CPT2))
W1_lpmeta_result
rownames(W1_lpmeta_result)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")
colnames(W1_lpmeta_result)<-"logFC"

GPL13667_lipidmeta_gene_probe
rownames(GPL13667_lipidmeta_gene_probe)<-GPL13667_lipidmeta_gene_probe[,1]
W1_lpmeta_raw<-W1_limma_result[GPL13667_lipidmeta_gene_probe$probe,]
W1_lpmeta_raw<-merge(W1_lpmeta_raw,GPL13667_lipidmeta_gene_probe,by="row.names")
W1_lpmeta_raw<-W1_lpmeta_raw[,c(8,2,9)]
W1_lpmeta_raw

#Incorporate OVCAR3 and TOV21G
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism/GSE172016")
GSE172016_exp<-read.csv("GSE172016_gene.count.annotated.csv",sep=",")
head(GSE172016_exp)
rownames(GSE172016_exp)<-GSE172016_exp[,1]
GSE172016_note<-GSE172016_exp[,c(14,15)]

#OVCAR3 voom
OVCAR3_for_voom<-GSE172016_exp[,c(5,6,7,2,3,4)]
head(OVCAR3_for_voom)

all(is.na(OVCAR3_for_voom)==F)
max(OVCAR3_for_voom)
min(OVCAR3_for_voom)

OVCAR3_group_list<-as.factor(c("t","t","t","unt","unt","unt"))
OVCAR3_group_list
OVCAR3_design<-model.matrix(~factor(OVCAR3_group_list))
colnames(OVCAR3_design)=levels(factor(OVCAR3_group_list))
rownames(OVCAR3_design)=colnames(OVCAR3_for_voom)
OVCAR3_design

OVCAR3_voom<-voom(OVCAR3_for_voom,OVCAR3_design,normalize="quantile")

#check data distribution
OVCAR3_voom_exprSet_new<-OVCAR3_voom$E
OVCAR3_voom_sample<-ncol(OVCAR3_for_voom)
if(OVCAR3_voom_sample>40) par(cex = 0.5)
cols <- rainbow(OVCAR3_voom_sample*1.2)
par(mfrow=c(1,1))
boxplot(OVCAR3_for_voom, col = cols,main="expression value",las=2)
boxplot(OVCAR3_voom_exprSet_new, col = cols,main="expression value",las=2)

#continue voom

OVCAR3_fit<-lmFit(OVCAR3_voom,OVCAR3_design)
OVCAR3_fit2<-eBayes(OVCAR3_fit)
OVCAR3_tempOutput<-topTable(OVCAR3_fit2, coef=2, n=Inf)
OVCAR3_voom_result<-na.omit(OVCAR3_tempOutput)
head(OVCAR3_voom_result)
OVCAR3_voom_result$logFC<-(OVCAR3_voom_result$logFC)*(-1)
head(OVCAR3_voom_result)

OVCAR3_for_voom["ENSG00000259207",]
log(OVCAR3_for_voom["ENSG00000259207",],base = 2)

#note
head(GSE172016_note)
OVCAR3_voom_result_note<-merge(OVCAR3_voom_result,GSE172016_note,by="row.names")
OVCAR3_voom_result_note<-OVCAR3_voom_result_note[,c("Row.names","logFC","hgnc_symbol")]
head(OVCAR3_voom_result_note)

#c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

OVCAR3_lpmeta_result<-rbind(
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="ACLY",],
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="ACACA",],
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="ACACB",],
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="FASN",],
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="SCD",],
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="CPT1A",],
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="CPT1B",],
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="CPT1C",],
  OVCAR3_voom_result_note[OVCAR3_voom_result_note$hgnc_symbol=="CPT2",]
)

rownames(OVCAR3_lpmeta_result)<-OVCAR3_lpmeta_result$hgnc_symbol
OVCAR3_lpmeta_result<-OVCAR3_lpmeta_result[,c("logFC"),drop=F]
OVCAR3_lpmeta_result

#TOV21G voom

TOV21G_for_voom<-GSE172016_exp[,c(11,12,13,8,9,10)]
head(TOV21G_for_voom)

all(is.na(TOV21G_for_voom)==F)
max(TOV21G_for_voom)
min(TOV21G_for_voom)

TOV21G_group_list<-as.factor(c("t","t","t","unt","unt","unt"))
TOV21G_group_list
TOV21G_design<-model.matrix(~factor(TOV21G_group_list))
colnames(TOV21G_design)=levels(factor(TOV21G_group_list))
rownames(TOV21G_design)=colnames(TOV21G_for_voom)
TOV21G_design

TOV21G_voom<-voom(TOV21G_for_voom,TOV21G_design,normalize="quantile")

TOV21G_voom_exprSet_new<-TOV21G_voom$E
TOV21G_voom_sample<-ncol(TOV21G_for_voom)
if(TOV21G_voom_sample>40) par(cex = 0.5)
cols <- rainbow(TOV21G_voom_sample*1.2)
par(mfrow=c(1,1))
boxplot(TOV21G_for_voom, col = cols,main="expression value",las=2)
boxplot(TOV21G_voom_exprSet_new, col = cols,main="expression value",las=2)

TOV21G_fit<-lmFit(TOV21G_voom,TOV21G_design)
TOV21G_fit2<-eBayes(TOV21G_fit)
TOV21G_tempOutput<-topTable(TOV21G_fit2, coef=2, n=Inf)
TOV21G_voom_result<-na.omit(TOV21G_tempOutput)
head(TOV21G_voom_result)
TOV21G_voom_result$logFC<-(TOV21G_voom_result$logFC)*(-1)
head(TOV21G_voom_result)

TOV21G_for_voom["ENSG00000227954",]
log(TOV21G_for_voom["ENSG00000227954",],base = 2)

#note
head(GSE172016_note)
TOV21G_voom_result_note<-merge(TOV21G_voom_result,GSE172016_note,by="row.names")
TOV21G_voom_result_note<-TOV21G_voom_result_note[,c("Row.names","logFC","hgnc_symbol")]
head(TOV21G_voom_result_note)

#c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

TOV21G_lpmeta_result<-rbind(
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="ACLY",],
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="ACACA",],
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="ACACB",],
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="FASN",],
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="SCD",],
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="CPT1A",],
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="CPT1B",],
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="CPT1C",],
  TOV21G_voom_result_note[TOV21G_voom_result_note$hgnc_symbol=="CPT2",]
)

rownames(TOV21G_lpmeta_result)<-TOV21G_lpmeta_result$hgnc_symbol
TOV21G_lpmeta_result<-TOV21G_lpmeta_result[,c("logFC"),drop=F]
TOV21G_lpmeta_result


#Incorporate SKOV3
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism/GSE60335")
GSE60335_exp<-read.csv("GSE60335_exp.csv",sep=",")
head(GSE60335_exp)

colnames(GSE60335_exp)<-c("ID_REF","unt1","unt2","t50_1","t50_2","t600_1","t600_2")
GSE60335_exp[is.na(GSE60335_exp)==T]<-0

max(GSE60335_exp[,-1])
min(GSE60335_exp[,-1])

#SKOV3 voom

SKOV3_for_voom<-GSE60335_exp
rownames(SKOV3_for_voom)<-SKOV3_for_voom$ID_REF
SKOV3_for_voom<-SKOV3_for_voom[,c(4,5,6,7,2,3)]
head(SKOV3_for_voom)

SKOV3_group_list<-as.factor(c("t","t","t","t","unt","unt"))
SKOV3_group_list
SKOV3_design<-model.matrix(~factor(SKOV3_group_list))
colnames(SKOV3_design)=levels(factor(SKOV3_group_list))
rownames(SKOV3_design)=colnames(SKOV3_for_voom)
SKOV3_design

SKOV3_voom<-voom(SKOV3_for_voom,SKOV3_design,normalize="quantile")

#check data distribution
SKOV3_voom_exprSet_new<-SKOV3_voom$E
SKOV3_voom_sample<-ncol(SKOV3_for_voom)
if(SKOV3_voom_sample>40) par(cex = 0.5)
cols <- rainbow(SKOV3_voom_sample*1.2)
par(mfrow=c(1,1))
boxplot(SKOV3_for_voom, col = cols,main="expression value",las=2)
boxplot(SKOV3_voom_exprSet_new, col = cols,main="expression value",las=2)

SKOV3_fit<-lmFit(SKOV3_voom,SKOV3_design)
SKOV3_fit2<-eBayes(SKOV3_fit)
SKOV3_tempOutput<-topTable(SKOV3_fit2, coef=2, n=Inf)
SKOV3_voom_result<-na.omit(SKOV3_tempOutput)
head(SKOV3_voom_result)
SKOV3_voom_result$logFC<-(SKOV3_voom_result$logFC)*(-1)
head(SKOV3_voom_result)

SKOV3_for_voom["PH_hs_0042649",]
log(SKOV3_for_voom["PH_hs_0042649",],base = 2)

#note
GPL13287<-read.csv("GPL13287-12092.csv",sep=",")
head(GPL13287)

#ACLY
GPL13287[grep("ACLY",GPL13287$Gene_symbol),][,c(1,8,11)]

#ACC1->ACACA, ACC2->ACACB
GPL13287[grep("ACACA",GPL13287$Gene_symbol),][,c(1,8,11)]
GPL13287[grep("ACACB",GPL13287$Gene_symbol),][,c(1,8,11)]

#FASN
GPL13287[grep("FASN",GPL13287$Gene_symbol),][,c(1,8,11)]

#SCD
GPL13287[grep("SCD",GPL13287$Gene_symbol),][,c(1,8,11)]
#SCD only
GPL13287[grep("^SCD$",GPL13287$Gene_symbol),][,c(1,8,11)]

#CPT
GPL13287[grep("CPT1A",GPL13287$Gene_symbol),][,c(1,8,11)]
GPL13287[grep("CPT1B",GPL13287$Gene_symbol),][,c(1,8,11)]
GPL13287[grep("CPT1C",GPL13287$Gene_symbol),][,c(1,8,11)]

GPL13287[grep("CPT2",GPL13287$Gene_symbol),][,c(1,8,11)]

GPL13287_lipidmeta_gene_probe<-rbind(
  GPL13287[grep("ACLY",GPL13287$Gene_symbol),][,c(1,8,11)],
  GPL13287[grep("ACACA",GPL13287$Gene_symbol),][,c(1,8,11)],
  GPL13287[grep("ACACB",GPL13287$Gene_symbol),][,c(1,8,11)],
  GPL13287[grep("FASN",GPL13287$Gene_symbol),][,c(1,8,11)],
  GPL13287[grep("^SCD$",GPL13287$Gene_symbol),][,c(1,8,11)],
  GPL13287[grep("CPT1A",GPL13287$Gene_symbol),][,c(1,8,11)],
  GPL13287[grep("CPT1B",GPL13287$Gene_symbol),][,c(1,8,11)],
  GPL13287[grep("CPT1C",GPL13287$Gene_symbol),][,c(1,8,11)],
  GPL13287[grep("CPT2",GPL13287$Gene_symbol),][,c(1,8,11)]
)
GPL13287_lipidmeta_gene_probe
GPL13287_lipidmeta_gene_probe<-GPL13287_lipidmeta_gene_probe[,c(1,2)]
GPL13287_lipidmeta_gene_probe
colnames(GPL13287_lipidmeta_gene_probe)<-c("probe","gene")

#select probe
#ACLY
GPL13287[grep("ACLY",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("ACLY",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_ACLY<-mean(SKOV3_voom_result[GPL13287[grep("ACLY",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_ACLY

#ACACA
GPL13287[grep("ACACA",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("ACACA",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_ACACA<-mean(SKOV3_voom_result[GPL13287[grep("ACACA",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_ACACA

#ACACB
GPL13287[grep("ACACB",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("ACACB",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_ACACB<-mean(SKOV3_voom_result[GPL13287[grep("ACACB",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_ACACB

#FASN
GPL13287[grep("FASN",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("FASN",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_FASN<-mean(SKOV3_voom_result[GPL13287[grep("FASN",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_FASN

#SCD
GPL13287[grep("^SCD$",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("^SCD$",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_SCD<-mean(SKOV3_voom_result[GPL13287[grep("^SCD$",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_SCD

#CPT1A
GPL13287[grep("CPT1A",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("CPT1A",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_CPT1A<-mean(SKOV3_voom_result[GPL13287[grep("CPT1A",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_CPT1A

#CPT1B
GPL13287[grep("CPT1B",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("CPT1B",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_CPT1B<-mean(SKOV3_voom_result[GPL13287[grep("CPT1B",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_CPT1B

#CPT1C
GPL13287[grep("CPT1C",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("CPT1C",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_CPT1C<-mean(SKOV3_voom_result[GPL13287[grep("CPT1C",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_CPT1C

#CPT2
GPL13287[grep("CPT2",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_voom_result[GPL13287[grep("CPT2",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_CPT2<-mean(SKOV3_voom_result[GPL13287[grep("CPT2",GPL13287$Gene_symbol),][,c(1)],][,1])
SKOV3_CPT2

#union
SKOV3_lpmeta_result<-as.data.frame(c(SKOV3_ACLY,SKOV3_ACACA,SKOV3_ACACB,SKOV3_FASN,SKOV3_SCD,SKOV3_CPT1A,SKOV3_CPT1B,SKOV3_CPT1C,SKOV3_CPT2))
SKOV3_lpmeta_result
rownames(SKOV3_lpmeta_result)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")
colnames(SKOV3_lpmeta_result)<-"logFC"

GPL13287_lipidmeta_gene_probe
rownames(GPL13287_lipidmeta_gene_probe)<-GPL13287_lipidmeta_gene_probe[,1]
SKOV3_lpmeta_raw<-SKOV3_voom_result[GPL13287_lipidmeta_gene_probe$probe,]
SKOV3_lpmeta_raw<-merge(SKOV3_lpmeta_raw,GPL13287_lipidmeta_gene_probe,by="row.names")
SKOV3_lpmeta_raw<-SKOV3_lpmeta_raw[,c(8,2,9)]
SKOV3_lpmeta_raw

#Incorporate OV90
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism/GSE26465")
GSE26465<-read.csv("GSE26465_non-normalized.csv",sep=",")
head(GSE26465)
GSE26465$Probe.Id<-paste("I",GSE26465$Probe.Id,sep="")

GSE26465_exp<-GSE26465[,c("Probe.Id","AVG_Signal.OV90.D6.1","AVG_Signal.OV90.D6.2","AVG_Signal.OV90.D7.1","AVG_Signal.OV90.D7.2","AVG_Signal.OV90.P3.1","AVG_Signal.OV90.P3.2","AVG_Signal.OV90.P7.1","AVG_Signal.OV90.P7.2","AVG_Signal.OV90.1","AVG_Signal.OV90.2")]
head(GSE26465_exp)
GSE26465_note<-GSE26465[,c("Probe.Id","Symbol")]
head(GSE26465_note)

max(GSE26465_exp[,-1])
min(GSE26465_exp[,-1])

#voom

OV90_for_voom<-GSE26465_exp
head(OV90_for_voom)
rownames(OV90_for_voom)<-OV90_for_voom$Probe.Id
OV90_for_voom$Probe.Id<-NULL

#OV90 voom

all(is.na(OV90_for_voom)==F)
max(OV90_for_voom)
min(OV90_for_voom)

OV90_group_list<-as.factor(c("t","t","t","t","t","t","t","t","unt","unt"))
OV90_group_list
OV90_design<-model.matrix(~factor(OV90_group_list))
colnames(OV90_design)=levels(factor(OV90_group_list))
rownames(OV90_design)=colnames(OV90_for_voom)
OV90_design
#
OV90_voom<-voom(OV90_for_voom,OV90_design,normalize="quantile")
OV90_voom_exprSet_new<-OV90_voom$E
OV90_voom_sample<-ncol(OV90_for_voom)
if(OV90_voom_sample>40) par(cex = 0.5)
cols <- rainbow(OV90_voom_sample*1.2)
par(mfrow=c(1,1))
boxplot(OV90_for_voom, col = cols,main="expression value",las=2)
boxplot(OV90_voom_exprSet_new, col = cols,main="expression value",las=2)

OV90_fit<-lmFit(OV90_voom,OV90_design)
OV90_fit2<-eBayes(OV90_fit)
OV90_tempOutput<-topTable(OV90_fit2, coef=2, n=Inf)
OV90_voom_result<-na.omit(OV90_tempOutput)
head(OV90_voom_result)
OV90_voom_result$logFC<-(OV90_voom_result$logFC)*(-1)
head(OV90_voom_result)

OV90_for_voom["I4590398",]
log(OV90_for_voom["I4590398",],base = 2)
OV90_voom_result["I4590398",]

#note
head(OV90_voom_result)
head(GSE26465_note)

#ACLY
GSE26465[grep("ACLY",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("ACLY",GSE26465$Symbol),][,c(2)],]
OV90_ACLY<-mean(OV90_voom_result[GSE26465[grep("ACLY",GSE26465$Symbol),][,c(2)],][,1])
OV90_ACLY

#ACACA
GSE26465[grep("ACACA",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("ACACA",GSE26465$Symbol),][,c(2)],]
OV90_ACACA<-mean(OV90_voom_result[GSE26465[grep("ACACA",GSE26465$Symbol),][,c(2)],][,1])
OV90_ACACA

#ACACB
GSE26465[grep("ACACB",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("ACACB",GSE26465$Symbol),][,c(2)],]
OV90_ACACB<-mean(OV90_voom_result[GSE26465[grep("ACACB",GSE26465$Symbol),][,c(2)],][,1])
OV90_ACACB

#FASN
GSE26465[grep("FASN",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("FASN",GSE26465$Symbol),][,c(2)],]
OV90_FASN<-mean(OV90_voom_result[GSE26465[grep("FASN",GSE26465$Symbol),][,c(2)],][,1])
OV90_FASN

#SCD
GSE26465[grep("SCD",GSE26465$Symbol),][,c(2,6)]
#SCD only
GSE26465[grep("^SCD$",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("^SCD$",GSE26465$Symbol),][,c(2)],]
OV90_SCD<-mean(OV90_voom_result[GSE26465[grep("^SCD$",GSE26465$Symbol),][,c(2)],][,1])
OV90_SCD

#CPT1A
GSE26465[grep("CPT1A",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("CPT1A",GSE26465$Symbol),][,c(2)],]
OV90_CPT1A<-mean(OV90_voom_result[GSE26465[grep("CPT1A",GSE26465$Symbol),][,c(2)],][,1])
OV90_CPT1A

#CPT1B
GSE26465[grep("CPT1B",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("CPT1B",GSE26465$Symbol),][,c(2)],]
OV90_CPT1B<-mean(OV90_voom_result[GSE26465[grep("CPT1B",GSE26465$Symbol),][,c(2)],][,1])
OV90_CPT1B

#CPT1C
GSE26465[grep("CPT1C",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("CPT1C",GSE26465$Symbol),][,c(2)],]
OV90_CPT1C<-mean(OV90_voom_result[GSE26465[grep("CPT1C",GSE26465$Symbol),][,c(2)],][,1])
OV90_CPT1C

#CPT2
GSE26465[grep("CPT2",GSE26465$Symbol),][,c(2,6)]
OV90_voom_result[GSE26465[grep("CPT2",GSE26465$Symbol),][,c(2)],]
OV90_CPT2<-mean(OV90_voom_result[GSE26465[grep("CPT2",GSE26465$Symbol),][,c(2)],][,1])
OV90_CPT2

#union
OV90_lpmeta_result<-as.data.frame(c(OV90_ACLY,OV90_ACACA,OV90_ACACB,OV90_FASN,OV90_SCD,OV90_CPT1A,OV90_CPT1B,OV90_CPT1C,OV90_CPT2))
OV90_lpmeta_result
rownames(OV90_lpmeta_result)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")
colnames(OV90_lpmeta_result)<-"logFC"

head(GSE26465_note)

GSE26465_lipidmeta_gene_probe<-rbind(
  GSE26465[grep("ACLY",GSE26465$Symbol),][,c(2,6)],
  GSE26465[grep("ACACA",GSE26465$Symbol),][,c(2,6)],
  GSE26465[grep("ACACB",GSE26465$Symbol),][,c(2,6)],
  GSE26465[grep("FSAN",GSE26465$Symbol),][,c(2,6)],
  GSE26465[grep("SCD",GSE26465$Symbol),][,c(2,6)],
  GSE26465[grep("CPT1A",GSE26465$Symbol),][,c(2,6)],
  GSE26465[grep("CPT1B",GSE26465$Symbol),][,c(2,6)],
  GSE26465[grep("CPT1C",GSE26465$Symbol),][,c(2,6)],
  GSE26465[grep("CPT2",GSE26465$Symbol),][,c(2,6)]
)
GSE26465_lipidmeta_gene_probe
rownames(GSE26465_lipidmeta_gene_probe)<-GSE26465_lipidmeta_gene_probe$Probe.Id
GSE26465_lipidmeta_gene_probe<-merge(GSE26465_lipidmeta_gene_probe,OV90_voom_result,by="row.names")
GSE26465_lipidmeta_gene_probe<-GSE26465_lipidmeta_gene_probe[,c(2,3,4)]

SKOV3_lpmeta_raw<-GSE26465_lipidmeta_gene_probe

SKOV3_lpmeta_raw

#ADR union

OV90_lpmeta_result
SKOV3_lpmeta_result
OVCAR3_lpmeta_result
TOV21G_lpmeta_result
A2780_lpmeta_result
W1_lpmeta_result

acquired_drug_resistance<-cbind(
  OV90_lpmeta_result,
  SKOV3_lpmeta_result,
  OVCAR3_lpmeta_result,
  TOV21G_lpmeta_result,
  A2780_lpmeta_result,
  W1_lpmeta_result
)

acquired_drug_resistance
colnames(acquired_drug_resistance)<-c("OV90","SKOV3","OVCAR3","TOV21G","A2780","W1")
acquired_drug_resistance

#heatmap plot
for_plot_acquired_drug_resistance<-acquired_drug_resistance
for_plot_acquired_drug_resistance<-cbind(rep(0,nrow(for_plot_acquired_drug_resistance)),for_plot_acquired_drug_resistance)
for_plot_acquired_drug_resistance<-rbind(rep(0,ncol(for_plot_acquired_drug_resistance)),for_plot_acquired_drug_resistance)

for_plot_acquired_drug_resistance<-cbind(for_plot_acquired_drug_resistance,rep(0,nrow(for_plot_acquired_drug_resistance)))
for_plot_acquired_drug_resistance<-rbind(for_plot_acquired_drug_resistance,rep(0,ncol(for_plot_acquired_drug_resistance)))


colnames(for_plot_acquired_drug_resistance)[1]<-""
rownames(for_plot_acquired_drug_resistance)[1]<-""

colnames(for_plot_acquired_drug_resistance)[ncol(for_plot_acquired_drug_resistance)]<-""
rownames(for_plot_acquired_drug_resistance)[nrow(for_plot_acquired_drug_resistance)]<-" "

for_plot_acquired_drug_resistance[1,1]<--3.5
for_plot_acquired_drug_resistance[nrow(for_plot_acquired_drug_resistance),ncol(for_plot_acquired_drug_resistance)]<-3.5

for_plot_acquired_drug_resistance

pheatmap(for_plot_acquired_drug_resistance,
         main = "Acquired drug resistance: changes in gene expression, log(FC,base=2)",
         color = colorRampPalette(c("Blue", "white","red"))(50),
         border="white",
         cluster_cols = F,
         cluster_rows = F,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize_row = 16,
         fontsize_col = 16,
         display_numbers = T,
         number_color = "black"
)

#PDR 
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism/GSE34615")
GSE34615_exp<-read.table("GSE34615_exp.txt",header = T)
head(GSE34615_exp)

#pac from https://www.cancerrxgene.org/compound/Paclitaxel/1080/overview/auc?tissue=OV
#doce from https://www.cancerrxgene.org/compound/Docetaxel/1007/overview/auc?tissue=OV

boxplot(GSE34615_exp[,-1])

GSE34615_note<-read.csv("GSE34615_note.csv",sep=",")
head(GSE34615_note)
colnames(GSE34615_note)[1]<-c("cell")
GSE34615_note$cell<-gsub("cell line: ","",GSE34615_note$cell)

AUC_pac<-read.csv("pac.csv",sep=",")
head(AUC_pac)
colnames(AUC_pac)[6]<-"AUC_pac"
AUC_doce<-read.csv("doce.csv",sep=",")
head(AUC_doce)
colnames(AUC_doce)[6]<-"AUC_doce"

length(AUC_doce$Cell.line)
length(AUC_pac$Cell.line)
length(intersect(AUC_doce$Cell.line,AUC_pac$Cell.line))
#AUC_doce$Cell.line == AUC_pac$Cell.line

AUC<-merge(AUC_doce,AUC_pac,by="Cell.line")
AUC<-AUC[,c("Cell.line","AUC_doce","AUC_pac")]
head(AUC)

AUC$Cell.line
GSE34615_note$cell

#Manual screening of available cell lines
AUC_cell_union<-c("OVCA420","OVCAR-4","OV-90","PEO1","SK-OV-3","FU-OV-1","OVCAR-5","OVCAR-3","TOV-21G","IGROV-1","TOV-112D","DOV13","OVCAR-8")   
exp_cell_union<-c("OVCA420","OVCAR4","OV90","PEO1","SKOV3","FUOV1","OVCAR5","OVCAR3","TOV21G","IGROV1","TOV112D","DOV13","OVCAR8")

cell_union<-as.data.frame(cbind(AUC_cell_union,exp_cell_union))
cell_union

#Primary drug resistance note
GPL15048<-getGEO("GPL15048",destdir = ".")

colnames(Table(GPL15048))
head(Table(GPL15048))
GPL15048<-Table(GPL15048)[,c(1,4)]
head(GPL15048)

#ACLY
GPL15048[grep("ACLY",GPL15048$GeneSymbol),][,c(1,2)]

#ACC1->ACACA, ACC2->ACACB
GPL15048[grep("ACACA",GPL15048$GeneSymbol),][,c(1,2)]
GPL15048[grep("ACACB",GPL15048$GeneSymbol),][,c(1,2)]

#FASN
GPL15048[grep("FASN",GPL15048$GeneSymbol),][,c(1,2)]

GPL15048[grep("SCD",GPL15048$GeneSymbol),][,c(1,2)]
#SCD only
GPL15048[grep("^SCD$",GPL15048$GeneSymbol),][,c(1,2)]

#CPT
GPL15048[grep("CPT1A",GPL15048$GeneSymbol),][,c(1,2)]
GPL15048[grep("CPT1B",GPL15048$GeneSymbol),][,c(1,2)]
#CPT1B only
GPL15048[grep("^CPT1B",GPL15048$GeneSymbol),][,c(1,2)]
GPL15048[grep("CPT1C",GPL15048$GeneSymbol),][,c(1,2)]

GPL15048[grep("CPT2",GPL15048$GeneSymbol),][,c(1,2)]

GPL15048_lipidmeta_gene_probe<-rbind(
  GPL15048[grep("ACLY",GPL15048$GeneSymbol),][,c(1,2)],
  GPL15048[grep("ACACA",GPL15048$GeneSymbol),][,c(1,2)],
  GPL15048[grep("ACACB",GPL15048$GeneSymbol),][,c(1,2)],
  GPL15048[grep("FASN",GPL15048$GeneSymbol),][,c(1,2)],
  GPL15048[grep("^SCD$",GPL15048$GeneSymbol),][,c(1,2)],
  GPL15048[grep("CPT1A",GPL15048$GeneSymbol),][,c(1,2)],
  GPL15048[grep("^CPT1B",GPL15048$GeneSymbol),][,c(1,2)],
  GPL15048[grep("CPT1C",GPL15048$GeneSymbol),][,c(1,2)],
  GPL15048[grep("CPT2",GPL15048$GeneSymbol),][,c(1,2)]
)
GPL15048_lipidmeta_gene_probe
colnames(GPL15048_lipidmeta_gene_probe)<-c("probe","gene")

head(GSE34615_exp)
cell_union<-as.data.frame(cbind(AUC_cell_union,exp_cell_union))
colnames(cell_union)[2]<-"cell"

GSE34615_note
cell_union

cell_union<-merge(cell_union,GSE34615_note,by="cell")
cell_union

colnames(GSE34615_exp)
cell_union

GSE34615_exp_sel<-GSE34615_exp[,c("ID_REF",cell_union$ID_REF)]
head(GSE34615_exp_sel)
rownames(GSE34615_exp_sel)<-GSE34615_exp_sel$ID_REF
GSE34615_exp_sel$ID_REF<-NULL

#ACLY
GPL15048[grep("ACLY",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("ACLY",GPL15048$GeneSymbol),][,c(1)],]
PDR_ACLY<-apply(GSE34615_exp_sel[GPL15048[grep("ACLY",GPL15048$GeneSymbol),][,c(1)],],
                2,mean)
PDR_ACLY

#ACACA
GPL15048[grep("ACACA",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("ACACA",GPL15048$GeneSymbol),][,c(1)],]
PDR_ACACA<-apply(GSE34615_exp_sel[GPL15048[grep("ACACA",GPL15048$GeneSymbol),][,c(1)],],
                 2,mean)
PDR_ACACA

#ACACB
GPL15048[grep("ACACB",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("ACACB",GPL15048$GeneSymbol),][,c(1)],]
PDR_ACACB<-apply(GSE34615_exp_sel[GPL15048[grep("ACACB",GPL15048$GeneSymbol),][,c(1)],],
                 2,mean)
PDR_ACACB

#FASN
GPL15048[grep("FASN",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("FASN",GPL15048$GeneSymbol),][,c(1)],]
PDR_FASN<-apply(GSE34615_exp_sel[GPL15048[grep("FASN",GPL15048$GeneSymbol),][,c(1)],],
                2,mean)
PDR_FASN

#SCD
GPL15048[grep("^SCD$",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("^SCD$",GPL15048$GeneSymbol),][,c(1)],]
PDR_SCD<-apply(GSE34615_exp_sel[GPL15048[grep("^SCD$",GPL15048$GeneSymbol),][,c(1)],],
               2,mean)
PDR_SCD

#CPT1A

#CPT1A
GPL15048[grep("CPT1A",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("CPT1A",GPL15048$GeneSymbol),][,c(1)],]
PDR_CPT1A<-apply(GSE34615_exp_sel[GPL15048[grep("CPT1A",GPL15048$GeneSymbol),][,c(1)],],
                 2,mean)
PDR_CPT1A

#CPT1B
GPL15048[grep("^CPT1B",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("^CPT1B",GPL15048$GeneSymbol),][,c(1)],]
PDR_CPT1B<-apply(GSE34615_exp_sel[GPL15048[grep("^CPT1B",GPL15048$GeneSymbol),][,c(1)],],
                 2,mean)
PDR_CPT1B

#CPT1C
GPL15048[grep("CPT1C",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("CPT1C",GPL15048$GeneSymbol),][,c(1)],]
PDR_CPT1C<-apply(GSE34615_exp_sel[GPL15048[grep("CPT1C",GPL15048$GeneSymbol),][,c(1)],],
                 2,mean)
PDR_CPT1C

#CPT2
GPL15048[grep("CPT2",GPL15048$GeneSymbol),][,c(1,2)]
GSE34615_exp_sel[GPL15048[grep("CPT2",GPL15048$GeneSymbol),][,c(1)],]
PDR_CPT2<-apply(GSE34615_exp_sel[GPL15048[grep("CPT2",GPL15048$GeneSymbol),][,c(1)],],
                2,mean)
PDR_CPT2

#union
PDR_lpmeta_result<-cbind(PDR_ACLY,PDR_ACACA,PDR_ACACB,PDR_FASN,PDR_SCD,PDR_CPT1A,PDR_CPT1B,PDR_CPT1C,PDR_CPT2)
PDR_lpmeta_result<-as.data.frame(PDR_lpmeta_result)
PDR_lpmeta_result

colnames(PDR_lpmeta_result)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")
#
rownames(cell_union)<-cell_union$ID_REF
cell_union
PDR_lpmeta_result<-merge(cell_union,PDR_lpmeta_result,by="row.names")

PDR_lpmeta_result<-PDR_lpmeta_result[,c("cell","ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")]
PDR_lpmeta_result
rownames(PDR_lpmeta_result)<-PDR_lpmeta_result$cell
PDR_lpmeta_result$cell<-NULL

PDR_lpmeta_result

#plot1
head(AUC)
cell_union
colnames(cell_union)[2]<-"Cell.line"

#Paclitaxel and docetaxel resistance are highly correlated
cor(AUC$AUC_doce,AUC$AUC_pac)

AUC<-AUC[order(AUC$AUC_pac),]

plot(AUC$AUC_pac,AUC$AUC_doce,xlab="AUC : Paclitaxel",ylab="AUC : Docetaxel",main="Drug resistance of ovarian cancer cell lines",cex.main=2,cex.lab=1.5,pch=10,cex=2,lwd=2,col="forestgreen")
text(AUC$AUC_pac,AUC$AUC_doce,AUC$Cell.line,cex = 0.6,pos = 1)
text(min(AUC$AUC_pac),max(AUC$AUC_doce),"Pearson correlation coefficient : 0.8408984",pos=4)

#plot2
for_plot_PDR_lpmeta_result<-PDR_lpmeta_result
cell_union
AUC
cell_union<-merge(AUC,cell_union,by="Cell.line")
for_plot_PDR_lpmeta_result$cell<-rownames(for_plot_PDR_lpmeta_result)
for_plot_PDR_lpmeta_result
for_plot_PDR_lpmeta_result<-merge(for_plot_PDR_lpmeta_result,cell_union,by="cell")

#plot2 PAC
for_plot_PDR_lpmeta_result_PAC<-for_plot_PDR_lpmeta_result[order(for_plot_PDR_lpmeta_result$AUC_pac),]
for_plot_PDR_lpmeta_result_PAC<-for_plot_PDR_lpmeta_result_PAC[,c(1:10,13)]
for_plot_PDR_lpmeta_result_PAC
rownames(for_plot_PDR_lpmeta_result_PAC)<-for_plot_PDR_lpmeta_result_PAC$cell
for_plot_PDR_lpmeta_result_PAC$cell<-NULL
for_plot_PDR_lpmeta_result_PAC

#Obtain Spearman correlation coefficient
Spearman_pac<-as.data.frame(rep(NA,9))
for_plot_PDR_lpmeta_result_PAC

for(i in 1:9){
  clum<-cor.test(for_plot_PDR_lpmeta_result_PAC[,i], sort(cell_union$AUC_pac),method = "spearman")
  a<-round(clum$p.value,3)
  Spearman_pac[i,1]<-a
}

Spearman_pac
colnames(Spearman_pac)<-"pac_p.val"
rownames(Spearman_pac)<-colnames(for_plot_PDR_lpmeta_result_PAC)[1:9]

#ACLY ACACA CPT1A cor

#all gene plot

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$ACLY,col="orchid",pch=8,cex=2,type="b",lwd=6,lty=3,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="ACLY : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$ACACA,col="forestgreen",pch=2,cex=2,type="b",lwd=6,lty=1,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="ACACA : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$ACACB,col="forestgreen",pch=2,cex=2,type="b",lwd=6,lty=1,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="ACACB : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$FASN,col="dodgerblue",pch=1,cex=2,type="b",lwd=6,lty=2,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="FASN : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$SCD,col="black",pch=10,cex=2,type="b",lwd=6,lty=3,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="SCD : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$CPT1A,col="red",pch=11,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="CPT1A : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$CPT1B,col="red",pch=11,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="CPT1B : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$CPT1C,col="red",pch=11,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="CPT1C : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$CPT2,col="blue",pch=12,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="CPT2 : Gene expression",cex.main=2,cex.lab=1.5)

#sel
plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$ACLY,col="orchid",pch=8,cex=2,type="b",lwd=6,lty=3,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="ACLY : Gene expression",cex.main=2,cex.lab=1.5)
plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$ACACA,col="forestgreen",pch=2,cex=2,type="b",lwd=6,lty=1,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="ACACA : Gene expression",cex.main=2,cex.lab=1.5)
plot(for_plot_PDR_lpmeta_result_PAC$AUC_pac,for_plot_PDR_lpmeta_result_PAC$CPT1A,col="red",pch=11,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Paclitaxel",ylab="Gene expression",
     main="CPT1A : Gene expression",cex.main=2,cex.lab=1.5)

#doce
for_plot_PDR_lpmeta_result_doce<-for_plot_PDR_lpmeta_result[order(for_plot_PDR_lpmeta_result$AUC_doce),]
for_plot_PDR_lpmeta_result_doce<-for_plot_PDR_lpmeta_result_doce[,c(1:10,12)]
for_plot_PDR_lpmeta_result_doce
rownames(for_plot_PDR_lpmeta_result_doce)<-for_plot_PDR_lpmeta_result_doce$cell
for_plot_PDR_lpmeta_result_doce$cell<-NULL
for_plot_PDR_lpmeta_result_doce


#Obtain Spearman correlation coefficient
Spearman_doce<-as.data.frame(rep(NA,9))
for_plot_PDR_lpmeta_result_doce

for(i in 1:9){
  clum<-cor.test(for_plot_PDR_lpmeta_result_doce[,i], sort(cell_union$AUC_doce),method = "spearman")
  a<-round(clum$p.value,3)
  Spearman_doce[i,1]<-a
}

Spearman_doce
colnames(Spearman_doce)<-"doce_p.val"
rownames(Spearman_doce)<-colnames(for_plot_PDR_lpmeta_result_doce)[1:9]

#FASN,CPT1A cor

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$ACLY,col="orchid",pch=8,cex=2,type="b",lwd=6,lty=3,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="ACLY : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$ACACA,col="forestgreen",pch=2,cex=2,type="b",lwd=6,lty=1,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="ACACA : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$ACACB,col="forestgreen",pch=2,cex=2,type="b",lwd=6,lty=1,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="ACACB : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$FASN,col="dodgerblue",pch=1,cex=2,type="b",lwd=6,lty=2,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="FASN : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$SCD,col="black",pch=10,cex=2,type="b",lwd=6,lty=3,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="SCD : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$CPT1A,col="red",pch=11,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="CPT1A : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$CPT1B,col="red",pch=11,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="CPT1B : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$CPT1C,col="red",pch=11,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="CPT1C : Gene expression",cex.main=2,cex.lab=1.5)

plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$CPT2,col="blue",pch=12,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="CPT2 : Gene expression",cex.main=2,cex.lab=1.5)

#sel
plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$FASN,col="dodgerblue",pch=1,cex=2,type="b",lwd=6,lty=2,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="FASN : Gene expression",cex.main=2,cex.lab=1.5)
plot(for_plot_PDR_lpmeta_result_doce$AUC_doce,for_plot_PDR_lpmeta_result_doce$CPT1A,col="red",pch=11,cex=2,type="b",lwd=6,lty=6,xlab="AUC : Docetaxel",ylab="Gene expression",
     main="CPT1A : Gene expression",cex.main=2,cex.lab=1.5)

#union
Spearman_doce
Spearman_pac

Spearman_PDR<-cbind(Spearman_pac,Spearman_doce)
Spearman_PDR

#Tracing the evolution process of drug resistance
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism/GSE159791")
GSE159791<-read.csv("GSE159791_raw_filtered.csv",sep=",")
head(GSE159791)

boxplot(log((GSE159791[,-1]+0.5),base=2))

#Batch effect needs to be removed
all(is.na(GSE159791)==F)
max(GSE159791[,-1])
min(GSE159791[,-1])
#need log 2

GSE159791[GSE159791$GeneID=="ACLY",]
GSE159791[GSE159791$GeneID=="ACACA",]
GSE159791[GSE159791$GeneID=="ACACB",]
GSE159791[GSE159791$GeneID=="FASN",]
GSE159791[GSE159791$GeneID=="SCD",]
GSE159791[GSE159791$GeneID=="CPT1A",]
#CPT1B no found
GSE159791[GSE159791$GeneID=="CPT1C",]
GSE159791[GSE159791$GeneID=="CPT2",]

#remove Batch effect
GSE159791[,2:ncol(GSE159791)]<-normalizeBetweenArrays(log((GSE159791[,-1]+0.5),base=2))
boxplot(GSE159791[,-1])

head(GSE159791)

evo<-rbind(
  GSE159791[GSE159791$GeneID=="ACLY",],
  GSE159791[GSE159791$GeneID=="ACACA",],
  GSE159791[GSE159791$GeneID=="ACACB",],
  GSE159791[GSE159791$GeneID=="FASN",],
  GSE159791[GSE159791$GeneID=="SCD",],
  GSE159791[GSE159791$GeneID=="CPT1A",],
  GSE159791[GSE159791$GeneID=="CPT1C",],
  GSE159791[GSE159791$GeneID=="CPT2",]
)

rownames(evo)<-evo$GeneID
evo$GeneID<-NULL
evo

evo_mean<-evo
evo_mean$C<-(evo$A2780a+evo$A2780b)/2
evo_mean$T4<-(evo_mean$A.4PTXa+evo_mean$A.4PTXb)/2
evo_mean$T8<-(evo_mean$A.8PTXa+evo_mean$A.8PTXb)/2
evo_mean$T16<-c(evo_mean$A.16PTXa+evo_mean$A.16PTXb)/2
evo_mean$T32<-c(evo_mean$A.32PTXa+evo_mean$A.32PTXb)/2
evo_mean$T64<-c(evo_mean$A.64PTXa+evo_mean$A.64PTXb)/2
evo_mean$T128<-c(evo_mean$A.128PTXa+evo_mean$A.128PTXb)/2

evo_mean<-evo_mean[,c("C","T4","T8","T16","T32","T64","T128")]
evo_mean

evo_mean_1<-evo_mean
colnames(evo_mean_1)<-c("A2780","A.4PTX","A.8PTX","A.16PTX","A.32PTX","A.64PTX","A.128PTX")
evo_mean_1

evo_mean_loged<-as.data.frame(t(evo_mean))
evo_mean_loged

evo_mean_loged$PTX<-c(0,4,8,16,32,64,128)

#Lipogenesis
plot(c(-1,128),c(min(evo_mean_loged[,1:4]),max(evo_mean_loged[,1:4])),col="white",pch=12,cex=2,type="b",lwd=6,lty=6,xlab="Tolerance to Paclitaxel (nM)",ylab="Gene expression",
     main="Lipid Biosynthesis related genes",cex.main=2,cex.lab=1.5)
lines(evo_mean_loged$PTX,evo_mean_loged$ACLY,col="blue",pch=12,cex=2,type="b",lwd=4,lty=4)
lines(evo_mean_loged$PTX,evo_mean_loged$ACACA,col="red",pch=1,cex=2,type="b",lwd=4,lty=1)
lines(evo_mean_loged$PTX,evo_mean_loged$ACACB,col="red",pch=1,cex=2,type="b",lwd=4,lty=1)
lines(evo_mean_loged$PTX,evo_mean_loged$FASN,col="forestgreen",pch=6,cex=2,type="b",lwd=4,lty=3)

text(64,evo_mean_loged[6,1],"ACLY",pos=3,cex=2)
text(64,evo_mean_loged[6,2],"ACACA",pos=3,cex=2)
text(64,evo_mean_loged[6,3],"ACACB",pos=3,cex=2)
text(64,evo_mean_loged[6,4],"FASN",pos=3,cex=2)

#Lipid Desaturation
plot(c(-1,128),c(min(evo_mean_loged[,5]),max(evo_mean_loged[,5])),col="white",pch=12,cex=2,type="b",lwd=6,lty=6,xlab="Tolerance to Paclitaxel (nM)",ylab="Gene expression",
     main="Lipid Desaturation related genes",cex.main=2,cex.lab=1.5)
lines(evo_mean_loged$PTX,evo_mean_loged$SCD,col="black",pch=4,cex=2,type="b",lwd=4,lty=4)
text(64,evo_mean_loged[6,5],"SCD",pos=3,cex=2)

#Lipid Catabolism related genes
plot(c(-1,128),c(min(evo_mean_loged[,6:8]),max(evo_mean_loged[,6:8])),col="white",pch=12,cex=2,type="b",lwd=6,lty=6,xlab="Tolerance to Paclitaxel (nM)",ylab="Gene expression",
     main="Lipid Catabolism related genes",cex.main=2,cex.lab=1.5)
lines(evo_mean_loged$PTX,evo_mean_loged$CPT1A,col="orchid",pch=8,cex=2,type="b",lwd=4,lty=4)
lines(evo_mean_loged$PTX,evo_mean_loged$CPT1C,col="orchid",pch=8,cex=2,type="b",lwd=4,lty=4)
lines(evo_mean_loged$PTX,evo_mean_loged$CPT2,col="dodgerblue",pch=11,cex=2,type="b",lwd=4,lty=1)
text(64,evo_mean_loged[6,6],"CPT1A",pos=3,cex=2)
text(64,evo_mean_loged[6,7],"CPT1C",pos=1,cex=2)
text(64,evo_mean_loged[6,8],"CPT2",pos=3,cex=2)

#clinic data
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism/GSE172016")
clinic<-read.table("GSE172016_organoid.gene.count.txt",header = T)
head(clinic)

boxplot(log((clinic[,2:8]+0.5),base=2))

#Batch effect needs to be removed
rownames(clinic)<-clinic$gene_id
head(GSE172016_note)

#T : UK1225,1226,1236
#C : UK1254,1393,2238,2326

GSE172016_note[grep("ACLY",GSE172016_note$hgnc_symbol),]
GSE172016_note[grep("ACACA",GSE172016_note$hgnc_symbol),]
GSE172016_note[grep("ACACB",GSE172016_note$hgnc_symbol),]
GSE172016_note[grep("FASN",GSE172016_note$hgnc_symbol),]
GSE172016_note[grep("SCD",GSE172016_note$hgnc_symbol),]
#SCD only
GSE172016_note[grep("^SCD$",GSE172016_note$hgnc_symbol),]

GSE172016_note[grep("CPT1A",GSE172016_note$hgnc_symbol),]
GSE172016_note[grep("CPT1B",GSE172016_note$hgnc_symbol),]
#CPT1B only
GSE172016_note[grep("^CPT1B$",GSE172016_note$hgnc_symbol),]

GSE172016_note[grep("CPT1C",GSE172016_note$hgnc_symbol),]
GSE172016_note[grep("CPT2",GSE172016_note$hgnc_symbol),]

#Batch effect needs to be removed
head(clinic)
clinic[,2:ncol(clinic)]<-normalizeBetweenArrays(log((clinic[,-1]+0.5),base=2))
boxplot(clinic[,-1])

head(clinic)

#Batch effect removed

clinic_sel<-rbind(
  clinic[grep(GSE172016_note[grep("ACLY",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),],
  clinic[grep(GSE172016_note[grep("ACACA",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),],
  clinic[grep(GSE172016_note[grep("ACACB",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),],
  clinic[grep(GSE172016_note[grep("FASN",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),],
  clinic[grep(GSE172016_note[grep("^SCD$",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),],
  clinic[grep(GSE172016_note[grep("CPT1A",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),],
  clinic[grep(GSE172016_note[grep("^CPT1B$",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),],
  clinic[grep(GSE172016_note[grep("CPT1C",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),],
  clinic[grep(GSE172016_note[grep("CPT2",GSE172016_note$hgnc_symbol),][,1],clinic$gene_id),])

clinic_sel
rownames(clinic_sel)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

clinic_sel$gene_id<-NULL

clinic_sel_loged<-as.data.frame(t(clinic_sel))

clinic_sel_loged

clinic_sel_loged$order<-c(6,5,7,3,4,2,1)
clinic_sel_loged<-clinic_sel_loged[order(clinic_sel_loged$order),]


clinic_sel_loged_T<-clinic_sel[,1:3]
clinic_sel_loged_T
apply(clinic_sel_loged_T,1,mean)

clinic_sel_loged_C<-clinic_sel[,4:7]
clinic_sel_loged_C
apply(clinic_sel_loged_C,1,mean)


#result
apply(clinic_sel_loged_T,1,mean)-apply(clinic_sel_loged_C,1,mean)



#for neu
head(A2780_for_limma)
head(W1_for_limma)

A2780_for_voom<-2^A2780_for_limma
W1_for_voom<-2^W1_for_limma

max(A2780_for_voom)
max(W1_for_voom)

head(A2780_for_voom)
head(W1_for_voom)
head(OVCAR3_for_voom)
head(TOV21G_for_voom)
head(SKOV3_for_voom)
head(OV90_for_voom)

boxplot(c(A2780_for_voom,W1_for_voom,OVCAR3_for_voom,TOV21G_for_voom,SKOV3_for_voom,OV90_for_voom))

# all the code is written by Liu-ZhanAo
# Liu-ZhanAo creates and conceived all the processes
A2780_for_neu<-A2780_for_voom
for(i in 1:ncol(A2780_for_voom)){
  A2780_for_neu[,i]<-(A2780_for_neu[,i]-mean(A2780_for_neu[,i]))/sd(A2780_for_neu[,i])
}
head(A2780_for_neu)

W1_for_neu<-W1_for_voom
for(i in 1:ncol(W1_for_voom)){
  W1_for_neu[,i]<-(W1_for_neu[,i]-mean(W1_for_neu[,i]))/sd(W1_for_neu[,i])
}
head(W1_for_neu)

OVCAR3_for_neu<-OVCAR3_for_voom
for(i in 1:ncol(OVCAR3_for_voom)){
  OVCAR3_for_neu[,i]<-(OVCAR3_for_neu[,i]-mean(OVCAR3_for_neu[,i]))/sd(OVCAR3_for_neu[,i])
}
head(OVCAR3_for_neu)

TOV21G_for_neu<-TOV21G_for_voom
for(i in 1:ncol(TOV21G_for_voom)){
  TOV21G_for_neu[,i]<-(TOV21G_for_neu[,i]-mean(TOV21G_for_neu[,i]))/sd(TOV21G_for_neu[,i])
}
head(TOV21G_for_neu)

SKOV3_for_neu<-SKOV3_for_voom
for(i in 1:ncol(SKOV3_for_voom)){
  SKOV3_for_neu[,i]<-(SKOV3_for_neu[,i]-mean(SKOV3_for_neu[,i]))/sd(SKOV3_for_neu[,i])
}
head(SKOV3_for_neu)

OV90_for_neu<-OV90_for_voom
for(i in 1:ncol(OV90_for_voom)){
  OV90_for_neu[,i]<-(OV90_for_neu[,i]-mean(OV90_for_neu[,i]))/sd(OV90_for_neu[,i])
}
head(OV90_for_neu)

boxplot(c(A2780_for_neu,W1_for_neu,OVCAR3_for_neu,TOV21G_for_neu,SKOV3_for_neu,OV90_for_neu))

head(A2780_for_neu)

#ACLY
GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_ACLY<-apply(A2780_for_neu[GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_ACLY

#ACACA
GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_ACACA<-apply(A2780_for_neu[GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_ACACA

#ACACB
GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_ACACB<-apply(A2780_for_neu[GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_ACACB

#FASN
GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_FASN<-apply(A2780_for_neu[GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_FASN

#SCD
GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_SCD<-apply(A2780_for_neu[GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_SCD

#CPT1A
GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_CPT1A<-apply(A2780_for_neu[GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_CPT1A

#CPT1B
GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_CPT1B<-apply(A2780_for_neu[GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_CPT1B

#CPT1C
GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_CPT1C<-apply(A2780_for_neu[GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_CPT1C

#CPT2
GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
A2780_for_neu[GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1)],]
A2780_neu_CPT2<-apply(A2780_for_neu[GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
A2780_neu_CPT2

A2780_neu<-rbind(
  A2780_neu_ACLY,
  A2780_neu_ACACA,
  A2780_neu_ACACB,
  A2780_neu_FASN,
  A2780_neu_SCD,
  A2780_neu_CPT1A,
  A2780_neu_CPT1B,
  A2780_neu_CPT1C,
  A2780_neu_CPT2
)

A2780_neu

#W1
head(W1_for_neu)

#ACLY
GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_ACLY<-apply(W1_for_neu[GPL13667[grep("ACLY",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_ACLY

#ACACA
GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_ACACA<-apply(W1_for_neu[GPL13667[grep("ACACA",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_ACACA

#ACACB
GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_ACACB<-apply(W1_for_neu[GPL13667[grep("ACACB",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_ACACB

#FASN
GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_FASN<-apply(W1_for_neu[GPL13667[grep("FASN",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_FASN

#SCD
GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_SCD<-apply(W1_for_neu[GPL13667[grep("^SCD$",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_SCD

#CPT1A
GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_CPT1A<-apply(W1_for_neu[GPL13667[grep("CPT1A",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_CPT1A

#CPT1B
GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_CPT1B<-apply(W1_for_neu[GPL13667[grep("CPT1B",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_CPT1B

#CPT1C
GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_CPT1C<-apply(W1_for_neu[GPL13667[grep("CPT1C",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_CPT1C

#CPT2
GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1,2,15,16,17)]
W1_for_neu[GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1)],]
W1_neu_CPT2<-apply(W1_for_neu[GPL13667[grep("CPT2",GPL13667$Gene.Symbol),][,c(1)],],2,mean)
W1_neu_CPT2

W1_neu<-rbind(
  W1_neu_ACLY,
  W1_neu_ACACA,
  W1_neu_ACACB,
  W1_neu_FASN,
  W1_neu_SCD,
  W1_neu_CPT1A,
  W1_neu_CPT1B,
  W1_neu_CPT1C,
  W1_neu_CPT2
)

A2780_neu
W1_neu

#
head(GSE172016_note)
head(TOV21G_for_neu)

TOV21G_neu<-rbind(
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="ACLY",][1,1],],
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="ACACA",][1,1],],
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="ACACB",][1,1],],
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="FASN",][1,1],],
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="SCD",][1,1],],
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="CPT1A",][1,1],],
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="CPT1B",][1,1],],
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="CPT1C",][1,1],],
  TOV21G_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="CPT2",][1,1],]
)

TOV21G_neu
rownames(TOV21G_neu)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

#
head(GSE172016_note)
head(OVCAR3_for_neu)

OVCAR3_neu<-rbind(
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="ACLY",][1,1],],
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="ACACA",][1,1],],
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="ACACB",][1,1],],
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="FASN",][1,1],],
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="SCD",][1,1],],
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="CPT1A",][1,1],],
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="CPT1B",][1,1],],
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="CPT1C",][1,1],],
  OVCAR3_for_neu[GSE172016_note[GSE172016_note$hgnc_symbol=="CPT2",][1,1],]
)

OVCAR3_neu
rownames(OVCAR3_neu)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

TOV21G_neu
OVCAR3_neu

#
#ACLY
GPL13287[grep("ACLY",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("ACLY",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_ACLY<-apply(SKOV3_for_neu[GPL13287[grep("ACLY",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_ACLY

#ACACA
GPL13287[grep("ACACA",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("ACACA",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_ACACA<-apply(SKOV3_for_neu[GPL13287[grep("ACACA",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_ACACA

#ACACB
GPL13287[grep("ACACB",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("ACACB",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_ACACB<-apply(SKOV3_for_neu[GPL13287[grep("ACACB",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_ACACB

#FASN
GPL13287[grep("FASN",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("FASN",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_FASN<-apply(SKOV3_for_neu[GPL13287[grep("FASN",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_FASN

#SCD
GPL13287[grep("^SCD$",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("^SCD$",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_SCD<-apply(SKOV3_for_neu[GPL13287[grep("^SCD$",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_SCD

#CPT1A
GPL13287[grep("CPT1A",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("CPT1A",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_CPT1A<-apply(SKOV3_for_neu[GPL13287[grep("CPT1A",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_CPT1A

#CPT1B
GPL13287[grep("CPT1B",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("CPT1B",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_CPT1B<-apply(SKOV3_for_neu[GPL13287[grep("CPT1B",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_CPT1B

#CPT1C
GPL13287[grep("CPT1C",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("CPT1C",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_CPT1C<-apply(SKOV3_for_neu[GPL13287[grep("CPT1C",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_CPT1C

#CPT2
GPL13287[grep("CPT2",GPL13287$Gene_symbol),][,c(1,8,11)]
SKOV3_for_neu[GPL13287[grep("CPT2",GPL13287$Gene_symbol),][,c(1)],]
SKOV3_neu_CPT2<-apply(SKOV3_for_neu[GPL13287[grep("CPT2",GPL13287$Gene_symbol),][,c(1)],],2,mean)
SKOV3_neu_CPT2

SKOV3_neu<-rbind(
  SKOV3_neu_ACLY,
  SKOV3_neu_ACACA,
  SKOV3_neu_ACACB,
  SKOV3_neu_FASN,
  SKOV3_neu_SCD,
  SKOV3_neu_CPT1A,
  SKOV3_neu_CPT1B,
  SKOV3_neu_CPT1C,
  SKOV3_neu_CPT2
)

SKOV3_neu

#OV90

#ACLY
GSE26465[grep("ACLY",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("ACLY",GSE26465$Symbol),][,c(2)],]
OV90_neu_ACLY<-apply(OV90_for_neu[GSE26465[grep("ACLY",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_ACLY

#ACACA
GSE26465[grep("ACACA",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("ACACA",GSE26465$Symbol),][,c(2)],]
OV90_neu_ACACA<-apply(OV90_for_neu[GSE26465[grep("ACACA",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_ACACA

#ACACB
GSE26465[grep("ACACB",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("ACACB",GSE26465$Symbol),][,c(2)],]
OV90_neu_ACACB<-apply(OV90_for_neu[GSE26465[grep("ACACB",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_ACACB

#FASN
GSE26465[grep("FASN",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("FASN",GSE26465$Symbol),][,c(2)],]
OV90_neu_FASN<-apply(OV90_for_neu[GSE26465[grep("FASN",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_FASN

#SCD
GSE26465[grep("^SCD$",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("^SCD$",GSE26465$Symbol),][,c(2)],]
OV90_neu_SCD<-apply(OV90_for_neu[GSE26465[grep("^SCD$",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_SCD

#CPT1A
GSE26465[grep("CPT1A",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("CPT1A",GSE26465$Symbol),][,c(2)],]
OV90_neu_CPT1A<-apply(OV90_for_neu[GSE26465[grep("CPT1A",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_CPT1A

#CPT1B
GSE26465[grep("CPT1B",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("CPT1B",GSE26465$Symbol),][,c(2)],]
OV90_neu_CPT1B<-apply(OV90_for_neu[GSE26465[grep("CPT1B",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_CPT1B

#CPT1C
GSE26465[grep("CPT1C",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("CPT1C",GSE26465$Symbol),][,c(2)],]
OV90_neu_CPT1C<-apply(OV90_for_neu[GSE26465[grep("CPT1C",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_CPT1C

#CPT2
GSE26465[grep("CPT2",GSE26465$Symbol),][,c(2,6)]
OV90_for_neu[GSE26465[grep("CPT2",GSE26465$Symbol),][,c(2)],]
OV90_neu_CPT2<-apply(OV90_for_neu[GSE26465[grep("CPT2",GSE26465$Symbol),][,c(2)],],2,mean)
OV90_neu_CPT2

OV90_neu<-rbind(
  OV90_neu_ACLY,
  OV90_neu_ACACA,
  OV90_neu_ACACB,
  OV90_neu_FASN,
  OV90_neu_SCD,
  OV90_neu_CPT1A,
  OV90_neu_CPT1B,
  OV90_neu_CPT1C,
  OV90_neu_CPT2
)

OV90_neu

#CPT1B does not exist in evolutionary data, so it is not included

A2780_neu
colnames(A2780_neu)<-c(rep("A2780_t",6),rep("A2780_c",3))
rownames(A2780_neu)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

W1_neu
colnames(W1_neu)<-c(rep("W1_t",3),rep("W1_c",3))
rownames(W1_neu)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

OVCAR3_neu
colnames(OVCAR3_neu)<-c(rep("OVCAR3_t",3),rep("OVCAR3_c",3))

TOV21G_neu
colnames(TOV21G_neu)<-c(rep("TOV21G_t",3),rep("TOV21G_c",3))

SKOV3_neu
colnames(SKOV3_neu)<-c(rep("SKOV3_t",4),rep("SKOV3_c",2))
rownames(SKOV3_neu)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

OV90_neu
colnames(OV90_neu)<-c(rep("OV90_t",8),rep("OV90_c",2))
rownames(OV90_neu)<-c("ACLY","ACACA","ACACB","FASN","SCD","CPT1A","CPT1B","CPT1C","CPT2")

neu_union<-cbind(
  A2780_neu,
  W1_neu,
  OVCAR3_neu,
  TOV21G_neu,
  SKOV3_neu,
  OV90_neu
)
neu_union

neu_union<-as.data.frame(t(neu_union))
neu_union$label<-c(1,1,1,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,0)

neu_union
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism/for_neural_network")
write.csv(neu_union,"ADR_neu_union.csv")

#An Attempt at Multivariate Regression Model
par(mfrow=c(1,1))
par(mfrow=c(2,2))

set.seed(1234)

shrinkage10<-function(fit,k=10){
  require(bootstrap)
  theta.fit<-function(x,y){lsfit(x,y)}
  theta.predict<-function(fit,x){cbind(1,x)%*%fit$coef}
  x<-fit$model[,2:ncol(fit$model)]
  y<-fit$model[,1]
  results<-crossval(x,y,theta.fit,theta.predict,ngroup = k)
  r2<-cor(y,fit$fitted.values)^2
  r2cv<-cor(y,results$cv.fit)^2
  cat("Original R-square =",r2,"\n")
  cat(k,"Fold Cross-Valiated R-square =",r2cv,"\n")
  cat("Change =",r2-r2cv,"\n")
}

shrinkage5<-function(fit,k=5){
  require(bootstrap)
  theta.fit<-function(x,y){lsfit(x,y)}
  theta.predict<-function(fit,x){cbind(1,x)%*%fit$coef}
  x<-fit$model[,2:ncol(fit$model)]
  y<-fit$model[,1]
  results<-crossval(x,y,theta.fit,theta.predict,ngroup = k)
  r2<-cor(y,fit$fitted.values)^2
  r2cv<-cor(y,results$cv.fit)^2
  cat("Original R-square =",r2,"\n")
  cat(k,"Fold Cross-Valiated R-square =",r2cv,"\n")
  cat("Change =",r2-r2cv,"\n")
}


#fit1
fit1<-lm(label~.,data=neu_union)
summary(fit1)
AIC(fit1)
confint(fit1)
plot(fit1)

#fit2
fit2<-lm(label~ACLY+ACACA+SCD+CPT1A+CPT1B,data=neu_union)
summary(fit2)
AIC(fit2)
confint(fit2)
plot(fit2)

#fit3
fit3<-glm(label~ACLY+ACACA+SCD+CPT1A+CPT1B,data=neu_union,family = binomial())
summary(fit3)
AIC(fit3)

#fit4
fit4<-glm(label~ACLY+ACACA+SCD+CPT1B,data=neu_union,family = binomial())
summary(fit4)
AIC(fit4)
plot(fit4)
exp(coef(fit4))

#fit5
fit5<-glm(label~ACLY+ACACA+SCD,data=neu_union,family = binomial())
summary(fit5)
AIC(fit5)
#fit5 bad

#fit6
fit6<-glm(label~.,data=neu_union,family = binomial())
summary(fit6)
AIC(fit6)
#fit6 bad

shrinkage5(fit1)
shrinkage5(fit2)
shrinkage5(fit3)
shrinkage5(fit4)
shrinkage5(fit5)
shrinkage5(fit6)

#use fit3
summary(fit3)

library(vcd)

clum3<-rep(NA,100)

for(i in 1:100){
  set.seed(i)
  allframe<-neu_union
  allframe<-allframe[sample(length(allframe$ACLY)),]
  trainframe<-allframe[1:30,]
  verifyframe<-allframe[31:43,]
  fit3_train<-glm(label~ACLY+ACACA+SCD+CPT1A+CPT1B,data=trainframe,family = binomial())
  
  verifyframe_1<-verifyframe[,c(1,2,5,6,7)]
  verifyframe_1$result<-predict(fit3_train,verifyframe_1,type="response")
  verifyframe_1<-cbind(verifyframe_1,verifyframe$label)
  verifyframe_1$result_1<-NA
  verifyframe_1$result_1[verifyframe_1$result>0.5]<-1
  verifyframe_1$result_1[verifyframe_1$result<=0.5]<-0
  
  k<-Kappa(table(verifyframe_1$`verifyframe$label`,verifyframe_1$result_1),weights = "Equal-Spacing")
  clum3[i]<-k$Unweighted[1]
}

clum3
mean(clum3)

summary(fit3)

acc_all<-rep(NA,100)

for(i in 1:100){
  set.seed(i)
  allframe<-neu_union
  allframe<-allframe[sample(length(allframe$ACLY)),]
  trainframe<-allframe[1:30,]
  verifyframe<-allframe[31:43,]
  fit3_train<-glm(label~ACLY+ACACA+SCD+CPT1A+CPT1B,data=trainframe,family = binomial())
  
  verifyframe_1<-verifyframe[,c(1,2,5,6,7)]
  verifyframe_1$result<-predict(fit3_train,verifyframe_1,type="response")
  verifyframe_1<-cbind(verifyframe_1,verifyframe$label)
  verifyframe_1$result_1<-NA
  verifyframe_1$result_1[verifyframe_1$result>0.5]<-1
  verifyframe_1$result_1[verifyframe_1$result<=0.5]<-0
  
  t<-table(verifyframe_1$`verifyframe$label`,verifyframe_1$result_1)
  t1<-(t[1,1]+t[2,2])/(t[1,1]+t[2,2]+t[1,2]+t[2,1])
  acc_all[i]<-t1
}

acc_all

mean(acc_all)
setwd("C:/Users/86186/Desktop/RLM_Reprogramming_of_lipid_metabolism")

save.image()


