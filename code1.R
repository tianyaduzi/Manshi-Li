TCGA_PAAD_mRNA_FPKM <- read.delim("TCGA_PAAD_mRNA_FPKM.txt",check.names = F,row.names = 1)
TCGA_PAAD_mRNA_counts <- read.delim("TCGA_PAAD_mRNA_counts.txt",header = T,row.names = 1,check.names = F)
sample_ID <- colnames(TCGA_PAAD_mRNA_counts)
tumor_sample_ID=grep(".*-0\\d[ABC]-.*",sample_ID,value = T)
normal_sample_ID=grep(".*-1\\d[ABC]-.*",sample_ID,value = T)
TCGA_PAAD_mRNA_counts <- TCGA_PAAD_mRNA_counts[,c(tumor_sample_ID,normal_sample_ID)]
library("edgeR")
group_list<-c(rep("tumor",179),rep("normal",4))
mRNA_DGE<-DGEList(counts=TCGA_PAAD_mRNA_counts,group =group_list)
keep.exprs<- edgeR::filterByExpr(mRNA_DGE, group=group_list)
mRNA_DGE<- mRNA_DGE[keep.exprs,,keep.lib.sizes=FALSE]
mRNA_DGE <- calcNormFactors(mRNA_DGE)
mRNA_DGE_bcv <- mRNA_DGE
bcv <- 0.4
et_mRNA <- exactTest(mRNA_DGE_bcv, dispersion = bcv^2)
res <- et_mRNA$table
res$FDR <- p.adjust(res$PValue,method = "BH",length(res$PValue))#p.adjust进行
mydata <- data.frame(rownames(res),res[,c(1,3,4)])
colnames(mydata) <- c("mRNA","logFC","p","FDR")
rownames(mydata) <- 1:nrow(mydata)
mydata$condition <- ifelse(mydata$logFC>=1&mydata$FDR<0.05,"up",
                                ifelse(mydata$logFC<=-1&mydata$FDR<0.05,"down","normal"))
library("ggplot2")
ggplot(data=mydata, aes(x=logFC, y=-log10(FDR), colour=condition))+
  geom_point(alpha=0.8, size=1.5)+  
  labs(x="log2 Fold Change",y="-log10 FDR")+
  geom_hline(yintercept=-log10(0.05),linetype=4,linewidth=1.2,colour="grey")+
  geom_vline(xintercept=c(-1,1),linetype=4,linewidth=1.2,colour="grey")+
  scale_color_manual(values=c('up'='red','down'='green','normal'='gray'))+ 
  scale_x_continuous(breaks=seq(-7, 9, 2))+
  theme(axis.title = element_text(face = "bold",size = 20), 
        axis.text  = element_text(face = "bold",size = 13), 
        #plot.title = element_text(size = 20,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 15,face = "bold",), 
        legend.text = element_text(size = 15,face = "bold",),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,linewidth=2, colour = "black"),
        axis.title.x = element_text(margin = margin(0.4,1,0,1,'cm')),
        axis.title.y = element_text(margin = margin(0,0.4,0,0,'cm')),
        axis.ticks.length = unit(10,"pt"),
        axis.ticks = element_line(linewidth = 1.2),
        text =element_text(family = "serif"),
  )
DEmRNA <- mydata[-which(mydata$condition=="normal"),]
String_protein_interaction <- read.delim("String_protein_interaction.txt")
String_protein_interaction <- String_protein_interaction[,c(3,6,7)]
gene1 <- data.frame(DEmRNA$mRNA)
colnames(gene1) <- c("gene_name1")
gene2 <- data.frame(DEmRNA$mRNA)
colnames(gene2) <- c("gene_name2")
DEmRNA_network1 <- merge(gene1,String_protein_interaction,by="gene_name1")
DEmRNA_network2 <- merge(gene2,DEmRNA_network1,by="gene_name2")
DEmRNA_DEmRNA_unique <- unique(as.data.frame(t(apply(DEmRNA_network2,1,sort))))
DEmRNA_DEmRNA_unique <- DEmRNA_DEmRNA_unique[,c(2,3,1)]
DEmRNA_DEmRNA_unique <- DEmRNA_DEmRNA_unique[!duplicated(DEmRNA_DEmRNA_unique),]
colnames(DEmRNA_DEmRNA_unique) <- c("gene_name1","gene_name2","score")
res <- data.frame()
for(i in 1:nrow(DEmRNA_DEmRNA_unique)){
  a <- as.numeric(TCGA_PAAD_mRNA_FPKM[DEmRNA_DEmRNA_unique[i,1],])
  b <- as.numeric(TCGA_PAAD_mRNA_FPKM[DEmRNA_DEmRNA_unique[i,2],])
  p=cor.test(as.numeric(a),as.numeric(b),method = "pearson")$p.value
  x=c(DEmRNA_DEmRNA_unique[i,1],DEmRNA_DEmRNA_unique[i,2],DEmRNA_DEmRNA_unique[i,3],as.numeric(p))
  res=rbind(res,x)
  a=c()
  b=c()
  p=c()
  print(i)
}
colnames(res) <- c("gene_name1","gene_name2","score","p")
DEmRNA_DEmRNA_network <- res[which(res$p<0.05),]
gene_last <- unique(c(DEmRNA_DEmRNA_network$gene_name1,DEmRNA_DEmRNA_network$gene_name2))
weight_values_TCGA_mRNA <- read.delim("weight_values_TCGA_mRNA.txt")
rownames(weight_values_TCGA_mRNA) <- weight_values_TCGA_mRNA$X
sample <- sample(weight_values_TCGA_mRNA$X)
weight_values_TCGA_mRNA <- weight_values_TCGA_mRNA[sample,]
weight_values_TCGA_mRNA$lable <- c(1:nrow(weight_values_TCGA_mRNA))
library("ggplot2")
library("ggpubr")
library(ggExtra)
library(cowplot)
p <- ggplot(weight_values_TCGA_mRNA,aes(weight,lable))+geom_point(color="#F7903D",size=3)+
  geom_vline(xintercept=0.0017,linetype=4,lty=3,lwd=1,col="black")+
  labs(
    x="weight_value",
    y="mRNA"
  )
p
p1 <- p+theme(legend.position = "none")+scale_x_continuous(breaks=seq(0,0.05,0.0025))+scale_y_continuous(breaks=seq(0,1074,40))+
  theme(axis.title = element_text(face = "bold",size = 30), 
        text =element_text(family = "serif"),
        panel.border = element_rect(fill = NA,linewidth=2, colour = "black"),
        axis.text.y = element_text(face = "bold",size = 15),
        axis.text.x = element_text(face = "bold",size = 20,angle = 90,vjust = 0.5,hjust = 0.5))
ggMarginal(p1, 
           type = "histogram",
           xparams=list(color="black",fill= "#4D85BD"),
           yparams = list(color="black",fill= "#59A95A")
)
top_gene1 <- weight_values_TCGA_mRNA[which(weight_values_TCGA_mRNA$weight>=0.0017),1]
library(stringr)
library(ggplot2) 
library(ggrepel)
top_gene1 <- weight_values_TCGA_mRNA[which(weight_values_TCGA_mRNA$weight>=0.0017),]
ggplot(data = top_gene1, 
                   aes(x = reorder(X,weight), y = weight,fill=c("blue")))+ 
  geom_bar(stat = "identity",width = 0.9)+ 
  coord_flip()+ 
  ylim(0,0.05)+
  labs(x="",y = "weight-value")+ 
  theme(axis.title = element_text(face = "bold",size = 20),
        axis.text = element_text(face = "bold",size = 13),
        plot.title = element_text(size = 20,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 15,face = "bold",),
        legend.text = element_text(size = 15,face = "bold",),
        panel.background = element_blank(),
        axis.line.x.bottom = element_line(size=1, colour = "black"),
        axis.line.y.left = element_line(size=1, colour = "black"),
        axis.title.x = element_text(margin = margin(0.2,0,0,0,'cm')),
        axis.title.y = element_text(margin = margin(0,0.4,0,0,'cm')),
        axis.ticks.length = unit(10,"pt"),
        axis.ticks = element_line(size = 1.2),
        text =element_text(family = "serif"),
        legend.position = "none",
  )
rownames(network_analysis) <- network_analysis$node_name
sample <- sample( network_analysis$node_name)
network_analysis <- network_analysis[sample,]
network_analysis$lable <- c(1:nrow(network_analysis))
p <- ggplot(network_analysis,aes(MCC,lable))+geom_point(color="#F7903D",size=3)+
  geom_vline(xintercept=2,linetype=4,lty=3,lwd=1,col="black")+
  labs(
    x="MCC",
    y="mRNA"
  )
p
p1 <- p+theme(legend.position = "none")+scale_x_continuous(breaks=seq(0,54,4))+scale_y_continuous(breaks=seq(0,1074,40))+
  theme(axis.title = element_text(face = "bold",size = 30), 
        axis.text.y = element_text(face = "bold",size = 15),
        text =element_text(family = "serif"),
        axis.text.x = element_text(face = "bold",size = 20,angle = 90,vjust = 0.5,hjust = 0.5))
ggMarginal(p1, 
           type = "histogram",
           xparams=list(color="black",fill= "#4D85BD"),
           yparams = list(color="black",fill= "#59A95A")
)
top_gene2 <-network_analysis[which(network_analysis$MCC>=2),1]
cadidata_mRNA_sum <- intersect(top_gene1,top_gene2)
clinical = read.csv(file = "clinical.tsv", sep = "\t" , header = T)
clinical$days_to_death=as.numeric(clinical$days_to_death)
clinical$days_to_death[is.na(clinical$days_to_death)]=0
clinical$days_to_last_follow_up=as.numeric(clinical$days_to_last_follow_up)
clinical$days_to_last_follow_up[is.na(clinical$days_to_last_follow_up)]=0
clinical$days=as.numeric(clinical$days_to_death)+as.numeric(clinical$days_to_last_follow_up)
clinical$OS_time=round(clinical$days/30,2)
clinical$OS_statue=ifelse(clinical$vital_status=='Alive',0,1)
surv_data <- clinical[,c("case_submitter_id","OS_time","OS_statue")]
surv_data <- surv_data[!duplicated(surv_data),]
sample_id <- read.csv(file = "gdc_sample_sheet.2022-05-28.tsv", sep = "\t" , header = T)
sample_id <- sample_id[,c(6,7)]
colnames(sample_id) <- c("case_submitter_id","sample_id")
surv_data <- merge(surv_data,sample_id,by="case_submitter_id")
surv_data <- surv_data[,c(4,2,3)]
top_gene_exp <- TCGA_PAAD_mRNA_FPKM[c("CEL","PRSS1","CPB1","ELAPOR1"),tumor_sample_ID]
top_gene_exp <- data.frame(t(top_gene_exp))
top_gene_exp$sample_id <- substr(rownames(top_gene_exp),1,16)
surv_data <- merge(surv_data,top_gene_exp,by="sample_id")
library(survival)
library(survminer)
sur.cat = surv_cutpoint(surv_data, time = "OS_time",
                        event = "OS_statue",
                        variables = c("CEL","PRSS1","CPB1","ELAPOR1"))
print(sur.cat)
sur.cat = surv_categorize(sur.cat)
log_rank_res <- data.frame() 
for(i in 3:ncol(sur.cat)){
  fit <- survfit(Surv(OS_time,OS_statue) ~sur.cat[,i], data = sur.cat) 
  pValue = surv_pvalue(fit)$pval
  x <- c(colnames(sur.cat)[i],pValue)
  log_rank_res <- rbind(log_rank_res,x)
  print(i)
}
colnames(log_rank_res) <- c("mRNA","p")
sfit <- surv_fit(Surv(OS_time,OS_statue)~ELAPOR1,data = sur.cat)
ggsurvplot(sfit,
           xlab=("Times(months)"),
           legend = c(0.8,0.75),
           legend.labs=c("High","Low"), 
           legend.title="Classifaction",
           xlim=c(0,60),
           break.time.by = 12,
           pval = TRUE, 
           pval.size=6,
           risk.table = TRUE, 
           risk.table.col = "black", 
           surv.median.line = "hv",
           ggtheme = theme(axis.title = element_text(face = "bold",size = 20), 
                           axis.text  = element_text(face = "bold",size = 13), 
                           plot.title = element_text(size = 20,hjust = 0.5,face = "bold"),
                           legend.title = element_text(size = 15,face = "bold",),
                           legend.text = element_text(size = 15,face = "bold",),
                           panel.background = element_blank(),
                           axis.line.x = element_line(linewidth=1, colour = "black"),
                           axis.line.y = element_line(linewidth=1, colour = "black"),
                           text=element_text(family="serif")
           tables.theme = theme(axis.title = element_text(face = "bold",size = 15), 
                                axis.text  = element_text(face = "bold",size = 10), 
                                axis.line = element_line(linewidth=1, colour = "black"),
                                plot.title = element_text(size = 15,hjust = 0.5,face = "bold") 
           ),
           title="ELAPOR1 Survival Analysis")

clinical1 <- read.delim("TCGA-PAADclinc.txt")
radiation_sample <- clinical1[which(clinical1$treatment_or_therapy=="yes"&clinical1$treatment_type=="Radiation Therapy, NOS"),]
noradiation_sample <- clinical1[which(clinical1$treatment_or_therapy=="yes"&clinical1$treatment_type=="Pharmaceutical Therapy, NOS"),]
colnames(TCGA_PAAD_mRNA_counts) <- substr(colnames(TCGA_PAAD_mRNA_counts),1,16)
library("edgeR")
one_year_pre <- radiation_sample[which(radiation_sample$OS_year<=1),1]
one_year_past <- radiation_sample[which(radiation_sample$OS_year>1),1]
exp <- TCGA_PAAD_mRNA_counts[,c(one_year_pre,one_year_past)]
group_list<-c(rep("pre",length(one_year_pre)),rep("past",length(one_year_past)))
mRNA_DGE<-DGEList(counts=exp,group =group_list)
keep.exprs<- edgeR::filterByExpr(mRNA_DGE, group=group_list)
mRNA_DGE<- mRNA_DGE[keep.exprs,,keep.lib.sizes=FALSE]
mRNA_DGE <- calcNormFactors(mRNA_DGE)
mRNA_DGE_bcv <- mRNA_DGE
bcv <- 0.4
et_mRNA <- exactTest(mRNA_DGE_bcv, dispersion = bcv^2)
one_year_res <- et_mRNA$table
one_year_res$FDR <- p.adjust(one_year_res$PValue,method = "BH",length(one_year_res$PValue))
one_year_DEmRNA <- data.frame(rownames(one_year_res),one_year_res[,c(1,3,4)])
colnames(one_year_DEmRNA) <- c("mRNA","logFC","p","FDR")
rownames(one_year_DEmRNA) <- 1:nrow(one_year_DEmRNA)
one_year_DEmRNA$condition <- ifelse(one_year_DEmRNA$logFC>=1&one_year_DEmRNA$FDR<0.05,"up",
                           ifelse(one_year_DEmRNA$logFC<=-1&one_year_DEmRNA$FDR<0.05,"down","normal"))
one_year_DEmRNA <- one_year_DEmRNA[-which(one_year_DEmRNA$condition=="normal"),]
two_year_pre <- radiation_sample[which(radiation_sample$OS_year<=2),1]
two_year_past <- radiation_sample[which(radiation_sample$OS_year>2),1]
exp <- TCGA_PAAD_mRNA_counts[,c(two_year_pre,two_year_past)]
group_list<-c(rep("pre",length(two_year_pre)),rep("past",length(two_year_past)))
mRNA_DGE<-DGEList(counts=exp,group =group_list)
keep.exprs<- edgeR::filterByExpr(mRNA_DGE, group=group_list)
mRNA_DGE<- mRNA_DGE[keep.exprs,,keep.lib.sizes=FALSE]
mRNA_DGE <- calcNormFactors(mRNA_DGE)
mRNA_DGE_bcv <- mRNA_DGE
bcv <- 0.4
et_mRNA <- exactTest(mRNA_DGE_bcv, dispersion = bcv^2)
two_year_res <- et_mRNA$table
two_year_res$FDR <- p.adjust(two_year_res$PValue,method = "BH",length(two_year_res$PValue))
two_year_DEmRNA <- data.frame(rownames(two_year_res),two_year_res[,c(1,3,4)])
colnames(two_year_DEmRNA) <- c("mRNA","logFC","p","FDR")
rownames(two_year_DEmRNA) <- 1:nrow(two_year_DEmRNA)
two_year_DEmRNA$condition <- ifelse(two_year_DEmRNA$logFC>=1&two_year_DEmRNA$FDR<0.05,"up",
                                    ifelse(two_year_DEmRNA$logFC<=-1&two_year_DEmRNA$FDR<0.05,"down","normal"))
two_year_DEmRNA <- two_year_DEmRNA[-which(two_year_DEmRNA$condition=="normal"),]
three_year_pre <- radiation_sample[which(radiation_sample$OS_year<=3),1]
three_year_past <- radiation_sample[which(radiation_sample$OS_year>3),1]
exp <- TCGA_PAAD_mRNA_counts[,c(three_year_pre,three_year_past)]
group_list<-c(rep("pre",length(three_year_pre)),rep("past",length(three_year_past)))
mRNA_DGE<-DGEList(counts=exp,group =group_list)
keep.exprs<- edgeR::filterByExpr(mRNA_DGE, group=group_list)
mRNA_DGE<- mRNA_DGE[keep.exprs,,keep.lib.sizes=FALSE]
mRNA_DGE <- calcNormFactors(mRNA_DGE)
mRNA_DGE_bcv <- mRNA_DGE
bcv <- 0.4
et_mRNA <- exactTest(mRNA_DGE_bcv, dispersion = bcv^2)
three_year_res <- et_mRNA$table
three_year_res$FDR <- p.adjust(three_year_res$PValue,method = "BH",length(three_year_res$PValue))
three_year_DEmRNA <- data.frame(rownames(three_year_res),three_year_res[,c(1,3,4)])
colnames(three_year_DEmRNA) <- c("mRNA","logFC","p","FDR")
rownames(three_year_DEmRNA) <- 1:nrow(three_year_DEmRNA)
three_year_DEmRNA$condition <- ifelse(three_year_DEmRNA$logFC>=1&three_year_DEmRNA$FDR<0.05,"up",
                                      ifelse(three_year_DEmRNA$logFC<=-1&three_year_DEmRNA$FDR<0.05,"down","normal"))
three_year_DEmRNA <- three_year_DEmRNA[-which(three_year_DEmRNA$condition=="normal"),]
mean_year_pre <- radiation_sample[which(radiation_sample$OS_year<=mean(radiation_sample$OS_year)),1]
mean_year_past <- radiation_sample[which(radiation_sample$OS_year>mean(radiation_sample$OS_year)),1]
exp <- TCGA_PAAD_mRNA_counts[,c(mean_year_pre,mean_year_past)]
group_list<-c(rep("pre",length(mean_year_pre)),rep("past",length(mean_year_past)))
mRNA_DGE<-DGEList(counts=exp,group =group_list)
keep.exprs<- edgeR::filterByExpr(mRNA_DGE, group=group_list)
mRNA_DGE<- mRNA_DGE[keep.exprs,,keep.lib.sizes=FALSE]
mRNA_DGE <- calcNormFactors(mRNA_DGE)
mRNA_DGE_bcv <- mRNA_DGE
bcv <- 0.4
et_mRNA <- exactTest(mRNA_DGE_bcv, dispersion = bcv^2)
mean_year_res <- et_mRNA$table
mean_year_res$FDR <- p.adjust(mean_year_res$PValue,method = "BH",length(mean_year_res$PValue))
mean_year_DEmRNA <- data.frame(rownames(mean_year_res),mean_year_res[,c(1,3,4)])
colnames(mean_year_DEmRNA) <- c("mRNA","logFC","p","FDR")
rownames(mean_year_DEmRNA) <- 1:nrow(mean_year_DEmRNA)
mean_year_DEmRNA$condition <- ifelse(mean_year_DEmRNA$logFC>=1&mean_year_DEmRNA$FDR<0.05,"up",
                                    ifelse(mean_year_DEmRNA$logFC<=-1&mean_year_DEmRNA$FDR<0.05,"down","normal"))

mean_year_DEmRNA1 <- mean_year_DEmRNA[-which(mean_year_DEmRNA$condition=="normal"),]
DEmRNA_last <- intersect(one_year_DEmRNA$mRNA,two_year_DEmRNA$mRNA)
DEmRNA_last <- intersect(DEmRNA_last,three_year_DEmRNA$mRNA)
gene1 <- data.frame(DEmRNA_last)
colnames(gene1) <- c("gene_name1")
gene2 <- data.frame(DEmRNA_last)
colnames(gene2) <- c("gene_name2")
DEmRNA_last_network1 <- merge(gene1,String_protein_interaction,by="gene_name1")
DEmRNA_last_network2 <- merge(gene2,DEmRNA_last_network1,by="gene_name2")
DEmRNA_last_network_unique <- unique(as.data.frame(t(apply(DEmRNA_last_network2,1,sort))))
DEmRNA_last_network_unique <- DEmRNA_last_network_unique[,c(2,3,1)]
DEmRNA_last_network_unique <- DEmRNA_last_network_unique[!duplicated(DEmRNA_last_network_unique),]
colnames(DEmRNA_last_network_unique) <- c("gene_name1","gene_name2","score")
DEmRNA_last_network_analysis <- read.csv("DEmRNA_last_network.csv")
rownames(DEmRNA_last_network_analysis) <- DEmRNA_last_network_analysis$node_name
sample <- sample(DEmRNA_last_network_analysis)
DEmRNA_last_network_analysis <- DEmRNA_last_network_analysis[sample,]
DEmRNA_last_network_analysis$lable <- c(1:nrow(DEmRNA_last_network_analysis))
library("ggplot2")
library("ggpubr")
library(ggExtra)
library(cowplot)
p <- ggplot(DEmRNA_last_network_analysis,aes(MCC,lable))+geom_point(color="#F7903D",size=3)+
  geom_vline(xintercept=9,linetype=4,lty=3,lwd=1,col="black")+
  labs(
    x="MCC",
    y="mRNA"
  )
p
p1 <- p+theme(legend.position = "none")+
  scale_x_continuous(breaks=seq(0,240,10))+
  scale_y_continuous(breaks=seq(0,120,10))+
  theme(axis.title = element_text(face = "bold",size = 30), 
        axis.text.y = element_text(face = "bold",size = 15),
        text =element_text(family = "serif"),
        axis.text.x = element_text(face = "bold",size = 20,angle = 90,vjust = 0.5,hjust = 0.5))
ggMarginal(p1, 
           type = "histogram",
           xparams=list(color="black",fill= "#4D85BD"),
           yparams = list(color="black",fill= "#59A95A")
)
cadidata_mRNA_radiation <- DEmRNA_last_network_analysis[which(DEmRNA_last_network_analysis$MCC>=9),1]
cadidata_mRNA <- intersect(cadidata_mRNA_radiation,cadidata_mRNA_sum)
surv_data1 <- merge(surv_data,radiation_sample,by="sample_id")
sur.cat1 = surv_cutpoint(surv_data1, time = "OS_time.x",
                        event = "OS_statue",
                        variables = c("CEL","PRSS1","CPB1","ELAPOR1"))
print(sur.cat1)
sur.cat1 = surv_categorize(sur.cat1)
log_rank_res1 <- data.frame() 
for(i in 3:ncol(sur.cat1)){
  fit <- survfit(Surv(OS_time.x,OS_statue) ~sur.cat1[,i], data = sur.cat1) 
  pValue = surv_pvalue(fit)$pval
  x <- c(colnames(sur.cat1)[i],pValue)
  log_rank_res1 <- rbind(log_rank_res1,x)
  print(i)
}
colnames(log_rank_res1) <- c("mRNA","p")
sfit <- surv_fit(Surv(OS_time.x,OS_statue)~ELAPOR1,data = sur.cat1)
ggsurvplot(sfit,
           xlab=("Times(months)"),
           legend = c(0.8,0.75),
           legend.labs=c("High","Low"), 
           legend.title="Classifaction",
           xlim=c(0,60),
           break.time.by = 12,
           pval = TRUE, 
           pval.size=6,
           risk.table = TRUE, 
           risk.table.col = "black", 
           surv.median.line = "hv",l
           ggtheme = theme(axis.title = element_text(face = "bold",size = 20), 
                           axis.text  = element_text(face = "bold",size = 13), 
                           plot.title = element_text(size = 20,hjust = 0.5,face = "bold"),
                           legend.title = element_text(size = 15,face = "bold",), 
                           legend.text = element_text(size = 15,face = "bold",),
                           panel.background = element_blank(),
                           axis.line.x = element_line(linewidth=1, colour = "black"),
                           axis.line.y = element_line(linewidth=1, colour = "black"),
                           text=element_text(family="serif")
           ), 
           tables.theme = theme(axis.title = element_text(face = "bold",size = 15), 
                                axis.text  = element_text(face = "bold",size = 10), 
                                axis.line = element_line(linewidth=1, colour = "black"),
                                plot.title = element_text(size = 15,hjust = 0.5,face = "bold") 
                                ),
           title="ELAPOR1 Survival Analysis"
)
gene_last_exp <- TCGA_PAAD_mRNA_FPKM[gene_last,]
CPB1_exp <- TCGA_PAAD_mRNA_FPKM["CPB1",]
CPB1_exp <- data.frame(t(CPB1_exp))
CPB1_exp$lable <- ifelse(CPB1_exp$CPB1>=mean(CPB1_exp$CPB1),"high","low")
CPB1_exp$sample_id <- rownames(CPB1_exp)
high_sample <- CPB1_exp[which(CPB1_exp$lable=="high"),]
low_sample <- CPB1_exp[which(CPB1_exp$lable=="low"),]
gene_last_exp <- gene_last_exp[,c(high_sample$sample_id,low_sample$sample_id)]

cls <- c(rep("high",length(high_sample$sample_id)),rep("low",length(low_sample$sample_id)))
library(ggrepel)
library(stringr)
GO_BP <- read.csv("GO-BP.csv",header = T)
GO_CC <- read.csv("GO-CC.csv",header = T)
GO_MF <- read.csv("GO-MF.csv",header = T)
GO_BP1 <- GO_BP[1:5,]
GO_CC1 <- GO_CC[1:5,]
GO_MF1 <- GO_MF[1:5,]
GO <- rbind(GO_BP1,GO_CC1,GO_MF1)
GO$lable <- c(rep("BP",5),
              rep("CC",5),
              rep("MF",5))
ggplot(data = GO, 
       aes(x = Counts,y = reorder(Description,Counts)))+ 
  geom_point(aes(color = -LogP),size=10)+ 
  facet_wrap(~lable,
             nrow=3,
             scales = "free_y",
             strip.position = "right")+
  theme_bw()+
  scale_colour_gradient(low = "red",high = "purple")+ 
       color = expression(-LogP),size="count")+ 
  scale_x_continuous(breaks=seq(3,11,1))+
  theme(axis.title = element_text(size = 20,face = "bold"),
        axis.text = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 11,face = "bold"),
        panel.border = element_rect(fill = NA,linewidth=2, colour = "black"),
        axis.ticks.length = unit(8,"pt"),
        axis.ticks = element_line(linewidth = 1.2),
        text=element_text(family="serif"),
        strip.background = element_rect(fill = "grey",
                                        linetype = "dotted"),
        strip.text = element_text(color = "black",
                                  face = "bold",
                                  size = 20),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
KEGG <- read.csv("放疗样本候选基因KEGG富集分析结果.csv",header = T)
ggplot(data = KEGG, 
       aes(x = Description, y = Counts,fill = c("blue")))+ 
  geom_bar(stat = "identity",width = 0.9)+ 
  coord_flip()+
  theme_bw()+ 
  scale_y_continuous(breaks=seq(0,12,2))+
  labs(y = "GeneNumber")+ 
  theme(axis.title = element_text(size = 13,face = "bold"), 
        axis.text = element_text(size = 11,face = "bold"), 
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,linewidth=2, colour = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = "none",
        text=element_text(family="serif"))
target_exp <- data.frame(t(TCGA_PAAD_mRNA_counts["CPB1",c(normal_sample_ID,tumor_sample_ID)]))
target_exp$lable <- c(rep("normal",4),rep("tumor",179))
target_exp$sample <- substr(rownames(target_exp),1,16)
target_exp$CPB1 <- log2(target_exp$CPB1)
target_exp <- target_exp[order(target_exp$CPB1),]
target_exp <- target_exp[-c(1,2),]
library("ggplot2")
ggplot(target_exp,aes(x=lable,y=CPB1,fill=lable))+
  stat_boxplot(geom="errorbar",width=0.4,size=1,position=position_dodge(0.8),color="black")+
  geom_boxplot(width=0.6,
               position = position_dodge(0.8))+
  theme_bw()+
  scale_fill_manual(values = c("#00bfc4","#f8766d"))+
  scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim=c(-1,25))
  )


