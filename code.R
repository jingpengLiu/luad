#######################FUNCTION################
suriv_function <- function(x){
  if(x == "Dead"){
    return(1)
  }else{
    return(0)
  }
}

gender_ismale <- function(x){
  if(x == "male"){
    return(1)
  }else{
    return(0)
  }
}
M_is <- function(x){
  if(x == "M1" | x == "M1a"  | x == "M1b"){
    return(1)
  }else if(x == "M0"){
    return(0)
  }else if(x == "MX"){
    return(2)
  }
}
N_is <- function(x){
  if(x == "N1"){
    return(1)
  }else if(x == "N0"){
    return(0)
  }else if(x == "NX"){
    return(4)
  }else if(x == "N2"){
    return(2)
  }else if(x=="N3"){
    return(3)
  }
}
T_is <- function(x){
  if(x == "T1" | x == "T1a" | x == "T1b"){
    return(1)
  }else if(x == "T2" | x == "T2a" | x == "T2b"){
    return(2)
  }else if(x == "T4"){
    return(4)
  }else if(x == "TX"){
    return(5)
  }else if(x=="T3"){
    return(3)
  }
}
stage_is <- function(x){
  if(x == "stage i" | x == "stage ia" | x == "stage ib"){
    return(1)
  }else if(x == "stage ii" | x == "stage iia" | x == "stage iib"){
    return(2)
  }else if(x == "stage iiia" | x == "stage iiib"){
    return(3)
  }else if(x == "stage iv"){
    return(4)
  }
}


######################一、收集TCGA-LUAD数据并处理###############
library(dplyr)
library(stringr)
library(data.table)
###--------------------设置数据文件的路径--------------------
filesNameToBarcode <- read.table(file = "./data/gdc_sample_sheet.2023-05-11.tsv",
                                 header = TRUE, sep="\t", check.names = F,
                                 quote="")
filesNameToBarcode$filepath <- paste0('./data/download/da916720c149404d8be4fdd5fb6c66c8/',
                                      filesNameToBarcode$`File ID`,
                                      '/',
                                      filesNameToBarcode$`File Name`)
###--------------------读入转录组表达数据--------------------
for(i in 1:nrow(filesNameToBarcode)){
  file_name <- filesNameToBarcode$`File Name`[i]
  sample_name <- filesNameToBarcode$`Sample ID`[i]
  sample_type <- filesNameToBarcode$`Sample Type`[i]
  print(paste0("提示:正在读入文件:",file_name,"----",sample_name,'---',sample_type))
  
  oneSampExp <- fread(filesNameToBarcode$filepath[i],header =T)
  
  oneSampCounts <- oneSampExp[-c(1:4),c(1,2,3,4)]
  oneSampTPM <- oneSampExp[-c(1:4),c(1,2,3,7)]
  oneSampFPKM <- oneSampExp[-c(1:4),c(1,2,3,8)]
  
  colnames(oneSampCounts)[4] <- sample_name
  colnames(oneSampTPM)[4] <- sample_name
  colnames(oneSampFPKM)[4] <- sample_name
  if (i == 1){
    RNAseq_Counts <- oneSampCounts
    RNAseq_FPKM <- oneSampFPKM
    RNAseq_TPM <- oneSampTPM
  }
  else{
    RNAseq_Counts <- left_join(RNAseq_Counts,oneSampCounts)
    RNAseq_FPKM <- left_join(RNAseq_FPKM,oneSampFPKM)
    RNAseq_TPM <- left_join(RNAseq_TPM,oneSampTPM)
  }
}
###--------------------重复基因取最大值--------------------
max_gene <- function(exp){
  index=order(rowMeans(exp[,-c(1,2,3)]),decreasing = T)
  expr_ordered=exp[index,]
  keep=!duplicated(expr_ordered$gene_name)
  expr_max=expr_ordered[keep,]
  return(expr_max)
}

RNAseq_Counts <- max_gene(RNAseq_Counts)
RNAseq_FPKM <- max_gene(RNAseq_FPKM)
RNAseq_TPM <- max_gene(RNAseq_TPM)

#根据行的和来过滤，也可以根据行的平均值来过滤
RNAseq_Counts <- RNAseq_Counts[rowSums(RNAseq_Counts[,-c(1:3)])>1,]
RNAseq_FPKM <- RNAseq_FPKM[rowSums(RNAseq_FPKM[,-c(1:3)])>1,]
RNAseq_TPM <- RNAseq_TPM[rowSums(RNAseq_TPM[,-c(1:3)])>1,]


###--------------------提取normal和tumor的表达数据--------------------
Normal_sample <- filesNameToBarcode$`Sample ID`[grepl(pattern = 'Normal',x = filesNameToBarcode$`Sample Type`)]

RNAseq_Counts <- as.data.frame(RNAseq_Counts)
RNAseq_count_normal <- RNAseq_Counts[,c(1:3,which(colnames(RNAseq_Counts) %in% Normal_sample))]
RNAseq_count_tumor <- RNAseq_Counts[,-which(colnames(RNAseq_Counts) %in% Normal_sample)]

RNAseq_FPKM <- as.data.frame(RNAseq_FPKM)
RNAseq_FPKM_normal <- RNAseq_FPKM[,c(1:3,which(colnames(RNAseq_FPKM) %in%  Normal_sample))]
RNAseq_FPKM_tumor <- RNAseq_FPKM[,-which(colnames(RNAseq_FPKM) %in%  Normal_sample)]

RNAseq_TPM <- as.data.frame(RNAseq_TPM)
RNAseq_TPM_normal <- RNAseq_TPM[,c(1:3,which(colnames(RNAseq_TPM) %in%  Normal_sample))]
RNAseq_TPM_tumor <- RNAseq_TPM[,-which(colnames(RNAseq_TPM) %in%  Normal_sample)]

###---------------------记录样本性别数据-------------------
gender <- read.csv("./data/clinical.cart.2023-05-11/clinical.tsv", header=T, sep="\t")
male <- gender[gender$gender == "male", 2]
male <- paste(male, "-01A", sep="")
female <- gender[gender$gender == "female", 2]
female <- paste(female, "-01A", sep="")
###--------------------保存数据--------------------
save(RNAseq_Counts,RNAseq_FPKM,RNAseq_TPM,
     RNAseq_count_normal,RNAseq_count_tumor,
     RNAseq_FPKM_normal,RNAseq_FPKM_tumor,
     RNAseq_TPM_normal,RNAseq_TPM_tumor,
     male,female,file = 'TCGA-BLCA-RNAseq.Rdata')

###--------------------加载数据--------------------
load('TCGA-BLCA-RNAseq.Rdata')




######################二、识别与细胞周期相关的差异表达基因###############
load('TCGA-BLCA-RNAseq.Rdata')
# 读取转录组数据
RNAseq_Counts <- cbind(RNAseq_count_normal, RNAseq_count_tumor[,c(4:ncol(RNAseq_count_tumor))])
rna_seq <- RNAseq_Counts[,-c(1,3)]
rownames(rna_seq) <- rna_seq$gene_name
rna_seq <- rna_seq[,-1]
condition <- factor(c(rep("Normal",length(colnames(RNAseq_count_normal))-3), 
                      rep("Tumor",length(colnames(RNAseq_count_tumor))-3))) #
coldata <- data.frame(row.names = colnames(rna_seq), condition)
####DEseq差异分析
# 读取并且处理数据
library(DESeq2)
# 差异分析
dds <- DESeqDataSetFromMatrix(countData = rna_seq, colData = coldata,
                              design = ~condition)
# 过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 4 #
dds <- dds[keep, ]
# 利用deseq函数标准化dds矩阵
dep <- DESeq(dds)
res <- results(dep)
res <- res[order(res$padj),]
diff <- res
diff <- na.omit(diff)
dim(diff)
write.csv(diff, './data/all_diff.csv') #

##############3、
diff <- read.csv("./data/all_diff.csv")

# 差异基因和周期基因交集火山图
library(readxl)
library(ggplot2)
diff <- read.csv("./data/all_diff.csv")
zhouqi <- read_excel("./data/GeneCards-细胞周期相关评分大于30.xlsx")
colnames(diff)[1] <- "id"
colnames(zhouqi)[1] <- "id"
diff_zhouqi <- merge(diff, zhouqi, by="id")
cut_off_pvalue = 0.05
cut_off_logFC = 1
diff_zhouqi$Sig = ifelse(diff_zhouqi$pvalue < cut_off_pvalue & 
                           abs(diff_zhouqi$log2FoldChange) >= cut_off_logFC, 
                         ifelse(diff_zhouqi$log2FoldChange> cut_off_logFC ,'Up','Down'),'None')
Up_num <- sum(diff_zhouqi$Sig == 'Up')
Down_num <- sum(diff_zhouqi$Sig == 'Down')
#绘制火山图
ggplot(diff_zhouqi,aes(log2FoldChange, -log10(padj)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=-log10(padj), color= -log10(padj)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+
  theme(panel.grid = element_blank()) +
  # 注释
  # 注释
  annotate("text", label = "bolditalic(Down)", parse = TRUE, 
           x = -3.5, y = 335, size = 4, colour = "black")+
  annotate("text", label = "bolditalic(Up)", parse = TRUE, 
           x = 3.0, y = 335, size = 4, colour = "black")+
  annotate("text", label = Down_num, parse = TRUE, 
           x = -3.5, y = 320, size = 3, colour = "black")+
  annotate("text", label = Up_num, parse = TRUE, 
           x = 3.0, y = 320, size = 3, colour = "black")
ggsave("./picture/volco.pdf", width = 8, height = 8)

#绘制基因表达量热图
express_data <- RNAseq_TPM[,c(-1,-3)]
gene_20 <- diff_zhouqi[order(abs(diff_zhouqi$log2FoldChange), decreasing = T),]#取前20的差异基因
gene_20 <- gene_20$id[1:20]
express_data <- express_data[sapply(express_data$gene_name, function(x) x %in% gene_20),]
Normal_sample <- colnames(RNAseq_count_normal)[4:length(RNAseq_count_normal)]
Tumor_sample <- colnames(RNAseq_count_tumor)[4:length(RNAseq_count_tumor)]
#plot
library(pheatmap)
df <- as.matrix(express_data)
row.names(df) <- df[,1]
df <- df[,-1]
dfs <- apply(df,2,function(x) as.numeric(x))
row.names(dfs) <- row.names(df)
annotation_col <- data.frame(
  row.names = c(Normal_sample, Tumor_sample),
  Type = c(rep("Normal",length(Normal_sample)), rep("Tumor", length(Tumor_sample)))
)
dfs <- dfs[,row.names(annotation_col)]

bk = c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))
color = c(colorRampPalette(colors = c("#8854d0","#ffffff"))(length(bk)/2),colorRampPalette(colors = c("#ffffff","#fa8231"))(length(bk)/2))
pheatmap(dfs,
         show_rownames = T,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows=T,
         filename='./picture/express.pdf',#输出文件的名称
         fontsize_row=6, #行字体的大小
         height=10,  #输出图片的高度
         weight=10,
         scale = "row",
         angle_col=45, #调整字的角度
         color = color,
         clustering_distance_rows = 'euclidean', 
         clustering_method = 'single',
         annotation_col = annotation_col,
         annotation_names_col = TRUE,
         breaks = bk
)
#GO和KEGG功能富集分析############
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
gene.df <- bitr(diff_zhouqi[diff_zhouqi$Sig == "Up" | diff_zhouqi$Sig == "Down",]$id,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene <- gene.df$ENTREZID                
GO <- enrichGO( gene,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pvalueCutoff = 0.05,#设定p值阈值
                minGSSize = 1,
                qvalueCutoff = 0.05,
                pAdjustMethod = 'BH',
                readable= FALSE)
GO_df <- as.data.frame(GO)
write.csv(GO_df, "./data/GO_df.csv")
KEGG <- enrichKEGG(gene, 
                   organism = 'hsa',  ## hsa为人的简写
                   keyType = 'kegg', 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH', 
                   minGSSize = 10,
                   maxGSSize = 500,
                   qvalueCutoff = 0.05,
                   use_internal_data = FALSE)
KEGG_df <- as.data.frame(KEGG)
write.csv(KEGG_df, "./data/KEGG_df.csv")
#绘图
pdf(file="./picture/go_bar.pdf",width=10,height=15)
barplot(GO, split="ONTOLOGY", showCategory =10)+facet_grid(ONTOLOGY~., scale="free")#柱状图
dev.off()
pdf(file="./picture/kegg_bar.pdf",width=10,height=5)
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')
dev.off()
#记录
write.csv(gene.df,"./data/gene_entriz.csv")




######################三、识别与细胞周期相关的差异表达基因###############
#单因素Cox回归分析之数据处理
siglecox_data <- RNAseq_TPM[sapply(RNAseq_TPM$gene_name, function(x) x %in% gene.df$SYMBOL),-c(1,3)]
filesNameToBarcode <- read.table(file = "./data/gdc_sample_sheet.2023-05-11.tsv",
                                 header = TRUE, sep="\t", check.names = F,
                                 quote="")
filesNameToBarcode <- filesNameToBarcode[,c(6,7)]
gender <- read.csv("./data/clinical.cart.2023-05-11/clinical.tsv", header=T, sep="\t")
colnames(gender)[2] <- colnames(filesNameToBarcode)[1]
gender <- merge(gender, filesNameToBarcode, by="Case ID")
gender <- gender[,c(159,1:158)]
gender <- gender[!duplicated(gender$`Sample ID`),]
row.names(siglecox_data) <- siglecox_data$gene_name
siglecox_data <- siglecox_data[,-1]
siglecox_data <- t(siglecox_data)
gender_live <- gender[,c(1, 17, 11, 51)]
x <- c()
for(i in (1:nrow(gender_live))){
  if(gender_live[i,3] == "'--"){
    x = append(x,gender_live[i,4])
  }else{
    x = append(x,gender_live[i,3])
  }
}
gender_live$age_at_diagnosis <- x
gender_live <- gender_live[,-c(3,4)]
gender_live$surtime <- sapply(gender_live$vital_status, function(x) suriv_function(x))
siglecox_data <- as.data.frame(siglecox_data)
siglecox_data$ID <- row.names(siglecox_data)
colnames(gender_live)[1] <- "ID"
siglecox_data <- merge(siglecox_data, gender_live, by="ID") #得到带有生存时间和是否存活的整合数据
siglecox_data <- siglecox_data[,c(1, 115, 117, 116,2:114)]
row.names(siglecox_data) <- siglecox_data$ID
siglecox_data <- siglecox_data[,-c(1,2)]
colnames(siglecox_data)[1:2] <- c("surstat", "surtime")
siglecox_data <- siglecox_data[siglecox_data$surtime != "'--" , ]#去除为na的样本数据
siglecox_data[,2] <- as.numeric(siglecox_data[,2]) 
siglecox_data <- siglecox_data[siglecox_data$surtime != 0,]
#开始单因素分析
library("survival")
library("survminer")
td <- siglecox_data
pFilter=0.05 #设一个p值标准，后面用
outResult=data.frame() #建一个空白数据框，后面for循环输出用
sigGenes=c("surstat","surtime") #建一个向量，后面for循环输出用，因为后面还要用到surstat及surtime，所以先放在向量里
for(i in colnames(td[,3:ncol(td)])){ #从第3列开始循环，因为1列2列不是gene，是surstat和surtime
  tdcox <- coxph(Surv(surtime, surstat) ~ td[,i], data = td)#开始逐一循环cox分析
  tdcoxSummary = summary(tdcox) #summary命令对tdcox总结，方面后面提取数据
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #提取p值，这个主要是后面提取有意义的gene用
  if(pvalue<pFilter){ # 这里我们不需要提取所有基因的数据，只需要有意义的gene的结果，所以设置如果pvalue<0.05才提取数据
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,#合并行，实际上是对循环结果的合并，前面设置的空白数据框outResult这里用，循环必须有个开始
                    cbind(id=i,#合并列，是每个基因的统计数据
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],#提取单个基因的HR
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],#提取单个基因的HR的95%CI低值
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],#提取单个基因的HR的95%CI高值
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#提取单个基因的p值
    )
  }
}
write.table(outResult,file="./data/UniCoxSurvival.csv")
write.table(outResult,file="./data/UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)
UniCoxSurSigGeneExp=td[,sigGenes] #以有统计学意义的gene取表达数据的子集
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)#以id也就是样品名命名行名
write.table(UniCoxSurSigGeneExp,file="./data/UniCoxSurSigGeneExp.csv")

#森林图
tducs <- read.table("./data/UniCoxSurvival.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(tducs)
hr <- sprintf("%.3f",tducs$"HR")
hrLow  <- sprintf("%.3f",tducs$"L95CI")
hrHigh <- sprintf("%.3f",tducs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))

pdf(file="./picture/UniCoxSurForestPlot.pdf", width = 6,height = 10)
n <- nrow(tducs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))

xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex);text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)
dev.off()




######################四五六八、LASSO+COX回归建立细胞周期相关基因的预后模型###############
######################四五六八、基于预后模型进行风险评分，将患者分为底危组和高危组。###############
######################四五六八、KM和ROC曲线分析两组的评估预后特征###############
######################四五六八、验证集进行验证###############
#####lasso+cox多因素回归分析筛选基因
library(glmnet)
library(survival)
x = UniCoxSurSigGeneExp[,c(4:ncol(UniCoxSurSigGeneExp))]
x = as.matrix(x)
y = data.matrix(Surv(time=UniCoxSurSigGeneExp$surtime, event=UniCoxSurSigGeneExp$surstat))
fit <- glmnet(x=UniCoxSurSigGeneExp[,c(4:ncol(UniCoxSurSigGeneExp))], y, family="cox", alpha=1)
pdf("./picture/lassocox.pdf", height=5, width=7)
plot(fit, xvar="lambda", label=TRUE)
dev.off()
lasso_fit <- cv.glmnet(x, y, family = "cox", type.measure = "deviance")
pdf("./picture/lasso_fit.pdf", height=5, width=7)
plot(lasso_fit, label=T)
dev.off()

coefficient <- coef(lasso_fit, 0.004)
Active.Index <- which(as.numeric(coefficient)!=0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
sig_gene_multi_cox

#基于预后模型进行风险评分，将患者分为高危组和底危组
trainFinalGeneExp <- UniCoxSurSigGeneExp[,sig_gene_multi_cox]
trainScore <- predict(lasso_fit, x, s=0.002, type="response")
outCol <- c("surtime", "surstat", sig_gene_multi_cox)
risk <- as.vector(ifelse(trainScore>median(trainScore),"high","low"))
train <- cbind(UniCoxSurSigGeneExp[,outCol], riskScore=as.vector(trainScore), risk)
names(trainScore) <- rownames(train)
write.csv(train,"./data/train.csv")
#绘图
fp <- trainScore
phe <- train
fp_dat <- data.frame(patientid=1:length(fp), fp=as.numeric(sort(fp)))
fp_dat$riskgroup = ifelse(fp_dat$fp>=median(fp_dat$fp), "high", "low")

sur_dat = data.frame(patientid=1:length(fp), time=phe[names(sort(fp)),'surtime'], event=phe[names(sort(fp)),'surstat']) 
sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
exp_dat=phe[names(sort(fp)),3:(ncol(phe)-2)]
###第一个图
library(ggplot2)
p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
  scale_colour_manual(values = c("red","green"))+
  theme_bw()+labs(x="Patient ID(increasing risk score)",y="Risk score")+
  geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
ggsave("./picture/risk1.pdf",height=6,width=10)
#第二个图
p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
  scale_colour_manual(values = c("red","green"))+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
ggsave("./picture/risk2.pdf",height=6,width=10)
#第三个图
library(pheatmap)
mycolors <- colorRampPalette(c("white", "green", "red"), bias = 1.2)(100)
tmp=t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
annotation_col <- data.frame(
  row.names = row.names(train),
  type = train$risk,
  test = c(1:nrow(train))
)
annotation_col <- annotation_col[colnames(tmp),]
annotation_col <- data.frame(
  row.names = row.names(annotation_col),
  type = annotation_col$type
)
pdf("./picture/risk3ss.pdf", height=6, width=10)
p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,
            annotation_col = annotation_col)
dev.off()

#拼图实现三图联动
library(ggplotify)
plots = list(p1,p2,as.ggplot(as.grob(p3)))
library(gridExtra)
lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7))) #布局矩阵
pdf("./picture/risk_all.pdf", height=20, width=10)
grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
dev.off()



#KM曲线---------------------------------------------------------------------
library("survival")
library("survminer")
fit <- survfit(Surv(surtime, surstat) ~ risk, data = train)
#按分层更改图形颜色，线型等
pdf('./picture/train_KM.pdf',w=7,h=6,onefile = F)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,legend.labs=c("High","Low" ),legend.title="Risk", title="TCGA-LUAD",
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           palette = "npg",
           font.x=15,          
           font.y=15,
           font.title=15,
           font.subtitle = 14,
           font.caption = 14,
           font.tickslab = 14,
           font.legend=17,     
           fontsize = 5,      
)
dev.off()

# ROC---------------------------------------------------------------------------
library(survivalROC)
ROC <- train
ROC$futime <- round(ROC$surtime)
ROC$fustat <- ROC$surstat
ROC$riskscore <- ROC$risk
cutoff_1 <- 365*1 
cutoff_3 <- 365*2
cutoff_5 <- 365*3

year_1= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_1,
                    method = 'KM')
year_3= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_3,
                    method = 'KM')
year_5= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_5,
                    method = 'KM')
if(T){
  for(i in 1:7){
    cutoff <- 365*i
    year= survivalROC(Stime=ROC$futime,
                      status=ROC$fustat,
                      marker = ROC$riskscore, 
                      predict.time = cutoff,
                      method = 'KM')
  }
}

pdf("./picture/01.train_ROC.pdf",w=7,h=7)
par(mar=c(6,6,4,3))
plot(year_1$FP, year_1$TP,
     type="l",col="red",xlim=c(0,1), ylim=c(0,1),
     xlab="FP",
     ylab="TP",
     main="TCGA-LUAD, Method=KM\n Year = 1,2,3",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.5
)
abline(0,1,col="gray",lty=2)
lines(year_3$FP, year_3$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1))
lines(year_5$FP, year_5$TP, type="l",col="blue",xlim=c(0,1), ylim=c(0,1))
legend(0.5,0.2,c(paste("AUC of 1 year =",round(year_1$AUC,3)),
                 paste("AUC of 2 year =",round(year_3$AUC,3)),
                 paste("AUC of 3 year =",round(year_5$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("red","green",'blue'),
       bty = "n",
       seg.len=1,cex=1.15)#
dev.off()



############################使用GSE68465验证集做验证#############
test_gse68465 <- read.table("./data/07.expr_GSE68465.txt", header=T, sep="\t", row.names= 1)
test_surv <- read.table("./data/08.survival_GSE68465.txt", header=T)
test_gse68465 <- t(test_gse68465)
test_gene <- intersect(colnames(x), colnames(test_gse68465))
#test_gene <- c(test_gene[1:10],"CDCA4",test_gene[11:19],"FANCC",test_gene[20:23])
test_gse68465 <- test_gse68465[,test_gene]
#colnames(test_gse68465) <- colnames(x)
test_gse68465 <- as.data.frame(test_gse68465)
test_gse68465$futime <- test_surv$futime
test_gse68465$fustat <- test_surv$fustat
test_gse68465$surtime <- test_surv$futime
test_gse68465$surstat <- test_surv$fustat
test <- test_gse68465



#基于预后模型进行风险评分，将患者分为高危组和底危组
x = as.matrix(test_gse68465[1:25])
trainScore <- predict(lasso_fit, x, s=0.002, type="response")
outCol <- c("surtime", "surstat", sig_gene_multi_cox)
risk <- as.vector(ifelse(trainScore>median(trainScore),"high","low"))
train <- cbind(test[,outCol], riskScore=as.vector(trainScore), risk)
names(trainScore) <- rownames(train)
write.csv(train,"./data/test.csv")
#绘图
fp <- trainScore
phe <- train
fp_dat <- data.frame(patientid=1:length(fp), fp=as.numeric(sort(fp)))
fp_dat$riskgroup = ifelse(fp_dat$fp>=median(fp_dat$fp), "high", "low")

sur_dat = data.frame(patientid=1:length(fp), time=phe[names(sort(fp)),'surtime'], event=phe[names(sort(fp)),'surstat']) 
sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
exp_dat=phe[names(sort(fp)),3:(ncol(phe)-2)]
###第一个图
library(ggplot2)
p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
  scale_colour_manual(values = c("red","green"))+
  theme_bw()+labs(x="Patient ID(increasing risk score)",y="Risk score")+
  geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
ggsave("./picture/test_risk1.pdf",height=6,width=10)
#第二个图
p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
  scale_colour_manual(values = c("red","green"))+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
ggsave("./picture/test_risk2.pdf",height=6,width=10)
#第三个图
library(pheatmap)
mycolors <- colorRampPalette(c("white", "green", "red"), bias = 1.2)(100)
tmp=t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
annotation_col <- data.frame(
  row.names = row.names(train),
  type = train$risk,
  test = c(1:nrow(train))
)
annotation_col <- annotation_col[colnames(tmp),]
annotation_col <- data.frame(
  row.names = row.names(annotation_col),
  type = annotation_col$type
)
pdf("./picture/test_risk3ss.pdf", height=6, width=10)
p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,
            annotation_col = annotation_col)
dev.off()

#拼图实现三图联动
library(ggplotify)
plots = list(p1,p2,as.ggplot(as.grob(p3)))
library(gridExtra)
lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7))) #布局矩阵
pdf("./picture/test_risk_all.pdf", height=20, width=10)
grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
dev.off()



#KM曲线---------------------------------------------------------------------
library("survival")
library("survminer")
fit <- survfit(Surv(surtime, surstat) ~ risk, data = train)
#按分层更改图形颜色，线型等
pdf('./picture/test_KM.pdf',w=7,h=6,onefile = F)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,legend.labs=c("High","Low" ),legend.title="Risk", title="TCGA-LUAD",
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           palette = "npg",
           font.x=15,          
           font.y=15,
           font.title=15,
           font.subtitle = 14,
           font.caption = 14,
           font.tickslab = 14,
           font.legend=17,     
           fontsize = 5,      
)
dev.off()

# ROC---------------------------------------------------------------------------
library(survivalROC)
ROC <- train
ROC$futime <- round(ROC$surtime)
ROC$fustat <- ROC$surstat
ROC$riskscore <- ROC$risk
cutoff_1 <- 365*1 
cutoff_3 <- 365*2
cutoff_5 <- 365*3

year_1= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_1,
                    method = 'KM')
year_3= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_3,
                    method = 'KM')
year_5= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_5,
                    method = 'KM')
if(T){
  for(i in 1:7){
    cutoff <- 365*i
    year= survivalROC(Stime=ROC$futime,
                      status=ROC$fustat,
                      marker = ROC$riskscore, 
                      predict.time = cutoff,
                      method = 'KM')
  }
}

pdf("./picture/test_01.train_ROC.pdf",w=7,h=7)
par(mar=c(6,6,4,3))
plot(year_1$FP, year_1$TP,
     type="l",col="red",xlim=c(0,1), ylim=c(0,1),
     xlab="FP",
     ylab="TP",
     main="TCGA-LUAD, Method=KM\n Year = 1,2,3",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.5
)
abline(0,1,col="gray",lty=2)
lines(year_3$FP, year_3$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1))
lines(year_5$FP, year_5$TP, type="l",col="blue",xlim=c(0,1), ylim=c(0,1))
legend(0.5,0.2,c(paste("AUC of 1 year =",round(year_1$AUC,3)),
                 paste("AUC of 2 year =",round(year_3$AUC,3)),
                 paste("AUC of 3 year =",round(year_5$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("red","green",'blue'),
       bty = "n",
       seg.len=1,cex=1.15)#
dev.off()
#1、预后特征和临床之间的关系
library(reshape2)
cdc20_stat <- data.frame(
  gene = "CDC20",
  group = sapply(train$surstat, function(x) ifelse(x==0,"r1","r2")),
  values =  train$CDC20,
  group2 = train$surstat,
  cancer = "Stat"
)
cdc20_risk <- data.frame(
  gene = "CDC20",
  group = sapply(train$risk, function(x) ifelse(x=="high","r3","r4")),
  values =  train$CDC20,
  group2 = train$risk,
  cancer = "Risk"
)
cdc20 <- rbind(cdc20_stat, cdc20_risk)


# 绘图：
#是否生存
ggplot(cdc20[cdc20$cancer=="Stat",],aes(x=group2, y=values))+
  geom_boxplot(outlier.shape = NA, fill = '#A4A4A4', color = "black")+
  # x轴标签：
  scale_x_discrete(labels = c("Alive","Dead"))+
  # x轴和y轴标签
  xlab("fustat")+
  ylab("Gene expression")+
  # 主题：
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
  ggtitle("CDC20")+
  geom_signif(y_position=c(250, 8.5),comparisons =  list(c("0", "1")))

ggsave("./picture/cdc20_sur.pdf", height=5, width=5)
#患者的高低风险
ggplot(cdc20[cdc20$cancer=="Risk",],aes(x=group2, y=values))+
  geom_boxplot(outlier.shape = NA, fill = '#A4A4A4', color = "black")+
  # x轴标签：
  scale_x_discrete(labels = c("High","Low"))+
  # x轴和y轴标签
  xlab("Risk")+
  ylab("Gene expression")+
  # 主题：
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
  ggtitle("CDC20")+
  geom_signif(y_position=c(250, 8.5),comparisons =  list(c("high", "low")))

ggsave("./picture/cdc20_risk.pdf", height=5, width=5)
#患者的性别
gender_df <- data.frame(
  id = row.names(train),
  sid = substr(row.names(train),1,12)
)
gender_1 <- data.frame(
  sid = gender$case_submitter_id,
  gender = gender$gender
)
gender_df <- merge(gender_df, gender_1, by="sid", all.x=T)
gender_df <- gender_df[!duplicated(gender_df$id),]
train$gender <- sapply(row.names(train),function(x) gender_df[gender_df$id == x,3])

gender_pl <- data.frame(
  gene = "cdc20",
  values = train$CDC20,
  group2 = train$gender
)
ggplot(gender_pl,aes(x=group2, y=values))+
  geom_boxplot(outlier.shape = NA, fill = '#A4A4A4', color = "black")+
  # x轴标签：
  scale_x_discrete(labels = c("Female","Male"))+
  # x轴和y轴标签
  xlab("Gender")+
  ylab("Gene expression")+
  # 主题：
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
  ggtitle("CDC20")+
  geom_signif(y_position=c(250, 8.5),comparisons =  list(c("female", "male")))

ggsave("./picture/cdc20_gender.pdf", height=5, width=5)
#肿瘤和非肿瘤样本
train$normal <- sapply(row.names(train),function(x) filesNameToBarcode[filesNameToBarcode$`Sample ID` == x,8][1])
train$normal <- sapply(train$normal, function(x) ifelse(grepl("Tumor", x),"Tumor","Normal"))
normal_pl <- data.frame(
  gene = "cdc20",
  values = train$CDC20,
  group2 = train$normal
)
ggplot(normal_pl,aes(x=group2, y=values))+
  geom_boxplot(outlier.shape = NA, fill = '#A4A4A4', color = "black")+
  # x轴标签：
  scale_x_discrete(labels = c("Normal","Tumor"))+
  # x轴和y轴标签
  xlab("Type")+
  ylab("Gene expression")+
  # 主题：
  theme_bw()+
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
  ggtitle("CDC20")+
  geom_signif(y_position=c(250, 8.5),comparisons =  list(c("Normal", "Tumor")))

ggsave("./picture/cdc20_normal.pdf", height=5, width=5)
#2、通过单因素Cox回归分析和多因素Cox回归分析，评估预后模型与年龄、性别、组织学分、病理分期和TNM分期等常见的临床因素
library(survival)
train <- read.csv("./data/train.csv", row.names = 1)
train <- train[,c(1,2,25,26)]
gender <- read.table("./data/xianyu/05.clinical_tcga_luad.txt", sep="\t", header=T, row.names = 1)
TNM_cox <- merge(train, gender, by="row.names")
TNM_cox <- TNM_cox[TNM_cox$pathologic_M != "" & TNM_cox$pathologic_N != "" & TNM_cox$pathologic_T != "" & TNM_cox$Stage != "",]
TNM_cox <- TNM_cox[TNM_cox$Stage != "not reported",]
TNM_cox$math_Gender <- as.vector(sapply(TNM_cox$Gender, function(x) gender_ismale(x))) #男性为1女性为0
TNM_cox$math_pathologic_M <- as.vector(sapply(TNM_cox$pathologic_M, function(x) M_is(x)))
TNM_cox$math_pathologic_N <- as.vector(sapply(TNM_cox$pathologic_N, function(x) N_is(x)))
TNM_cox$math_pathologic_T <- as.vector(sapply(TNM_cox$pathologic_T, function(x) T_is(x)))
TNM_cox$math_Stage <- as.vector(sapply(TNM_cox$Stage, function(x) stage_is(x)))

tdmultiCox <- coxph(Surv(surtime, surstat) ~ (Age + math_Gender + math_Stage + math_pathologic_T + math_pathologic_M + math_pathologic_N + riskScore), data = TNM_cox) #这里有个“.”，代表分析td数据中所有变量（列名）
tdmultiCoxSum=summary(tdmultiCox)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
# 可视化
library(survminer) #载入所需R包
ggforest(tdmultiCox, #直接用前面多因素cox回归分析的结果
         main = "Hazard ratio",
         cpositions = c(0.02,-0.15, 0.25), #前三列的位置，第二列是样品数，设了个负值，相当于隐藏了
         fontsize = 0.8, #字体大小
         refLabel = "reference", 
         noDigits = 2)
ggsave("./picture/snm_mutil_cox.pdf", width=8, height=6)
# 5、gsea富集分析
gsea_diff <- read.csv("./data/all_diff.csv", row.names = 1)
library(clusterProfiler)
# 基因id转换
library(org.Hs.eg.db)
entriz <- mapIds(org.Hs.eg.db, keys = row.names(gsea_diff), keytype = "SYMBOL", column="ENTREZID")
entriz <- as.data.frame(entriz)
write.csv(entriz, "./data/entriz.csv")
gsea_entriz_df <- merge(entriz, gsea_diff, by="row.names")
gsea_entriz_df <- gsea_entriz_df[is.na(gsea_entriz_df$entriz) == FALSE,]
# gsea_entriz_cox_df <- gsea_entriz_df[sapply(gsea_entriz_df$Row.names,function(x) x %in% colnames(train_new)[4:25]),]
GSEA_input <- gsea_entriz_df$log2FoldChange
names(GSEA_input) <- gsea_entriz_df$entriz
GSEA_input <- sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = 'hsa', pvalueCutoff = 0.05)#GSEA富集分析
gsea_df <- as.data.frame(GSEA_KEGG)
write.csv(gsea_df, "./data/gsea_df.csv")
# 可视化
library(enrichplot)
pdf(file="./picture/gsea_rich.pdf",width=10,height=10)
gseaplot2(GSEA_KEGG,1:5)
dev.off()


# 6、免疫浸润分析
mianyi_data <- train_new[,c(1,4:25)]
row.names(mianyi_data) <- mianyi_data$id
mianyi_data <- mianyi_data[,-1]
mianyi_data <- t(mianyi_data)
mianyi_data <- as.matrix(mianyi_data)
#处理免疫细胞数据表
marker_gene <- read.csv("./data/CellReports.csv", header = F)
marker_gene <- marker_gene[,-2]
library("reshape2")
marker_gene_melt <- melt(marker_gene,id.vars="V1")
marker_gene_melt <- marker_gene_melt[,-2]
marker_gene_melt <- marker_gene_melt[,c(2,1)]
colnames(marker_gene_melt) <- c("Metagene", "Cell type")
marker_gene_melt <- marker_gene_melt[marker_gene_melt$Metagene != "",]
marker_gene_melt <- marker_gene_melt[order(marker_gene_melt$`Cell type`),]
write.csv(marker_gene_melt, "./data/cell_marker.csv")
#处理样本是高风险还是低风险
sample_type <- train_new[,c(1,27)]
#ssgsea分析
library(GSVA)
library("pheatmap")
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(tinyarray)
library(dplyr)
geneset <- split(marker_gene_melt$Metagene, marker_gene_melt$`Cell type`)
re <- gsva(mianyi_data, geneset, method="ssgsea")
write.csv(re, "./data/mianyi_cell.csv")
group <- sample_type$risk
pdf("./picture/mianyi_cell.pdf", width=6, height=6)
draw_boxplot(re,group,color = c("#1d4a9b","#e5171a"))
dev.off()

# 7、差异免疫细胞浸润水平之间的相关性
library(ggplot2)
risk_cell <- as.data.frame(t(re))
risk_cell$riskScore <- train_new$riskScore
write.csv(risk_cell, "./data/risk_cell.csv")
p <- ggplot(data=risk_cell, aes(x=riskScore, y=`Type 2 T helper cell`)) +
  geom_point(size=1) +
  geom_smooth(method=lm, color="red") +
  scale_color_manual('#FF7400') +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_cor(method="pearson", aes(x=riskScore, y=`Type 2 T helper cell`), label.x=4, label.y=0.6)
ggsave("./picture/mianyin_cell3.pdf",width=6, height=6)

# 8、KM曲线，免疫细胞的生存曲线-------
library("survival")
library("survminer")
mianyi_cell <- read.csv("./data/mianyi_cell.csv", row.names = 1)
mianyi_cell <- t(mianyi_cell)
row.names(mianyi_cell) <- gsub("\\.","-",row.names(mianyi_cell))
new_km <- cbind(train, mianyi_cell)
new_km <- new_km[,c(1,2,25:29)]
write.csv(new_km, "./data/cell_km.csv")
a <- ifelse(new_km$`Type 2 T helper cell` <= median(new_km$`Type 2 T helper cell`), "Low", "High")
fit <- survfit(Surv(surtime, surstat) ~ a, data = new_km)
#按分层更改图形颜色，线型等
pdf('./picture/mianyi_KM3.pdf',w=7,h=6,onefile = F)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,legend.labs=c("High","Low" ),legend.title="Risk", title="Type 2 T helper cell",
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           palette = "npg",
           font.x=15,          
           font.y=15,
           font.title=15,
           font.subtitle = 14,
           font.caption = 14,
           font.tickslab = 14,
           font.legend=17,     
           fontsize = 5,      
)
dev.off()






######################七、构建列线图评估患者1-5年的生存率###############
# 3、列线图预测一年三年五年的生存率
library(survival)
library(rms)
train <- read.csv("./data/train.csv", row.names = 1)
dd<-datadist(train)
options(datadist="dd")
f <- cph(Surv(surtime, surstat) ~ riskScore, data=train,
         x=T, y=T, surv=T)
survival <- Survival(f)
survival1 <- function(x) survival(365, x)
survival2 <- function(x) survival(1095, x)
survival3 <- function(x) survival(1825, x)
nom <- nomogram(f, fun = list(survival1, survival2, survival3),
                fun.at = c(0.05, seq(0.1, 0.9, by=0.05), 0.95),
                funlabel = c('1 year survival', "3 year survival", "5 year survival"))
pdf("./picture/liexiantu1.pdf",height=5,width=8)
plot(nom, cex.axis = 0.5)
dev.off()

# 3.1、校准曲线
f2 <- psm(Surv(surtime, surstat) ~ riskScore, data=train,
          x=T,y=T, dist="lognormal")
call <- calibrate(f2,
                  cmethod='KM',
                  method="boot",
                  u=1825,
                  m=76,
                  B=1000)
pdf("./picture/jiaozhun3.pdf", width=5, height=5)
plot(call, lwd=2, lty=1,
     conf.int=T,
     errbar.col="blue",
     col="red",
     xlab="Nomogram-Predicted Probability of 5-Year OS",
     ylab="Actual 5-Year OS (proportion)",
     subtitles=F)
dev.off()
