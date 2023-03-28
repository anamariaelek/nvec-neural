setwd("~/Dropbox/Cris_RNAseq/")

library('RColorBrewer')
library('scales')
library('viridis')
source("./utility_fxns.r")
library(xlsx)
library("statmod")
library('DESeq2')


Data=read.table("Data_table.txt",header=T,row.names=1,sep="\t")

Data2=Data

############################PLOT EXPRESSION DISTIRBUTIONS########################################
x=Data2
x=x[which(rowSums(x)>0),]
col_dist=brewer.pal(6,"Dark2")
plot(density(log10(x[,1]+1)),col=col_dist[1],lwd=2)
for (i in 2:ncol(x)){
  plot(density(log10(x[,i]+1)),col=col_dist[i],lwd=2)
}
legend('topright',legend=colnames(x),col=col_dist,lwd=2)

png("Exp_distribution.png",h=2500,w=4000)
par(mfrow=c(7,6))
for (i in 1:ncol(x)){
  plot(density(log10(x[,i]+1)),col="black",lwd=6,main=colnames(x)[i],cex.main=5)
}
dev.off()

barplot(apply(x,2,function(x) sum(x>10)),las=2)


#########################GLMFIT to select to variable genes####################################
ff=Data2
cv2<-apply(ff[,],1,function(x) var(x)/(mean(x)^2))
means<-apply(ff,1,mean)
vars<-apply(ff,1,var)


smoothScatter(log(cv2)~log(means))
#minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
useForFit <- means >= 50 # & spikeins
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients

# repeat previous plot
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2));
xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
vfit <- a1/xg + a0
# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(x) - 1
# add confidence interval
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")

afit <- a1/means+a0
varFitRatio <- vars/(afit*means^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- ff[varorder,]

# repeat previous plot
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# add top 100 genes
points(log(means[varorder[1:500]]),log(cv2[varorder[1:500]]),col=2)




############################PLOT ALL SAMPLE-SAMPLE  CORR########################################
Data_caca=cbind(Data2[,which(grepl("Treg",colnames(Data2),perl = T))])#,
                #Data2[,which(grepl("Treg_d7",colnames(Data2),perl=T))])
                #Data2[,which(grepl("early_d",colnames(Data2),perl=T))],
                #Data2[,which(grepl("Tconv_d",colnames(Data2),perl=T))])
Data_caca=Data2


Data_transf=as.matrix(Data_caca[which(rowSums(Data_caca)>0 & rowMeans(Data_caca)<1000000),])
#Data_transf=quantile_normalization(Data_transf)
Data_transf_norm=(Data_transf[,]+10)/apply(Data_transf[,]+10,1,median)

vm_cor_analysis=apply(Data_transf[,], 1, function(x) var(x)/mean(x))
vm_cor_analysis[is.na(vm_cor_analysis)]=0
top_varmean_genes=names(sort(vm_cor_analysis,decreasing=TRUE))[1:7000]


##Option 1: Gene markers by high FC
genes_high_FC=names(which(apply(Data_transf_norm,1,function(x) any(x>1))))
genes_high_FC=names(which(apply(Data_transf_norm, 1, function(x) sort(x,decreasing = T)[2])>1.8))
gene_markers=genes_high_FC

##Option 2: Gene markers by varmean order
gene_markers=names(means[varorder[1:500]])


##Option 3: Gene markers by high FC+significant kruskall wallis
factor=c(rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2),rep(6,2),rep(7,2),rep(8,2),rep(9,2))    ########CHANGE according to your experimental design

kw_pvals=apply(Data_transf,1,function(x) kruskal.test(x~factor)$p.value)
kw_qvals=p.adjust(kw_pvals,method="BH")
kw_genes=names(which(kw_qvals<0.01))
genes_high_FC=names(which(apply(Data_transf_norm, 1, function(x) sort(x,decreasing = T)[2])>2))
gene_markers=intersect(genes_high_FC,kw_genes)


#Option 4: Gene markers by high FC+significant DEseq2 LTR test (see below)
deseq_qval=deseq_res$padj
names(deseq_qval)=rownames(deseq_res)
deseq_genes=names(which(deseq_qval<1e-2))
#genes_high_FC=names(which(apply(Data_transf_norm, 1, function(x) sort(x,decreasing = T)[4])>2))
#gene_markers=intersect(genes_high_FC,deseq_genes)
gene_markers=deseq_genes


##Option 0. CAREFUL. This is only to enforce specific pairwise comparisons. USE WITH CARE AND CHANGE x indices accordingly
pvals=apply(Data_caca,1,function(x) wilcox.test(x[1:16],x[17:32])$p.value)
qvals=p.adjust(pvals,method="BH")
gene_markers=names(which(qvals < 0.02))
######

cor_data=cor(Data_transf[gene_markers,],method="pearson")
diag(cor_data)=NA
hclust_data=hclust(as.dist(1-cor_data),method="ward.D2")
nms=rownames(cor_data[hclust_data$order,hclust_data$order])
plot(hclust_data)

png("Clustering_corr_Treg_samples.png",h=3000,w=3000)
par(fig=c(0.215,0.795,0.7,0.98))
plot(hclust_data,xaxt='n',yaxt='n', xaxs="i",sub="",ylab="",cex=1,main="",xlab="",
     lwd = 5,hang=-1,labels=FALSE)

par(fig=c(0.2,0.8,0.2,0.75),new=TRUE)
cor_shades = colorRampPalette(c("white","white","darkgoldenrod1","brown3","black"))(1000)
#cor_shades=rev(viridis(1000,option="magma"))
diag(cor_data)=1
image(cor_data[nms,nms],col=cor_shades, xaxt='n',yaxt='n',zlim=c(min(cor_data,na.rm=TRUE),max(cor_data,na.rm=TRUE)))
mtext(nms, side=2,at=seq(0,1,length.out=length(nms)),las=1,adj=1,cex=4)
mtext(nms, side=1,at=seq(0,1,length.out=length(nms)), las=2,cex=4)
dev.off()

###################Now show gene expression bietch

#Option 1. Cluster genes Hclust
cor_genes=cor(t(Data_transf_norm[gene_markers,]))  ##############WHICH GENES YOU WANT TO LOOK AT BITCH???????
#diag(cor_genes)=NA
hclust_genes=hclust(as.dist(1-cor(cor_genes)),method="ward.D2")
gene_order=rownames(cor_genes[hclust_genes$order,hclust_genes$order])

#Option 2. Kmeans genes clsut
K=15
kmeans_genes=TGLKMeans_wrapper(data = log2(Data_transf_norm[gene_markers,]+1),fn="km.markers",k=K)
gene_module_order=as.character(hclust(dist(cor(cor(t(kmeans_genes$centers)))),method="ward.D2")$order)

gene_order=c()
for(module in gene_module_order){ gene_order=c(gene_order,names(kmeans_genes$cluster[which(kmeans_genes$cluster==module)]))   }

OUTPUT="Treg_d2_vs_Treg_d7_DESeqe-5"
plot_vertical=F

write.xlsx(rev(kmeans_genes$cluster[gene_order]),file=paste0(OUTPUT,".xlsx"))

shades2=colorRampPalette(brewer.pal(11,'BrBG'))(1000)
###CHECKOUT FTP NORMALIZATION
#tf_fp_to_plot=pmin((0.1+tf_fp[nms,niche_order])/(0.1+apply(tf_fp[nms,],1,median)),8)
matrix_to_plot=log2(Data_transf_norm[gene_order,])
#matrix_to_plot=pmin(pmax(matrix_to_plot-rowMeans(matrix_to_plot),-3),3)
matrix_to_plot=pmin(pmax(matrix_to_plot,-min(abs(range(matrix_to_plot))),-min(abs(range(matrix_to_plot)))))
png(paste0(OUTPUT,".png"),h=4000,w=800)
par(mar=c(1,15,10,1))
image(t(matrix_to_plot[,]), xaxt='n',yaxt='n', col=shades2)

if(plot_vertical==T){
  factor_sizes=as.numeric(table(factor))
  current.line = 0
  for (i in 1:length(factor_sizes)){
   current.line = current.line + factor_sizes[i]
   abline(v = current.line/(sum(factor_sizes)-1) - 1/(2*sum(factor_sizes)), col = "gray50", lwd = 4, lty = 3)
  }
}
module_sizes=as.numeric(table(kmeans_genes$cluster[gene_order])[as.character(unique(kmeans_genes$cluster[gene_order]))])
current.line = 0
for (i in 1:length(module_sizes)){
  current.line = current.line + module_sizes[i]
  abline(h = current.line/(sum(module_sizes)-1) - 1/(2*sum(module_sizes)), col = "gray50", lwd = 2, lty = 3)
  mtext(as.character(unique(kmeans_genes$cluster[gene_order]))[i], side=4, at=current.line/(sum(module_sizes)-1) - 1/(2*sum(module_sizes)), las=2, line=0,cex=1)
}
mtext(colnames(matrix_to_plot), side=3,at=seq(0,1,length.out=length(colnames(matrix_to_plot))), las=2,cex=1)
mtext(gene_order, side=2,at=seq(0,1,length.out=length(gene_order)), las=2,cex=1)
dev.off()





image(x=seq(-max(matrix_to_plot),max(matrix_to_plot),(2*max(matrix_to_plot))/(length(shades2)-1)), 
      y=c(0,1), col=shades2, yaxt="n", z=matrix(nrow=length(shades2),ncol=1,data=c(1:length(shades2))),
      ylab="",xlab="",cex.axis=2)
cor_shades = colorRampPalette(c("darkblue","blue","white", "gold","brown"))(1000)
image(cor_genes[gene_order,gene_order],col=cor_shades, xaxt='n',yaxt='n',zlim=c(-1,1))
mtext("log2 FC",side=1,line=5,cex=3)



#----------------Quantile norm function-----------------------------------------------------------------------
quantile_normalization <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, function(x) { sort(x, na.last=TRUE)}))
  df_mean <- apply(df_sorted, 1, mean, na.rm=TRUE)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}






####################################################################################################################
####################################################################################################################


deseq_factors=read.table("Data_factors_DESeq.txt")
cts=read.table("Data_Table_RAW_COUNTS.txt")

counts_sub=cts

counts_sub=cts[,colnames(Data_caca)]
deseq_factors_sub=deseq_factors[colnames(counts_sub),]

dds_treated <- DESeqDataSetFromMatrix(countData = counts_sub,colData = deseq_factors_sub,design = ~ treated)
dds_celltype <- DESeqDataSetFromMatrix(countData = counts_sub,colData = deseq_factors_sub,design = ~ cell_type)
dds_time <- DESeqDataSetFromMatrix(countData = counts_sub,colData = deseq_factors_sub,design = ~ time)
dds_cd73 <- DESeqDataSetFromMatrix(countData = counts_sub,colData = deseq_factors_sub,design = ~ CD73)
dds_multifactor <- DESeqDataSetFromMatrix(countData = counts_sub,colData = deseq_factors_sub,design = ~ multifactor)

#Two-sided test
dds <- DESeq(dds_cd73)
deseq_res <- results(dds)



#Factor with multiple levels (ANOVA-like):
dds <- DESeq(dds_CD73, test="LRT", reduced=~1)
deseq_res <- results(dds)



######################################Barplot gene lists############################################################
####################################################################################################################
glist=scan("Specific_barplots/hypofree_highinLate.txt",what="")

spaces=c(0.1,0.1,0.1,0.1,0.5,0.1,0.1,0.1,0.5,0.1,0.5,0.1,0.1,0.1,0.5,0.1,0.1,0.1,0.5,0.1,0.1,
  0.1,0.5,0.1,0.1,0.1,0.5,0.1,0.1,0.1,0.5,0.1,0.1,0.1,0.51,0.1,0.1,0.1,0.5,0.1,0.1,0.1)


for(gene in glist){
  png(paste0("./Specific_barplots/",gene,".png"),h=600,w=1000)
  par(mar=c(12,7,6,6))
  barplot(as.matrix(Data[gene,]),las=2,col="gray33",main=gene,cex.main=3,cex.axis=1.5,space = spaces,border="gray33")
  dev.off()
}




