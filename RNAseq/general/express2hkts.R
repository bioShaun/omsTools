# HK TS gene selected and visualization
# zxchen 2016-07-11 

library("entropy")
library("ggplot2")
library("reshape2")

geneExpress = read.table("genes.TMM.EXPR.matrix",header = T,sep = "\t")
geneCounts = read.table("genes.counts.matrix",header = T,sep = "\t")
sampleInf =  read.table("sample_inf",header = T,sep = "\t")

numberGene = dim(geneExpressLog2)[1]
numberSample = dim(geneExpressLog2)[2]


# calculate folder change for each gene
cal_folder <- function(gene_array)
{
  S = length(gene_array)
  folder_array = c()
  gene_array = gene_array+0.1
  for (i in seq(1,S-1))
  {
    for (j in seq(i+1,S))
    {
      folder_array = c(folder_array,max(gene_array[c(i,j)])/min(gene_array[c(i,j)]))
      return(max(folder_array))
    }
  }
}



hclust = hclust(dist(express_data))
Cluster = cutree(hclust, k = 13)
express_data = as.data.frame(express_data)
express_data$cluster = Cluster
ggData = melt(express_data,id="cluster")
ggplot(ggData, aes(x=variable, y=value)) + geom_point(size=1)  + 
  geom_boxplot()  + facet_wrap(~cluster, scales = "free",ncol = 2)
  




# Normalization  read conut
geneCounts = log2(1+geneCounts/rep(apply(geneCounts,2,sum),dim(geneCounts)[1],1)*100000000)

#------------------------
# Modify sample name and transfer to log2
geneExpress = geneExpress[,as.character(sampleInf$OrignalID)]
colnames(geneExpress) <- as.character(sampleInf$Sample)
geneExpressLog2 = log2(geneExpress+1)


cutoff = 1.5
TissueGroup = unique(sampleInf$Tissue)
RepeatGroup = unique(sampleInf$Repeat)
folder_list = c()
zscore_list = c()
shanno_list = c()
Tissue_expr = matrix(0,numberGene,length(TissueGroup))
colnames(Tissue_expr) <- TissueGroup
rownames(Tissue_expr) <- rownames(geneExpressLog2)

for (ts in TissueGroup)
{
  Tissue = subset(sampleInf,Tissue==ts,Sample)
  TissueIndex = as.numeric(rownames(Tissue))
  
  group_a = geneExpressLog2[,TissueIndex]
  group_b = geneExpressLog2[,-TissueIndex]
  
  # Tissue group mean value
  Tissue_expr[,ts] = apply(group_a,1,mean)
  
  # TS selected by expression cutoff
  folder = apply(group_a+0.1,1,mean)/apply(group_b+0.1,1,max)
  folder_list = c(folder_list,rownames(geneExpressLog2[folder>cutoff,]))
  
  # TS selected by Zscore
  Zscore = scale(geneExpressLog2,scale = T,center = F)
  Zscore_a = Zscore[,TissueIndex]
  Zscore_b = Zscore[,-TissueIndex]
  Zscore_s = (apply(Zscore_a,1,min)>3) & (apply(Zscore_b,1,max)<2)
  zscore_list = c(zscore_list,rownames(geneExpressLog2[Zscore_s,]))
}

# Golble variant 
expr_mean = apply(geneExpressLog2,2,mean)
expr_avr = apply(geneExpressLog2,1,mean)
expr_max = apply(geneExpressLog2,1,max)
expr_min = apply(geneExpressLog2,1,min)

# HK/TS by Shanno Entropy
matix_expr = as.matrix(geneExpress)
ShannoEntropy = apply(matix_expr,1,entropy)
SH_ts_expr = Tissue_expr[(expr_max > 5 & ShannoEntropy <4 ),]
SH_hk_expr = Tissue_expr[(expr_avr > 3 & ShannoEntropy >5 & expr_min > 1),]

# HK by folder
expr_std = apply(geneExpressLog2,1,sd)
expr_folder = apply(geneExpressLog2,1,cal_folder)
FD_hk_expr = Tissue_expr[(expr_folder<1.2 & expr_avr > 3 & expr_min > 1),]

FD_ts_expr = geneExpressLog2[folder_list,]
ZS_ts_expr = geneExpressLog2[zscore_list,]

folder_mean = Tissue_expr[folder_list,]
zscore_mean = Tissue_expr[zscore_list,]

Hcluster = hclust(dist(t(Tissue_expr))) 

save(Tissue_expr,express_data,ggData,Hcluster,sampleInf,zscore_data,folder_data,zscore_mean,folder_mean,file = "heatmap.RDATA")

sample_cluster =  hclust(dist(t(geneExpressLog2)))



figure_data = t(selected)
sample_id = data.frame(Sample = seq(1,203))
rownames(sample_id) = sample_id$Sample
figure_data = cbind(sample_id,figure_data)
figData = melt(figure_data,id = "Sample")


hc <- hclust(dist(scale(SH_ts_expr,center = T, scale = F),method = "euclidean"),method = "complete")
cluster_id = cutree(hc,)
pheatmap(SH_ts_expr,scale = "row",show_rownames = F,cluster_rows = T,cluster_cols = F,cellwidth = 15)








