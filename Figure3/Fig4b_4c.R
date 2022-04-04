##deep sequence comparison of lung data
##Figure 3b and 3c

library(scales)
library(rDBEC)
library(ggplot2)
library(Seurat)
library(rstatix)
library(dplyr)
library(cowplot)

system("wget -L -O matrix_inflection_10Xv2-LungP7-8.txt.gz 
       https://tus.box.com/shared/static/61wvgbvgirdz3ej3vy7pm3dmnh5xmwi2.gz")
system("wget -L -O matrix_inflection_10Xv2-LungP7-9.txt.gz
https://tus.box.com/shared/static/xlwflls1qhup7hzis74kvp2oz6nhpeyq.gz")

system("wget -L -O matrix_inflection_lung_SmartSeq2.txt.gz 
       https://tus.box.com/shared/static/bsn7pbqy93oh9t2xgx6b7k8y7ibbdz6y.gz")

system("wget -L -O matrix_inflection_TASSeq.deep-1.txt.gz 
       https://tus.box.com/shared/static/gft7tuo3nlyf01ce41uj5kjokkrmig9v.gz")
system("wget -L -O matrix_inflection_TASSeq.deep-2.txt.gz
https://tus.box.com/shared/static/gcqjlcksnr0bdppvvoa2ou5qars971sx.gz")
system("wget -L -O matrix_inflection_TASSeq.deep-3.txt.gz
https://tus.box.com/shared/static/m3n6lnouqjw7qcfjwv0aelpx5x2qihea.gz")

fnames = dir(pattern = "matrix_inflection")
tablelist = lapply(fnames, tableread_fast_sparse)

fnames = c("Smart-seq2", "10Xv2-P7-8", "10Xv2-P7-9", "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3")
names(tablelist)=fnames

for(i in c(1:length(fnames))){
  colnames(tablelist[[i]])=paste(fnames[i], colnames(tablelist[[i]]), sep="_")
}

#remove ERCC spike-in
ERCC = grep(pattern = "^ERCC-", x = rownames(tablelist[[1]]), value = TRUE)
tablelist[[1]] = tablelist[[1]][!rownames() %in% ERCC,]

hoge = as.matrix(tablelist[[1]])
hoge = as.data.frame(hoge)
fwrite(hoge, file = "matrix_inflection_lung_SmartSeq2_noERCC.txt.gz", quote=F, sep="\t", row.names=T, col.names=T)


tablelist = lapply(tablelist, CreateSeuratObject, min.cells = 5, min.features = 200)

mBC = tablelist[[1]]
  for (i in c(2:length(tablelist))){
    mBC = merge(mBC, tablelist[[i]])
  }

#calculate percent mito
mito.genes = grep(pattern = "^mt.", x = rownames(x = mBC@assays$RNA@counts), value = TRUE)
percent.mito = Matrix::colSums(mBC@assays$RNA@counts[mito.genes, ])/Matrix::colSums(mBC@assays$RNA@counts)
names(percent.mito)=rownames(mBC@meta.data)
mBC = AddMetaData(object = mBC, metadata = percent.mito, col.name = "percent.mito")

nReads_log = log10(Matrix::colSums(mBC@assays$RNA@counts))
mBC = AddMetaData(object = mBC, metadata = nReads_log, col.name = "nReads.log")

#filter out mitochondrial high cells
mBC = subset(mBC, subset = percent.mito < 0.25) 

##draw nReads/cell distribution
p = VlnPlot(mBC, features="nReads.log", pt.size=0)

p_legend = get_legend(p)

p = VlnPlot(mBC, features="nReads.log", pt.size=0) + 
    geom_boxplot(outlier.shape=NA, width=0.15) + 
    theme_linedraw() + ylab(expression(paste("Read number (", {log[10]}, ")", sep=""))) +
    ggtitle("")+xlab("")+
    theme(plot.title=element_text(hjust = 0.5), 
        legend.position = 'none', text=element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))

ggsave(plot_grid(p, p_legend, rel_widths = c(0.5,0.4)), filename="nReads_log_deep.png", units="in", dpi=300, 
       width=5, height=3.5, limitsize=FALSE)

##draw nGene/cell distribution
p = VlnPlot(mBC, features="nFeature_RNA", pt.size=0)

p_legend = get_legend(p)

p = VlnPlot(mBC, features="nFeature_RNA", pt.size=0) + 
  geom_boxplot(outlier.shape=NA, width=0.15) + 
  theme_linedraw() + ylab("nGene") +
  ggtitle("")+xlab("")+
  theme(plot.title=element_text(hjust = 0.5), 
        legend.position = 'none', text=element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))

ggsave(plot_grid(p, p_legend, rel_widths = c(0.5,0.4)), filename="nFeature_deep.png", units="in", dpi=300, 
       width=5, height=3.5, limitsize=FALSE)

##reads vs gene plot
hoge = mBC@meta.data
hoge$orig.ident = as.character(hoge$orig.ident)
colnames(hoge)[2:3]=c("nReads", "nGene")

xbreaks = c(0, 1e+03, 1e+04, 1e+05, 1e+06, 3162278)
xbreaks_lab = c(0, 3,4,5,6,6.5)
p = ggplot(hoge, aes(nReads, nGene, group=orig.ident, color=orig.ident))+
  geom_point(data=hoge, aes(nReads, nGene, group=orig.ident), size=0.3, alpha=0.5, shape=20) +
  geom_smooth(method="loess", aes(nReads, nGene), se = FALSE) +
  scale_color_manual(values=c(hue_pal()(6))) +
  scale_fill_manual(values=c(hue_pal()(6))) +
  scale_x_log10(breaks=xbreaks, labels=xbreaks_lab) +
  theme_linedraw() +
  ggtitle("")+
  xlab(expression(paste("Read number (", {log[10]}, ")", sep=""))) + 
  theme(plot.title=element_text(hjust = 0.5), 
        text=element_text(size=12)) +
  guides(color = guide_legend(override.aes = list(size=4, alpha=1)))


ggsave(file="genes_counts.png", plot=p, device="png", units="in",
       dpi=300, width=5, height=4, limitsize = FALSE)

#wilcox test of gene number
res = hoge %>% wilcox_test(nGene ~ orig.ident, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "orig.ident")
colnames(res)[6]="W_statistic"
res=as.data.frame(res)
res=res[,c(2:8), drop=F]
write.table(res, file="stats_nGene.txt", quote=F, sep="\t", row.names=F, col.names=T)

#wilcox test of read number
res = hoge %>% wilcox_test(nReads.log ~ orig.ident, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "orig.ident")
colnames(res)[6]="W_statistic"
res=as.data.frame(res)
res=res[,c(2:8), drop=F]
write.table(res, file="stats_nReads.log.txt", quote=F, sep="\t", row.names=F, col.names=T)

#calculate the number of highly-variable genes
for (i in c(1,2,4:6)){
  mito.genes = grep(pattern = "^mt.", x = rownames(x = tablelist[[i]]@assays$RNA@counts), value = TRUE)
  percent.mito = Matrix::colSums(tablelist[[i]]@assays$RNA@counts[mito.genes, ])/Matrix::colSums(tablelist[[i]]@assays$RNA@counts)
  names(percent.mito)=rownames(tablelist[[i]]@meta.data)
  tablelist[[i]] = AddMetaData(object = tablelist[[i]], metadata = percent.mito, col.name = "percent.mito")
  tablelist[[i]]=subset(tablelist[[i]], subset = percent.mito < 0.25)
}

tablelist = lapply(tablelist, NormalizeData, scale.factor=1000000,verbose = FALSE)
tablelist = lapply(tablelist, FindVariableFeatures, selection.method = "mvp",
                           mean.cutoff=c(0.1, Inf), dispersion.cutoff=c(0.5, Inf))

tmp = NULL
for (i in c(1:6)){
  tmp = c(tmp, length(tablelist[[i]]@assays$RNA@var.features))
}

names(tmp)=fnames
tmp = as.data.frame(tmp)
tmp$orig.ident = fnames
colnames(tmp)[1]="var.genes"

p = ggplot(tmp, aes(y=var.genes, x=orig.ident, color=orig.ident, fill=orig.ident))+
  geom_bar(stat="identity", width=0.8) +
  scale_color_manual(values=c(hue_pal()(6))) +
  scale_fill_manual(values=c(hue_pal()(6))) +
  theme_linedraw() +
  ggtitle("")+
  ylab("Number of highly-variable genes") + 
  xlab("") + 
  theme(plot.title=element_text(hjust = 0.5), legend.position = "none",
        text=element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))

ggsave(file="hvg_number.png", plot=p, device="png", units="in",
       dpi=300, width=2.5, height=4, limitsize = FALSE)

#calculate mean gene number
tmp1 = NULL
for (i in c(1:6)){
  tmp1 = c(tmp1, mean(hoge[hoge$orig.ident==fnames[i],3]))
}
names(tmp1)=fnames
tmp1 = as.data.frame(tmp1)
colnames(tmp1)[1]="mean_gene_number"

tmp2 = NULL
for (i in c(1:6)){
  tmp2 = c(tmp2, mean(hoge[hoge$orig.ident==fnames[i],2]))
}
tmp1$mean_Read_number = tmp2

tmp2 = NULL
for (i in c(1:6)){
  tmp2 = c(tmp2, median(hoge[hoge$orig.ident==fnames[i],3]))
}
tmp1$median_Gene_number = tmp2

write.table(tmp1, file="stats_genes.txt", quote=F, sep="\t", row.names=F, col.names=T)


