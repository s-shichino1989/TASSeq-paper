##deep sequence comparison of spleen data
##Supplementary Figure 3d
library(scales)
library(rDBEC)
library(ggplot2)
library(Seurat)
library(rstatix)
library(dplyr)
library(cowplot)

##download annotated Seurat objects for mouse lung data
system("wget -L -O matrix_inflection_MouseSpleen1_deep.txt.gz https://tus.box.com/shared/static/d4b73innhjxa78tp6xwlgusehtuq3upi.gz")
system("wget -L -O matrix_inflection_MouseSpleen2_deep.txt.gz https://tus.box.com/shared/static/aah2fg88xex5y2sfmxtfwys8x9hqlhq0.gz")
system("wget -L -O matrix_inflection_MouseSpleen3_deep.txt.gz https://tus.box.com/shared/static/vltu1q6356ubwq27k64k31xczt9lblzf.gz")
system("wget -L -O matrix_inflection_SpleenSmartSeq2.txt.gz https://tus.box.com/shared/statzic/i25lv5ujxk3hmp38gvdvuj3cnwe30drw.gz")


fnames = dir(pattern = "matrix_")
tablelist = lapply(fnames, tableread_fast_sparse)

fnames = c("TAS-Seq.deep-1","TAS-Seq.deep-2","TAS-Seq.deep-3", "Smart-seq2")
names(tablelist)=fnames

for(i in c(1:length(fnames))){
  colnames(tablelist[[i]])=paste(fnames[i], colnames(tablelist[[i]]), sep="_")
}

tablelist = lapply(tablelist, CreateSeuratObject, min.cells = 5, min.features = 200)

#####

seu1 = tablelist[[1]]
seu2 = tablelist[[2]]
seu3 = tablelist[[3]]
seu = merge(seu1, seu2)
seu = merge(seu, seu3)
seu1 = seu
raw.data.seurat = seu1@assays$RNA@counts

mito.genes = grep(pattern = "^mt.", x = rownames(x = raw.data.seurat), value = TRUE)
percent.mito = Matrix::colSums(raw.data.seurat[mito.genes, ])/Matrix::colSums(raw.data.seurat)
names(percent.mito)=rownames(seu1@meta.data)
seu1 = AddMetaData(object = seu1, metadata = percent.mito, col.name = "percent.mito")

fuga = data.frame(nGene = seu1@meta.data$nFeature_RNA, mito.proportion = seu1@meta.data$percent.mito)
p.mito = ggplot(fuga) +
  geom_hex(aes(nGene, mito.proportion),
           bins = 100,
           show.legend = TRUE) +
  ylim(0, 1) +
  xlim(0, 14000) +
  scale_fill_gradientn("cell number", colours = rev(rainbow(10, end = 4/6))) +
  geom_hline(aes(yintercept=0.25), colour="magenta", size=0.5) +
  ggtitle("Mitochondrial gene proportion") +
  theme_linedraw() +
  theme(plot.title=element_text(hjust = 0.5), text=element_text(size=10)) +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10))

ggsave(p.mito, filename="mito_TASSeq_SpleenTAS-Seq.png", units="in", dpi=300,
       width=5, height=3, limitsize=FALSE)

seu1 = tablelist[[4]]
raw.data.seurat = seu1@assays$RNA@counts

mito.genes = grep(pattern = "^mt.", x = rownames(x = raw.data.seurat), value = TRUE)
percent.mito = Matrix::colSums(raw.data.seurat[mito.genes, ])/Matrix::colSums(raw.data.seurat)
names(percent.mito)=rownames(seu1@meta.data)
seu1 = AddMetaData(object = seu1, metadata = percent.mito, col.name = "percent.mito")

fuga = data.frame(nGene = seu1@meta.data$nFeature_RNA, mito.proportion = seu1@meta.data$percent.mito)
p.mito = ggplot(fuga) +
  geom_hex(aes(nGene, mito.proportion),
           bins = 100,
           show.legend = TRUE) +
  ylim(0, 1) +
  xlim(0, 14000) +
  scale_fill_gradientn("cell number", colours = rev(rainbow(10, end = 4/6))) +
  geom_hline(aes(yintercept=0.25), colour="magenta", size=0.5) +
  ggtitle("Mitochondrial gene proportion") +
  theme_linedraw() +
  theme(plot.title=element_text(hjust = 0.5), text=element_text(size=10)) +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10))

ggsave(p.mito, filename="mito_smartseq2_Spleen.png", units="in", dpi=300,
       width=5, height=3, limitsize=FALSE)

##merge data
mBC = tablelist[[1]]
for (i in c(2:length(tablelist))){
  mBC = merge(mBC, tablelist[[i]])
}

colnames(mBC@meta.data)[3]="nGene"
colnames(mBC@meta.data)[2]="nReads"

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

ggsave(plot_grid(p, p_legend, rel_widths = c(0.5,0.4)), filename="nReads_log_spleen.png", units="in", dpi=300,
       width=4, height=3.5, limitsize=FALSE)

##draw nGene/cell distribution
p = VlnPlot(mBC, features="nFeature_RNA", pt.size=0)

p_legend = get_legend(p)

p = VlnPlot(mBC, features="nFeature_RNA", pt.size=0) +
  geom_boxplot(outlier.shape=NA, width=0.15) +
  theme_linedraw() + ylab("nGene") +
  ggtitle("")+xlab("")+
  theme(plot.title=element_text(hjust = 0.5),
        legend.position = 'none', text=element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))

ggsave(plot_grid(p, p_legend, rel_widths = c(0.5,0.4)), filename="nFeature_spleen.png", units="in", dpi=300,
       width=4, height=3.5, limitsize=FALSE)

##reads vs gene plot
hoge = mBC@meta.data
hoge$orig.ident = as.character(hoge$orig.ident)
colnames(hoge)[2:3]=c("nReads", "nGene")

xbreaks = c(0, 1e+03, 1e+04, 1e+05, 1e+06)
xbreaks_lab = c(0, 3,4,5,6)
p = ggplot(hoge, aes(nReads, nGene, group=orig.ident, color=orig.ident))+
  geom_point(data=hoge, aes(nReads, nGene, group=orig.ident), size=0.3, alpha=0.5, shape=20) +
  geom_smooth(method="loess", aes(nReads, nGene), se = FALSE) +
  scale_color_manual(values=c(hue_pal()(4))) +
  scale_fill_manual(values=c(hue_pal()(4))) +
  scale_x_log10(breaks=xbreaks, labels=xbreaks_lab) +
  theme_linedraw() +
  ggtitle("")+
  xlab(expression(paste("Read number (", {log[10]}, ")", sep=""))) +
  theme(plot.title=element_text(hjust = 0.5),
        text=element_text(size=12)) +
  guides(color = guide_legend(override.aes = list(size=4, alpha=1)))

ggsave(file="genes_counts_speen.png", plot=p, device="png", units="in",
       dpi=300, width=5, height=4, limitsize = FALSE)

#wilcox test of gene number
res = hoge %>% wilcox_test(nGene ~ orig.ident, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "orig.ident")
colnames(res)[6]="W_statistic"
res=as.data.frame(res)
res=res[,c(2:8), drop=F]
write.table(res, file="stats_nGene_spleen.txt", sep="\t", quote=F, row.names=F, col.names=T)

#wilcox test of read number
res = hoge %>% wilcox_test(nReads.log ~ orig.ident, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "orig.ident")
colnames(res)[6]="W_statistic"
res=as.data.frame(res)
res=res[,c(2:8), drop=F]
write.table(res, file="stats_nReads_spleen.log.txt", quote=F,  sep="\t", row.names=F, col.names=T)

#calculate the number of highly-variable genes
for (i in c(1:4)){
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
for (i in c(1:4)){
  tmp = c(tmp, length(tablelist[[i]]@assays$RNA@var.features))
}

names(tmp)=fnames
tmp = as.data.frame(tmp)
tmp$orig.ident = fnames
colnames(tmp)[1]="var.genes"

p = ggplot(tmp, aes(y=var.genes, x=orig.ident, color=orig.ident, fill=orig.ident))+
  geom_bar(stat="identity", width=0.8) +
  scale_color_manual(values=c(hue_pal()(4))) +
  scale_fill_manual(values=c(hue_pal()(4))) +
  theme_linedraw() +
  ggtitle("")+
  ylab("Number of highly-variable genes") +
  xlab("") +
  theme(plot.title=element_text(hjust = 0.5), legend.position = "none",
        text=element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))

ggsave(file="hvg_number_spleen.png", plot=p, device="png", units="in",
       dpi=300, width=2.5, height=4, limitsize = FALSE)

#calculate mean gene number
tmp1 = NULL
for (i in c(1:4)){
  tmp1 = c(tmp1, mean(hoge[hoge$orig.ident==fnames[i],3]))
}
names(tmp1)=fnames
tmp1 = as.data.frame(tmp1)
colnames(tmp1)[1]="mean_gene_number"

tmp2 = NULL
for (i in c(1:4)){
  tmp2 = c(tmp2, mean(hoge[hoge$orig.ident==fnames[i],2]))
}
tmp1$mean_Read_number = tmp2

tmp2 = NULL
for (i in c(1:4)){
  tmp2 = c(tmp2, median(hoge[hoge$orig.ident==fnames[i],3]))
}
tmp1$median_Gene_number = tmp2

write.table(tmp1, file="stats_genes_spleen.txt", quote=F, sep="\t", row.names=F, col.names=T)

