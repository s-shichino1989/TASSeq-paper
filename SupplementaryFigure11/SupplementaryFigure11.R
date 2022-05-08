##genes-reads distribution of pseudo-bulk data 
##Supplementary Figure 11

library(scales)
library(rDBEC)
library(ggplot2)
library(Seurat)
library(rstatix)
library(dplyr)
library(cowplot)
library(ggExtra)
library(tibble)
library(qs)

##download annotated Seurat objects for mouse lung data
system("wget -L -O UTlung_deep_Seurat_annot.qs https://tus.box.com/shared/static/ohoohhgf3t12tcmupueajbbapu4xdxm6.qs")
system("wget -L -O UTlung_shallow_Seurat_annot.qs https://tus.box.com/shared/static/rhlcprxv2l3mtmulphsotwl2serl6u7o.qs")
system("wget -L -O UTlung_SmartSeq2_Seurat_annot.qs https://tus.box.com/shared/static/ej8tff4k7szdfqsw6niljnzrronnx910.qs")
system("wget -L -O UTlung_10Xv2_TabulaMuris_Seurat_annot.qs https://tus.box.com/shared/static/fqiabrfwxd7y15jgcwlovelhoz2tgz65.qs")

##download GSE110540 bulkl RNA-seq (3'SAGE-Seq) data
system("wget -L -O SAGE_data.txt https://tus.box.com/shared/static/zt0iy2eh4x1qxlwtn4gbff5ilsdehmwx.txt")

files = dir(pattern = ".qs")

seu_list = lapply(files, qread, nthreads=24)

#extract commonly-detected cell subset
ident = intersect(unique(as.character(seu_list[[1]]@meta.data$celltype)), unique(as.character(seu_list[[2]]@meta.data$celltype)))
ident = intersect(ident, unique(as.character(seu_list[[3]]@meta.data$celltype)))
ident = intersect(ident, unique(as.character(seu_list[[4]]@meta.data$celltype)))
ident = unique(ident)
ident = sort(ident)

#separate seurat object

#10X v2 Tabula Muris

seu1 = seu_list[[1]]
tmp = seu1@meta.data$orig.ident
names(tmp)=rownames(seu1@meta.data)
tmp = factor(tmp)
seu1@active.ident = tmp
seu2 = subset(seu1, idents = "10Xv2-LungP7-8")
seu3 = subset(seu1, idents = "10Xv2-LungP7-9")

tmp = as.character(seu2@meta.data$celltype)
names(tmp)=rownames(seu2@meta.data)
tmp = factor(tmp)
seu2@active.ident=tmp

tmp = as.character(seu3@meta.data$celltype)
names(tmp)=rownames(seu3@meta.data)
tmp = factor(tmp)
seu3@active.ident=tmp

#TAS-Seq.deep
seu1 = seu_list[[2]]
tmp = seu1@meta.data$orig.ident
names(tmp)=rownames(seu1@meta.data)
tmp = as.factor(tmp)
seu1@active.ident = tmp
seu4 = subset(seu1, idents = "TASSeq.deep-1")
seu5 = subset(seu1, idents = "TASSeq.deep-2")
seu6 = subset(seu1, idents = "TASSeq.deep-3")

tmp = as.character(seu4@meta.data$celltype)
names(tmp)=rownames(seu4@meta.data)
tmp = as.factor(tmp)
seu4@active.ident=tmp

tmp = as.character(seu5@meta.data$celltype)
names(tmp)=rownames(seu5@meta.data)
tmp = as.factor(tmp)
seu5@active.ident=tmp

tmp = as.character(seu6@meta.data$celltype)
names(tmp)=rownames(seu6@meta.data)
tmp = as.factor(tmp)
seu6@active.ident=tmp


#TAS-Seq.shallow
seu1 = seu_list[[3]]
tmp = seu1@meta.data$orig.ident
names(tmp)=rownames(seu1@meta.data)
tmp = as.factor(tmp)
seu1@active.ident = tmp
seu7 = subset(seu1, idents = "TASSeq.shallow-1")
seu8 = subset(seu1, idents = "TASSeq.shallow-2")
seu9 = subset(seu1, idents = "TASSeq.shallow-3")

tmp = as.character(seu7@meta.data$celltype)
names(tmp)=rownames(seu7@meta.data)
tmp = as.factor(tmp)
seu7@active.ident=tmp

tmp = as.character(seu8@meta.data$celltype)
names(tmp)=rownames(seu8@meta.data)
tmp = as.factor(tmp)
seu8@active.ident=tmp

tmp = as.character(seu9@meta.data$celltype)
names(tmp)=rownames(seu9@meta.data)
tmp = as.factor(tmp)
seu9@active.ident=tmp

#create new Seurat list
seu_list1 = list()
seu_list1[[1]] = seu2
seu_list1[[2]] = seu3
seu_list1[[3]] = seu_list[[4]]
seu_list1[[4]] = seu7
seu_list1[[5]] = seu8
seu_list1[[6]] = seu9
seu_list1[[7]] = seu4
seu_list1[[8]] = seu5
seu_list1[[9]] = seu6
seu_list = seu_list1

#create pseudo-bulk count data for each commonly-detected cell subset

res = list()

for (j in c(1:9)){
 for (i in c(1:length(ident))){
  seu_tmp = subset(seu_list[[j]], idents = ident[i])
  seu_tmp = seu_tmp@assays$RNA@counts
  seu_tmp = as.data.frame(rowSums(seu_tmp))
  colnames(seu_tmp)=ident[i]
  if(i > 1){
  res_tmp = dplyr::full_join(rownames_to_column(res_tmp), rownames_to_column(seu_tmp), by="rowname")
  } else {
  res_tmp = rownames_to_column(seu_tmp)
  }
  rownames(res_tmp)=res_tmp$rowname
  res_tmp = res_tmp[,c(2:ncol(res_tmp)), drop=F]
  res_tmp = as.data.frame(res_tmp)
 }
  res_tmp[is.na(res_tmp)]=0
  res[[j]] = res_tmp
}

#normalized count data to total 10M reads
for (i in c(1:9)){
  nf = 10000000 / colSums(res[[i]])
  res[[i]] = sweep(res[[i]], MARGIN=2, nf, FUN = "*")
}

#create cell type-separated matrix list
res_list = NULL
for (i in c(1:length(ident))){
  tmp = dplyr::full_join(rownames_to_column(res[[1]][,i, drop=F]),
                         rownames_to_column(res[[2]][,i,drop=F]),
                         by="rowname")
  for (j in c(3:9)){
  tmp = dplyr::full_join(tmp, rownames_to_column(res[[j]][,i,drop=F]),
                         by="rowname")
  }
  rownames(tmp)=tmp$rowname
  tmp = tmp[,c(2:10)]
  tmp[is.na(tmp)]=0
  colnames(tmp)=c("10Xv2-P7-8", "10Xv2-P7-9", "Smart-seq2",
                  "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3",
                  "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3")
  res_list[[i]]=tmp
}

#create ridgeline plot

library(ggplot2)

plot_res = list()

for(i in c(1:length(ident))){
tmp_df = res_list[[i]]
tmp_df$gene = rownames(tmp_df)
tmp_df = tmp_df %>% tidyr::gather(key = dataset, value = expression, -gene)
tmp_df = tmp_df[tmp_df$expression>0,]

tmp_df$expression = log10(tmp_df$expression)
p = ggplot(tmp_df, aes(x=factor(dataset, levels = c("10Xv2-P7-8", "10Xv2-P7-9", "Smart-seq2",
                                                   "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3",
                                                   "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3")),
                       y = expression, color=dataset, fill=dataset, alpha=0.8)) +
  geom_violin(aes(fill=dataset), color="black", alpha=1, scale="width") +
  geom_boxplot(aes(fill=dataset), width=0.2, color="black", alpha=0.2, outlier.shape=NA) +
  scale_color_manual(values=c(hue_pal()(9))) +
  scale_fill_manual(values=c(hue_pal()(9))) +
  theme_linedraw() +
  theme(legend.position="none") +
  ggtitle(ident[i]) +
  scale_y_continuous(limits=c(0, max(c(tmp_df$expression, 8))), label=label_number_si(accuracy = 0.1))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust=0.5,size=12),
        plot.margin = unit(c(1,2,0,1), "mm")) +
  xlab("") + ylab("")


plot_res[[i]]=p
}


#create figure legend
p = ggplot(tmp_df, aes(x=factor(dataset, levels = c("10Xv2-P7-8", "10Xv2-P7-9", "Smart-seq2",
                                                   "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3",
                                                   "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3")),
                       y = expression, color=dataset, fill=dataset, alpha=0.8)) +
  geom_violin(aes(fill=dataset), color="black", alpha=1, scale="width")

p = get_legend(p)
plot_res[[18]] = p

ggsave(file="gene_reads_distribution_pseudobulk.png", plot=plot_grid(plotlist = plot_res, ncol = 5), device="png", units="in",
       dpi=300, width=15, height=12, limitsize = FALSE, bg="white")


#read SAGE-Seq data

hoge = read.table("./SAGE_test.txt", row.names=1, header=T, sep="\t", quote="")
#extract UT fiboblast data
hoge = hoge[,c(1:3)]

#normalized count data to total 10M reads
nf = 10000000 / colSums(hoge)
hoge = sweep(hoge, MARGIN=2, nf, FUN = "*")

#rename columns
colnames(hoge)=c("bulk-1", "bulk-2", "bulk-3")

#merge expression table
tmp_df = res_list[[9]]
tmp_df$gene = rownames(tmp_df)
tmp_df = tmp_df %>% tidyr::gather(key = dataset, value = expression, -gene)
tmp_df = tmp_df[tmp_df$expression>0,]

tmp_df1 = hoge
tmp_df1$gene = rownames(tmp_df1)
tmp_df1 = tmp_df1 %>% tidyr::gather(key = dataset, value = expression, -gene)
tmp_df1 = tmp_df1[tmp_df1$expression>0,]

tmp_df = rbind(tmp_df, tmp_df1)

tmp_df$expression = log10(tmp_df$expression)

#count proportion of lower quantile
res=NULL
datasets = c("10Xv2-P7-8", "10Xv2-P7-9", "Smart-seq2",
             "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3",
             "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3",
             "bulk-1", "bulk-2", "bulk-3")
for (i in c(1:length(datasets))){
  tmp_df1 = tmp_df[tmp_df$dataset %in% datasets[i],]
  tmp1 = quantile(tmp_df1$expression, c(0.25, 0.75))
  prop = nrow(tmp_df1[tmp_df1$expression <= tmp1[1],]) / nrow(tmp_df1)
  res = rbind(res, tmp_df1)
}

tmp_df = res


p = ggplot(tmp_df, aes(x=factor(dataset, levels = c("10Xv2-P7-8", "10Xv2-P7-9", "Smart-seq2",
                                                    "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3",
                                                    "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3",
                                                    "bulk-1", "bulk-2", "bulk-3")),
                       y = expression, color=dataset, fill=dataset, alpha=0.8)) +
  geom_violin(aes(fill=dataset), color="black", alpha=1, scale="width") +
  geom_boxplot(aes(fill=dataset), width=0.2, color="black", alpha=0.2, outlier.shape=NA) +
  scale_color_manual(values=c(hue_pal()(12))) +
  scale_fill_manual(values=c(hue_pal()(12))) +
  theme_linedraw() +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, max(c(tmp_df$expression, 8))), label=label_number_si(accuracy = 0.1))+
  ggtitle(ident[9]) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=12),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust=0.5,size=12),
        plot.margin = unit(c(1,2,0,1), "mm")) +
  xlab("") +  ylab("")

ggsave(file="gene_reads_distribution_pseudobulk_fibro.png", plot=p, device="png", units="in",
       dpi=300, width=4, height=4, limitsize = FALSE, bg="white")

##extract kernel density information

library(philentropy)
p = ggplot(tmp_df, aes(x=factor(dataset, levels = c("10Xv2-P7-8", "10Xv2-P7-9", "Smart-seq2",
                                                    "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3",
                                                    "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3",
                                                    "bulk-1", "bulk-2", "bulk-3")),
                       y = expression, color=dataset, fill=dataset, alpha=0.8)) +
  geom_violin(aes(fill=dataset), color="black", alpha=1, scale="width")
q = ggplot_build(p)
q = q$data[[1]]

#calculate distance between scaled kernel densities
res = NULL
for(i in c(1:12)){
  tmp = q[q$group == i,"ndensity"]
  tmp = tmp / sum(tmp)
  res = rbind(res, tmp)
}

rownames(res) = c("10Xv2-P7-8", "10Xv2-P7-9","bulk-1", "bulk-2", "bulk-3", "Smart-seq2",
                  "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3",
                  "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3")

#calculate distance by kullback-leibler method
hoge = distance(res, method = "kullback-leibler")

rownames(hoge)=c("10Xv2-P7-8", "10Xv2-P7-9","bulk-1", "bulk-2", "bulk-3", "Smart-seq2",
                 "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3",
                 "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3")
colnames(hoge) = c("10Xv2-P7-8", "10Xv2-P7-9","bulk-1", "bulk-2", "bulk-3", "Smart-seq2",
                   "TAS-Seq.deep-1", "TAS-Seq.deep-2", "TAS-Seq.deep-3",
                   "TAS-Seq.shallow-1", "TAS-Seq.shallow-2", "TAS-Seq.shallow-3")

#export distance matrix
write.table(hoge, file="KL_distance_density.txt", row.names = T, col.names = T, sep="\t", quote=F)

#visualize distance by heatmap

hoge = hoge[,c(1,2,6:12,3,4,5)]
hoge = hoge[c(1,2,6:12,3,4,5),]
library(viridis)
library(pheatmap)
library(ggplotify)

p = pheatmap(
  hoge,
  scale="none",
  show_colnames=TRUE,
  show_rownames=TRUE,
  color=plasma(1024),
  silent=TRUE
)

p = as.ggplot(p[[4]])
ggsave(file = "KL_distance_Heatmap.png", plot = p,
       device="png", units="in", dpi = 300,
       width = 5, height = 4, limitsize=FALSE)



