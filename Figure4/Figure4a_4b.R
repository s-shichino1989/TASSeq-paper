###############
#draw figure 4a and 4b
library(qs)
library(rDBEC)
library(Seurat)
library(pheatmap)
library(dplyr)
library(tibble)
library(Matrix)
library(ggplot2)
library(cowplot)


custom_colors <- list()
colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)
colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)
colors_custom = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
                  '#49beaa', '#611c35', '#2708a0')

custom_colors$discrete <- unique(c(colors_dutch, colors_spanish, colors_custom))


##download annotated Seurat objects for mouse lung data
system("wget -L -O matrix_inflection_10Xv2-GSM3926540.txt.gz 
       https://tus.box.com/shared/static/wf4whgfdj2bx0kpbe05r66ntq1p7zltl.gz")

system("wget -L -O UTlung_deep_Seurat_annot.qs
       https://tus.box.com/shared/static/ohoohhgf3t12tcmupueajbbapu4xdxm6.qs")
system("wget -L -O UTlung_shallow_Seurat_annot.qs
       https://tus.box.com/shared/static/rhlcprxv2l3mtmulphsotwl2serl6u7o.qs")
system("wget -L -O UTlung_SmartSeq2_Seurat_annot.qs
       https://tus.box.com/shared/static/ej8tff4k7szdfqsw6niljnzrronnx910.qs")
system("wget -L -O UTlung_10Xv2_TabulaMuris_Seurat_annot.qs
       https://tus.box.com/shared/static/fqiabrfwxd7y15jgcwlovelhoz2tgz65.qs")



fnames = dir(pattern = ".qs")
seu_list = lapply(fnames, qread, nthreads=24)

#rename Smart-seq2 orig.ident
seu_list[[5]]@meta.data$orig.ident = rep("Smart-seq2", nrow(seu_list[[5]]@meta.data))

##################
#plotting Figure 4a left (FIt-SNE plot)
hoge = sort(unique(as.character(seu@meta.data$celltype)))
hoge1 = custom_colors$discrete[1:length(hoge)]
names(hoge1)=hoge

title = c("10Xv2-GSM3926540", "10Xv2-TabulaMuris", "TAS-Seq.deep",
          "TAS-Seq.shallow", "Smart-seq2")

i=1
p = DimPlot(object = seu_list[[i]], reduction = "FItSNE", label = TRUE, label.size=6,
            ncol=3, combine = FALSE, pt.size = 0.8, cols = hoge1, repel = TRUE)
p = CombinePlots(p) +   theme_void() + theme(legend.position='none')

ggsave(file=paste0(title[i], "_DimPlot.png"),
       plot=p,
       device="png", units="in",
       dpi=600, width=6, height=6, limitsize = FALSE)


i=2
p = DimPlot(object = seu_list[[i]], reduction = "FItSNE", label = TRUE, label.size=6,
            ncol=3, combine = FALSE, split.by = "orig.ident", pt.size = 0.8, cols = hoge1, repel = TRUE)
p = CombinePlots(p) +   theme_void() + theme(legend.position='none')

ggsave(file=paste0(title[i], "_DimPlot.png"),
       plot=p,
       device="png", units="in",
       dpi=600, width=12, height=6, limitsize = FALSE)

i=3
p = DimPlot(object = seu_list[[i]], reduction = "FItSNE", label = TRUE, label.size=6,
            ncol=3, combine = FALSE, split.by = "orig.ident", pt.size = 0.8, cols = hoge1, repel = TRUE)
p = CombinePlots(p) +   theme_void() + theme(legend.position='none')

ggsave(file=paste0(title[i], "_DimPlot.png"),
       plot=p,
       device="png", units="in",
       dpi=600, width=18, height=6, limitsize = FALSE)

i=4
p = DimPlot(object = seu_list[[i]], reduction = "FItSNE", label = TRUE, label.size=6,
            ncol=3, combine = FALSE, split.by = "orig.ident", pt.size = 0.8, cols = hoge1, repel = TRUE)
p = CombinePlots(p) +   theme_void() + theme(legend.position='none')

ggsave(file=paste0(title[i], "_DimPlot.png"),
       plot=p,
       device="png", units="in",
       dpi=600, width=18, height=6, limitsize = FALSE)

i=5
p = DimPlot(object = seu_list[[i]], reduction = "FItSNE", label = TRUE, label.size=6,
            combine = FALSE, pt.size = 0.8, cols = hoge1, repel = TRUE)
p = CombinePlots(p) +   theme_void() + theme(legend.position='none')

ggsave(file=paste0(title[i], "_DimPlot.png"),
       plot=p,
       device="png", units="in",
       dpi=600, width=6, height=6, limitsize = FALSE)

##################
#plot figure 4b (stacking plot)
#merge Seurat object
for(i in c(1:5)){
  if(i>1){
   seu = merge(x=seu, y=seu_list[[i]])
  } else {
    seu = seu_list[[1]]
  }
}

#draw % of total cell bar graph
tmp = table(seu@meta.data$celltype, seu@meta.data$orig.ident)
tmp = as.data.frame(tmp)

tmp1 = tmp %>% tidyr::spread(key = Var2, value=Freq)
colnames(tmp1)[1]="celltype"
hoge = c("doublet", "misc")
tmp1 = tmp1[,!colnames(tmp1) %in% hoge]

tmp_cellcount = tmp1[,c(1:5, 9:11,6:8)]

#cell composition analysis, calculate % of total
rownames(tmp_cellcount)=tmp_cellcount[,1]
tmp_cellcount = tmp_cellcount[,2:ncol(tmp_cellcount)]
nf = 1/colSums(tmp_cellcount)
tmp_cellcount = sweep(tmp_cellcount, 2, nf, "*")

temp_labels <- seu@meta.data %>%
  group_by(orig.ident) %>%
  tally()

  tmp2 = tmp_cellcount %>% tibble::rownames_to_column() %>%
    reshape2::melt(id.vars = 'rowname') %>%
    mutate(rowname = factor(rowname, levels = levels(seu@active.ident)))

colnames(tmp2)[1:2]=c("celltype", "Sample")

p_PercentOfTotal = tmp2 %>%
  ggplot(aes(Sample, value)) +
  geom_bar(aes(fill = celltype), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'celltype', values = custom_colors$discrete) +
  scale_y_continuous(name = '% of total cells', labels = scales::percent_format(), expand = c(0.01,0)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(file = "percentoftotal_comparison.png", plot = p_PercentOfTotal, device="png", units="in", dpi = 300,
       width = 6.5, height = 6, limitsize=FALSE)
