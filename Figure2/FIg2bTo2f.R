#create Figure 2b, 2c, 2d, 2f
#load libraries and data
options(stringsAsFactors = FALSE)
suppressWarnings(suppressMessages(source("/datadrive/Rhapsody_analysis/Rscripts/library_source_Seurat.R")))

#Please Place preprocessed Seurat v2.3.4 objects into working directory.

load("./Lung10X.rda")
load("Lung_day00_Rhapsody.rda")
load("Lung_SMARTSeq2.rda")
load("Lung10Xv3.rda")

#set custom color palette
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


#remove doublets
seu_list=list()

mBC = Lung_BDRhapsody
doublets = mBC@meta.data$celltype %in% c("doublet", "misc")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)
seu_list[[1]] = mBC

mBC = Lung_SMARTSeq2
doublets = mBC@meta.data$celltype %in% c("doublet", "misc")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)
seu_list[[2]] = mBC

mBC = Lung10X
doublets = mBC@meta.data$celltype %in% c("doublet", "misc")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)
seu_list[[3]] = mBC

mBC = Lung10Xv3
doublets = mBC@meta.data$celltype %in% c("doublet", "misc")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)
seu_list[[4]] = mBC

#Figure 2b
Lung = MergeSeurat(object1 = Lung10X, object2=Lung_SMARTSeq2, scale.factor = 1000000)
Lung = MergeSeurat(object1 = Lung, object2=Lung_BDRhapsody, scale.factor = 1000000)
Lung = MergeSeurat(object1 = Lung, object2=Lung10Xv3, scale.factor = 1000000)
p = RidgePlot(object = Lung, features.plot = "nGene", group.by="orig.ident", nCol = 1, do.sort = TRUE,
          cols.use = c("#F8766D", "#00BA38","#C77CFF","#619CFF"), do.return = TRUE) +
     scale_x_continuous(breaks=seq(0,12000,by=3000))
ggsave(file = "nGene_RidgePlot.png", plot = p, 
       dpi = 300, width = 4, height = 3)

#Figure 2c (Scatter plot)
Lung = MergeSeurat(object1 = Lung_BDRhapsody, object2=Lung_SMARTSeq2, scale.factor = 1000000)
hoge = Lung@meta.data
hoge$nReads = hoge$nReads/1000000
p = ggplot(hoge, aes(x=nReads, y=nGene, color=orig.ident)) +
  geom_point(size=0.4, alpha=0.3) +
  scale_color_manual(values=c("#F8766D", "#00BA38")) +
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(size=10)) +
  theme(axis.text.x = element_text(size=10, colour = 1), 
        axis.text.y = element_text(size=10, colour = 1),
        plot.title = element_text(hjust=0.5)) +
  xlab("") + ylab("")
ggsave(file = "nGene_nReads_scatter.png", plot = p, 
       dpi = 300, width = 2.5, height = 2.5)

#Figure 2d
hoge1 = Lung@meta.data
hoge1 = hoge1[,c(3,7)]
hoge1$orig.ident=factor(hoge1$orig.ident)
res = hoge1 %>% wilcox_test(nReads.log~orig.ident, exact = TRUE) %>% add_significance() %>% add_xy_position(x="orig.ident")
hoge1  = hoge1 %>% group_by(orig.ident)
p = ggplot(hoge1, aes(x=orig.ident, y=nReads.log)) +
  geom_violin(aes(fill=orig.ident)) +
  geom_boxplot(aes(fill=orig.ident), width=0.2, color="black", alpha=0.2, outlier.shape=NA) +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  ylim(c(3.5, round((max(hoge1$nReads.log)+2), 1)))+
  stat_pvalue_manual(res, tip.length = 0, 
                     y.position =  round((max(hoge1$nReads.log)+1), 1), size = 4) +
  theme_classic()+
  theme(legend.position="none") +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=8, colour = 1),
        plot.title = element_text(hjust=0.5,size=8)) +
  xlab("") + ylab("") 

ggsave(file = "nReadsLog.png", plot = p, 
       dpi = 300, width = 1, height = 2.5)

#Figure 2f
hoge = unique(c(as.character(unique(seu_list[[1]]@meta.data$celltype)),
                as.character(unique(seu_list[[2]]@meta.data$celltype)),
                as.character(unique(seu_list[[3]]@meta.data$celltype)),
                as.character(unique(seu_list[[4]]@meta.data$celltype))))
hoge1 = custom_colors$discrete[1:length(hoge)]
names(hoge1)=hoge

#plotting Figure 2f left (FIt-SNE plot)
names(seu_list)=c("TAS-Seq", "Smart-seq2", "10X Chromium v2", "10X Chromium v3")
p = list()
color = c("#F8766D", "#00BA38", "#619CFF", "#C77CFF")

for(i in c(1:4)){
p1 = DimPlot(object = seu_list[[i]], reduction.use = "FItSNE", do.label = TRUE,label.size = 7,
             do.return = TRUE, group.by = "celltype", vector.friendly = TRUE,
             pt.size = 1.0,
             no.legend=FALSE, cols.use =hoge1) +
  theme(axis.title.x = element_text(size=10, family = "Arial"),
        axis.title.y = element_text(size=10, family = "Arial"),
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"),
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial"))
p[[i]] = p1 + theme_void() + theme(legend.position='none') + ggtitle(names(seu_list[i])) +
  theme(plot.title = element_text(color = color[i])) +
  theme(plot.title = element_text(hjust=0.5))
}

ggsave(file="DimPlot.png",
       plot=plot_grid(plotlist=p, ncol=2, nrow=2),
       device="png", units="in",
       dpi=600, width=6, height=6, limitsize = FALSE)

#plotting Figure 2f right
#cell composition stacking plot

seu1 = MergeSeurat(seu_list[[1]], seu_list[[2]])
seu1 = MergeSeurat(seu1, seu_list[[3]])
seu1 = MergeSeurat(seu1, seu_list[[4]])

tmp = table(seu1@meta.data$celltype, seu1@meta.data$orig.ident)
tmp = as.data.frame(tmp)

tmp1 = tmp %>% tidyr::spread(key = Var2, value=Freq)
colnames(tmp1)[1]="Seurat_Clusters"
hoge = c("doublet", "not-detected", "not")
tmp1 = tmp1[,!colnames(tmp1) %in% hoge]

file.name_table=paste(dir.name.table, sample_name, "_subsetNumber.txt", sep='')
write.table(x = tmp1, file = file.name, row.names=F, sep="\t", quote=F)
tmp_cellcount = tmp1

rownames(tmp1) = tmp1[,1]
tmp1 = tmp1[,2:ncol(tmp1)]

#cell composition analysis, calculate % of total
rownames(tmp_cellcount)=tmp_cellcount[,1]
tmp_cellcount = tmp_cellcount[,2:ncol(tmp_cellcount)]
nf = 1/colSums(tmp_cellcount)
tmp_cellcount = sweep(tmp_cellcount, 2, nf, "*")

temp_labels <- seu1@meta.data %>%
  group_by(orig.ident) %>%
  tally()
colnames(tmp_cellcount)=c("TAS-Seq", "Smart-seq2", "10X Chromium v2", "10X Chromium v3")
temp_labels = c("TAS-Seq", "Smart-seq2", "10X Chromium v2", "10X Chromium v3")

tmp2 = tmp_cellcount %>% tibble::rownames_to_column() %>%
  reshape2::melt(id.vars = 'rowname') %>%
  mutate(rowname = factor(rowname, levels = levels(factor(seu1@meta.data$celltype))))

colnames(tmp2)[1:2]=c("Seurat_clusters", "Sample")

p_PercentOfTotal = tmp2 %>%
  ggplot(aes(Sample, value)) +
  geom_bar(aes(fill = Seurat_clusters), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'Seurat_clusters', values = hoge1) +
  scale_y_continuous(name = '% of total cells', labels = scales::percent_format(), expand = c(0.01,0)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(file = "Allcell_PercentOfTotal.png", plot = p_PercentOfTotal, device="png", units="in", dpi = 600,
       width = 4, height = 5, limitsize=FALSE)


