
library(Seurat)
library(dplyr)
library(rDBEC)
library(qs)

#rename idents of deep sequencing data

#rename MacDC object
load("./Mac_recluster.rda")
MacDC_Metadata = mBC@meta.data

load("./Endo_recluster.rda")
Endo_Metadata = mBC@meta.data

load("./Epi_recluster.rda")
Epi_Metadata = mBC@meta.data

load("./FB_recluster.rda")
FB_Metadata = mBC@meta.data


#rename entire data

mBC = qread("./Seurat/RAILD_combined_Seurat.qs")

new.ident = read.table("./RAILD_annot.txt", row.names=1, header = T, sep="\t", quote="")

tmp = new.ident[,3]
names(tmp)=as.numeric(rownames(new.ident))

mBC = RenameIdents(mBC, tmp)
hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")


#combine annotation
hoge = mBC@meta.data

hoge = hoge[!hoge$RNA_snn_res.2.5 %in% c("3", "6", "8", "11", "12", "14", "16", "20", "21", "22", "25"),]
hoge = hoge[,c(1,ncol(hoge))]
hoge1 = MacDC_Metadata[,c("orig.ident","celltype")]
hoge1_1 = Endo_Metadata[,c("orig.ident","celltype")]
hoge1_2 = Epi_Metadata[,c("orig.ident","celltype")]
hoge1_3 = FB_Metadata[,c("orig.ident","celltype")]

hoge2 = rbind(hoge, hoge1)
hoge2 = rbind(hoge2, hoge1_1)
hoge2 = rbind(hoge2, hoge1_2)
hoge2 = rbind(hoge2, hoge1_3)

hoge3 = hoge2$celltype
names(hoge3)=rownames(hoge2)
mBC = AddMetaData(mBC, metadata = hoge3, col.name="celltype")

#change active.ident
hoge3 = mBC@meta.data$celltype
names(hoge3)=rownames(mBC@meta.data)
hoge3 = factor(hoge3)
mBC@active.ident = hoge3

#save renamed object
qsave(mBC, file="./Seurat/RAILD_combined_Seurat_annot.qs", nthreads=12)

##generate DimPlot

p = DimPlot(object = mBC, group.by = "orig.ident", reduction = "FItSNE", label = FALSE,
            combine = TRUE,pt.size = 0.5) + theme_void() +
    theme(plot.title=element_blank(), legend.position = 'none')

ggsave(filename="RAILD_DimPlot_UMAP.png", plot=p, device="png", width=5, height=5, units="in", limitsize=F)


p = DimPlot(object = mBC, group.by = "celltype", reduction = "FItSNE", label = TRUE, repel = TRUE, label.size = 4,
            combine = TRUE,pt.size = 0.5) + theme_void() +
  scale_color_manual(values = custom_colors$discrete,  labels=levels(mBC@active.ident)) +
  theme(plot.title=element_blank(), legend.position = 'none')

ggsave(filename="RAILD_DimPlot_UMAP_celltype.png", plot=p, device="png", width=5, height=5, units="in", limitsize=F)


##generate stacking plot
tmp = table(mBC@active.ident, mBC@meta.data$orig.ident)
tmp = as.data.frame(tmp)

tmp1 = tmp %>% tidyr::spread(key = Var2, value=Freq) 
colnames(tmp1)[1]="celltype"
tmp1 = tmp1[,!colnames(tmp1) %in% hoge]

tmp_cellcount = tmp1

rownames(tmp1) = tmp1[,1]
tmp1 = tmp1[,2:ncol(tmp1)]

#cell composition analysis, calculate % of total
rownames(tmp_cellcount)=tmp_cellcount[,1]
tmp_cellcount = tmp_cellcount[,2:ncol(tmp_cellcount)]
nf = 1/colSums(tmp_cellcount)            
tmp_cellcount = sweep(tmp_cellcount, 2, nf, "*") 

write.table(tmp_cellcount, file="./RAILD_cellytpe_composition.txt", sep="\t", quote=F, row.names=T, col.names=T)

temp_labels <- mBC@meta.data %>%
  group_by(orig.ident) %>%
  tally()

  tmp2 = tmp_cellcount %>% tibble::rownames_to_column() %>%
    reshape2::melt(id.vars = 'rowname') %>% 
    mutate(rowname = factor(rowname, levels = levels(mBC@active.ident)))

colnames(tmp2)[1:2]=c("celltype", "Sample")

p_PercentOfTotal = tmp2 %>%
  ggplot(aes(Sample, value)) +
  geom_bar(aes(fill = celltype), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'celltype', values = custom_colors$discrete) +
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

ggsave(file = "RAILD_percentofTotal.png", plot = p_PercentOfTotal, device="png", units="in", dpi = 300,
       width = 5, height = 5, limitsize=FALSE)


