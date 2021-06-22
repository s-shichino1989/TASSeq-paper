fnames = dir(pattern = "matrix_inflection_")
tablelist = lapply(fnames, tableread_fast_sparse)

fnames = gsub("\\..+$", "", fnames) # remove ".txt" from file names
fnames = gsub("matrix_inflection_", "", fnames)
gc()

# set minimum read depth as 2^7 raw count for backrgound subtraction
# because the rate of mRNA leak (S/N ratio) is approximately 2^7~2^10
names(tablelist)=fnames
DBEC_filter = background_subtraction_Biex(tablelist, min.event=100, minimum.max.expr=7,
                                          species="hsa", min.ave=6.0, AutoThreshold = FALSE,
                                          min.diff=5.5, modelnames="E",
                                          uncert.thre=1, nthreads=32, sample.name=fnames)

DBEC_filter=cbind(DBEC_filter, DBEC_filter)

names(DBEC_filter) = fnames

DBEC_res = apply_DBEC_filter(tablelist, DBEC_filter=DBEC_filter, nthreads=24, sample.name = fnames)
names(DBEC_res) = fnames

sample.name="RAILD1_"

##create Seurat object (Seurat v2 workflow) and annotate by hashtag
for (i in 1:length(DBEC_res)){
  colnames(DBEC_res[[i]]) = paste(names(tablelist)[i], colnames(tablelist[[i]]), sep='_')
}

seu = lapply(DBEC_res, CreateSeuratObject, min.cells = 1, min.genes = 1)
mBC = seu[[1]]

if(length(seu)>1){
  for (i in c(2:length(DBEC_res))){
    mBC = MergeSeurat(mBC, seu[[i]], do.normalize=FALSE, min.cells = 5, min.genes = 500)
  }
}

#add metadata
colnames(mBC@meta.data)[2]="nReads"

mito.genes = grep(pattern = "^MT.", x = rownames(x = mBC@data), value = TRUE)
percent.mito = Matrix::colSums(mBC@raw.data[mito.genes, ])/Matrix::colSums(mBC@raw.data)
mBC = AddMetaData(object = mBC, metadata = percent.mito, col.name = "percent.mito")

ribo.genes = grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = mBC@data), value = TRUE)
percent.ribo = Matrix::colSums(mBC@raw.data[ribo.genes, ])/Matrix::colSums(mBC@raw.data)
mBC = AddMetaData(object = mBC, metadata = percent.ribo, col.name = "percent.ribo")

riboRNA.genes = grep(pattern = "^RNA[[:digit:]]", x = rownames(x = mBC@data), value = TRUE)
percent.riboRNA = Matrix::colSums(mBC@raw.data[riboRNA.genes, ])/Matrix::colSums(mBC@raw.data)
mBC = AddMetaData(object = mBC, metadata = percent.riboRNA, col.name = "percent.riboRNA")



#filter out outliers
mBC = FilterCells(object = mBC, subset.names = c("percent.mito", "nGene"),
                  low.thresholds = c(-Inf, 500), high.thresholds = c(0.4, Inf))

file.name=paste(dir.name, sample.name, "_nGene.png", sep='')
p.nGene = RidgePlot(object = mBC, features.plot = "nGene", group.by="orig.ident", nCol = 1, do.return=TRUE)


p.nGene = VlnPlot(object = mBC, features.plot = "nGene", 
                  group.by="orig.ident", point.size.use = 0, do.return=TRUE) +
  theme_classic() + theme(legend.position = "none")
  

ggsave(file = file.name, plot = p.nGene, device="png", units="in", dpi = 300,
       width = 3, height = 4, limitsize=FALSE)

file.name=paste(dir.name, sample_name, "_nReads_log.png", sep='')
p.nReads = RidgePlot(object = mBC, features.plot = "nReads.log", group.by="orig.ident", nCol = 1, do.return=TRUE)
ggsave(file = file.name, plot = p.nReads, device="png", units="in", dpi = 300,
       width = 8, height = 5, limitsize=FALSE)



# normalizing data
message("Normalizing data...")
mBC = NormalizeData(object = mBC, scale.factor=1000000,display.progress = FALSE)
message("Scaling data...")
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"), do.par=TRUE, num.cores=16, display.progress = FALSE)
mBC = FindVariableGenes(object = mBC, mean.function = ExpMean, dispersion.function = LogVMR,
                        x.low.cutoff = 0.1, x.high.cutoff = Inf, y.cutoff = 0.5, do.plot=FALSE)
hvg.number = length(x = mBC@var.genes)

mBC = RunPCA(object = mBC, pc.genes = mBC@var.genes, do.print = FALSE, pcs.compute = 100, overwrite=TRUE)
mBC = ProjectPCA(object = mBC, do.print = FALSE)
mBC = JackStraw(object = mBC, num.replicate = 100, num.pc = 100,do.par=TRUE, num.cores=8, display.progress = FALSE)
mBC = JackStrawPlot(object = mBC, PCs = 1:100)
quiet(dev.off())
quiet(gc())

#clustering and detection of cellular subsets
tmp = as.data.frame(mBC@dr$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05,1]
if(length(tmp1)<1){
  dims=1:100
} else {
  dims= c(1:(min(tmp1)-1))
}

mBC = FindClusters(object = mBC, reduction.type = "pca", dims.use = dims,
                   resolution = 5.0, print.output = 0, save.SNN = FALSE, force.recalc = TRUE)

py_set_seed(42)
mBC = DoFItSNE(mBC, reduction_use = "pca", dims_use = as.integer(dims),
               perplexity=100, n_jobs=as.integer(16), df = 0.9, random.seed=42)


p1 = DimPlot(object = mBC, reduction.use = "FItSNE", do.label = TRUE,label.size = 10,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1,
             no.legend=FALSE, cols.use =custom_colors$discrete ) +
  theme(axis.title.x = element_text(size=10, family = "Arial"),
        axis.title.y = element_text(size=10, family = "Arial"),
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"),
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial"))

p2 = DimPlot(object = mBC, reduction.use = "FItSNE", do.label = FALSE,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1,group.by = "orig.ident",
             no.legend=FALSE) +
  theme(axis.title.x = element_text(size=10, family = "Arial"),
        axis.title.y = element_text(size=10, family = "Arial"),
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"),
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial"))


p2 = p2 + theme_void() + theme(legend.position='none')

ggsave(file="DimPlot_RAILD1_origident.png",
       plot=plot_grid(p2),
       device="png", units="in",
       dpi=600, width=3, height=3, limitsize = FALSE)


mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 12,
                                               adj.p.val.threshold=0.05)

mBC.markers$cluster = as.numeric(mBC.markers$cluster)
mBC.markers = mBC.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
mBC.markers = as.data.frame(mBC.markers)

quiet(gc())

options(digits=3)
mBC.markers = mBC.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
mBC.markers$avg_logFC = round(mBC.markers$avg_logFC, digits = 3)
mBC.markers$within_avg_exp = round(mBC.markers$within_avg_exp, digits = 3)
mBC.markers$without_avg_exp = round(mBC.markers$without_avg_exp, digits = 3)
mBC.markers$p_val_adj = mBC.markers$p_val_adj
mBC.markers$p_val = mBC.markers$p_val

mBC.markers_write = mBC.markers
mBC.markers_write$cluster=paste("cluster", mBC.markers_write$cluster, sep="")

file.name.marker.table=paste(dir.name.table, sample_name, "_ALLmarkers_minpct0.2_Adj_p0.05.txt", sep='')
fwrite(x=mBC.markers_write, file = file.name.marker.table, row.names=F, col.names=T, sep="\t", quote=F)

markerGene.number = length(unique(mBC.markers$gene))

top20 = mBC.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
top20 = as.data.frame(top20)
top20  = top20 [!duplicated(top20$gene),]
top20 = top20 %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
top20 = as.data.frame(top20)

file.name.heatmap=paste(dir.name, sample.name, "_markers_heatmap.png", sep='')
p_heatmap = DoHeatmap2(object = mBC, genes.use = top20$gene, genes.ident = top20$cluster,
                       slim.col.label = TRUE, remove.key = FALSE, cex.row=3, disp.min = -2.5, disp.max = 2.5)
ggsave(file = file.name.heatmap, plot = p_heatmap, device="png", units="in", dpi = 300,
       width = 15, height = 25, limitsize=FALSE)

tmp = table(mBC@ident, mBC@meta.data$orig.ident)
tmp = as.data.frame(tmp)

tmp1 = tmp %>% tidyr::spread(key = Var2, value=Freq) 
colnames(tmp1)[1]="Seurat_Clusters"
hoge = c("doublet", "not-detected")
tmp1 = tmp1[,!colnames(tmp1) %in% hoge]

write.table(tmp1, "subsetNumber.txt", row.names=F, sep="\t", quote=F)

dir.name="./"
sample.name="TASSeq_"
#plot expression pattern of each marker genes on FIt-SNE plot
tmp = top20$gene
file.name=paste(dir.name, sample.name, "MarkerGenePlots/", sep='')
dir.create(file.name)

FeaturePlot2(object = mBC, features.plot = tmp, cols.use = c("grey", "red"),threads = 6,
             reduction.use = "FItSNE", no.legend = TRUE, 
             pt.size = 0.5, 
             do.return=FALSE, 
             plot.save = TRUE, dir.save = file.name)

hoge = mBC@meta.data$res.5
hoge = as.character(hoge)
hoge[hoge %in% c("1", "8", "12", "29")]="Tcell_CD4T"
hoge[hoge %in% c("3", "6", "7")]="Tcell_CD8T"
hoge[hoge %in% c("24")]="Tcell_Tgd"
hoge[hoge %in% c("30")]="proliferated"
hoge[hoge %in% c("32")]="Mast_cell"
hoge[hoge %in% c("33")]="Epi_ciliated"
hoge[hoge %in% c("13")]="Epi_AT2"
hoge[hoge %in% c("27")]="Epi_KRT5hi"
hoge[hoge %in% c("4", "11", "26", "31")]="Endo"
hoge[hoge %in% c("0", "28")]="FB"
hoge[hoge %in% c("34")]="FB_THBS2hi"
hoge[hoge %in% c("17","18","19", "21")]="Plasma_cell"
hoge[hoge %in% c("2", "5", "10")]="Bcell"
hoge[hoge %in% c("14")]="Mac_AM"
hoge[hoge %in% c("25")]="Mo_classical"
hoge[hoge %in% c("16")]="Mac_IM"
hoge[hoge %in% c("20")]="Mo_nonclassical"
hoge[hoge %in% c("22")]="Neutrophil"
hoge[hoge %in% c("23")]="misc"
hoge[hoge %in% c("9")]="SMC"
hoge[hoge %in% c("15")]="Pericyte"
hoge[hoge %in% c("35")]="Mac_CXCL5hi"
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name = "celltype")

p2 = DimPlot(object = mBC, reduction.use = "FItSNE", do.label = TRUE,
             do.return = TRUE, vector.friendly = TRUE, pt.size = 1,label.size = 6,
             group.by = "celltype",cols.use = custom_colors$discrete,
             no.legend=FALSE) +
  theme(axis.title.x = element_text(size=10, family = "Arial"),
        axis.title.y = element_text(size=10, family = "Arial"),
        axis.text.x = element_text(size=10, colour = 1, family = "Arial"),
        axis.text.y = element_text(size = 10, colour = 1, family = "Arial"))

p2 = p2 + theme_void() + theme(legend.position='none')

ggsave(file="DimPlot_RAILD1_celltype.png",
       plot=plot_grid(p2),
       device="png", units="in",
       dpi=600, width=3, height=3, limitsize = FALSE)


tmp = table(mBC@meta.data$celltype, mBC@meta.data$orig.ident)
tmp = as.data.frame(tmp)

tmp1 = tmp %>% tidyr::spread(key = Var2, value=Freq) 
colnames(tmp1)[1]="Celltype"
hoge = c("doublet", "not-detected")
tmp1 = tmp1[,!colnames(tmp1) %in% hoge]

tmp_cellcount = tmp1

rownames(tmp1) = tmp1[,1]
tmp1 = tmp1[,2:ncol(tmp1)]

#cell composition analysis, calculate % of total
rownames(tmp_cellcount)=tmp_cellcount[,1]
tmp_cellcount = tmp_cellcount[,2:ncol(tmp_cellcount)]
nf = 1/colSums(tmp_cellcount)            
tmp_cellcount = sweep(tmp_cellcount, 2, nf, "*") 

temp_labels <- mBC@meta.data %>%
  group_by(orig.ident) %>%
  tally()

tmp2 = tmp_cellcount %>% tibble::rownames_to_column() %>%
  reshape2::melt(id.vars = 'rowname') %>% 
  mutate(rowname = factor(rowname, levels = levels(mBC@meta.data$celltype)))

colnames(tmp2)[1:2]=c("Celltype", "Sample")

p_PercentOfTotal = tmp2 %>%
  ggplot(aes(Sample, value)) +
  geom_bar(aes(fill = Celltype), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'Celltype', values = custom_colors$discrete) +
  scale_y_continuous(name = '% of total cells', 
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(file = "composition.png", plot = p_PercentOfTotal, 
       device="png", units="in", dpi = 300,
       width = 3.8, height = 3.7, limitsize=FALSE)

save.pigz(mBC, file="RAILD1_reanalysis.rda", n.cores=8)

#plot library metrics of RA-ILD data (Supplementary Figure 3b)
#plot statistics (percent mitochondrial, percent ribosomal protein, percent.ribosomal RNA)
scientific_notation <- function(x) {
  x <- format(x, scientific = TRUE)
  x <- gsub("^(.*)e", "'\\1'e", x)
  x <- gsub("e", "%*%10^", x)
  x <- gsub('\\+', '', x)
  parse(text = x)
}

fuga = data.frame(nGene = mBC@meta.data$nGene, mito.proportion = mBC@meta.data$percent.mito)
p.mito = ggplot(fuga) + 
  geom_hex(aes(nGene, mito.proportion), 
           bins = 100, 
           show.legend = TRUE) +
  ylim(0, 0.8) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) + 
  scale_x_continuous(breaks=seq(0,round(max(mBC@meta.data$nGene), digits=-3),2000))+
  ggtitle("Mitochondrial gene proportion") +
  theme_linedraw() +
  theme(plot.title=element_text(hjust = 0.5), text=element_text(size=8)) + 
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8))

fuga = data.frame(nGene = mBC@meta.data$nGene, Ribosomal.protein.proportion = mBC@meta.data$percent.ribo)
p.ribo = ggplot(fuga) + 
  geom_hex(aes(nGene, Ribosomal.protein.proportion), 
           bins = 100, 
           show.legend = TRUE) +
  ylim(0, 0.4) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) + 
  scale_x_continuous(breaks=seq(0,round(max(mBC@meta.data$nGene), digits=-3),2000))+
  ggtitle("Ribosomal protein gene proportion") +
  theme_linedraw() +
  theme(plot.title=element_text(hjust = 0.5), text=element_text(size=8)) + 
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8))

fuga = data.frame(nGene = mBC@meta.data$nGene, RibosomalRNA.proportion = mBC@meta.data$percent.riboRNA)
p.riboRNA = ggplot(fuga) + 
  geom_hex(aes(nGene, RibosomalRNA.proportion), 
           bins = 100, 
           show.legend = TRUE) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) + 
  scale_x_continuous(breaks=seq(0,round(max(mBC@meta.data$nGene), digits=-3),2000))+
  scale_y_continuous(labels = scientific_notation)+
  ggtitle("RibosomalRNA proportion") +
  theme_linedraw() +
  theme(plot.title=element_text(hjust = 0.5), text=element_text(size=8)) + 
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8))

fuga = data.frame(nGene = mBC@meta.data$nGene, nReads = mBC@meta.data$nReads)
p.Gene = ggplot(fuga) + 
  geom_hex(aes(nReads, nGene), 
           bins = 100, 
           show.legend = TRUE) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) + 
  scale_x_continuous(breaks=seq(0,round(max(mBC@meta.data$nReads), digits=-3),300000))+
  scale_y_continuous(breaks=seq(0,round(max(mBC@meta.data$nGene), digits=-3),2000))+
  ggtitle("Genes-reads distribution") +
  theme_linedraw() +
  theme(plot.title=element_text(hjust = 0.5), text=element_text(size=8)) + 
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8))


p_stats = plot_grid(p.mito, p.ribo, p.riboRNA, p.Gene, ncol=2, nrow=2)
ggsave(file = "RAILD_stats_legend", plot = p_stats, 
       device="png", units="in", dpi = 300,
       width = 7, height = 4, limitsize=FALSE)

