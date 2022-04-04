#perform comparison of nGenes in each cell populations (for Supplementary Figure 6)
library(Seurat)
library(ggplot2)
library(rstatix)
library(dplyr)
library(ggpubr)

#load data

fnames = dir(pattern="metadata.txt")

tablelist = lapply(fnames, read.table, sep="\t", quote="", row.names=1, header=T)

#extract commonly detected cell subsets
cell1 = intersect(unique(tablelist[[1]]$celltype),
                  unique(tablelist[[2]]$celltype))
cell1 = intersect(cell1, unique(tablelist[[3]]$celltype))
cell1 = sort(cell1)
#remove doublets
cell1 = cell1[!cell1 %in% c("doublet", "misc")]

#concatenate metadata table

hoge = rbind(tablelist[[1]], tablelist[[2]], tablelist[[3]])

#extract cell annotation and detected gene number
hoge = hoge[hoge$celltype %in% cell1,]

p.val.res=NULL
plotlist=list()
cell1 = sort(cell1)
cell1 = cell1[c(4:10,15,1,16,17,14,11,12,2,3,13)]


#statistical test of nGene between TAS-Seq and Smart-seq2 for each cell subset by exact wicox_test
for (i in c(1:length(cell1))){
  hoge1 = hoge[hoge$celltype %in%cell1[i],]
  res = hoge1 %>% wilcox_test(nGene~orig.ident, exact = TRUE) %>% add_significance() %>% add_xy_position(x="orig.ident")
  res1 = as.data.frame(res)
  res1$celltype = rep(cell1[i], nrow(res1))
  p.val.res = rbind(p.val.res, res1)
}

#adjust multiple comparison
p.val.res$p.adj.BH = p.adjust(p.val.res$p, method="BH")
p.val.res = p.val.res %>% add_significance(p.col="p.adj.BH")
p.val.res = as.data.frame(p.val.res)

hoge$nGene = hoge$nGene/1000
hoge$nReads = round(log10(hoge$nReads), digits=4)

#create violin plot for each cell subsets
for (i in c(1:length(cell1))){
  hoge1 = hoge[hoge$celltype %in%cell1[i],]
  hoge1$orig.ident=factor(hoge1$orig.ident)
  hoge1  = hoge1 %>% group_by(orig.ident)
  p = ggplot(hoge1, aes(x=orig.ident, y=nGene, group=orig.ident, color=orig.ident)) +
    geom_violin(aes(fill=orig.ident), color="black", alpha=1, scale="width") +
    geom_boxplot(aes(fill=orig.ident), width=0.2, color="black", alpha=0.2, outlier.shape=NA) +
    scale_color_manual(values=c(hue_pal()(6))) +
    scale_fill_manual(values=c(hue_pal()(6))) +
    theme_linedraw() +
    theme(legend.position="none") +
    ggtitle(cell1[i]) +
    scale_y_continuous(limits=c(0, max(c(hoge1$nGene, 10))), label=label_number_si(accuracy = 0.1))+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          plot.title = element_text(hjust=0.5,size=12),
          plot.margin = unit(c(1,2,0,1), "mm")) +
    xlab("") + ylab("")
  plotlist[[i]]=p
}

p_res_violin = plot_grid(plotlist = plotlist, ncol = 5) +
  draw_label(expression(paste("Detected gene number (", {10^3}, ")", sep="")), x=0, y = 0.5, vjust=1, angle=90, size=12)

ggsave(file = "nGene_population.png", plot = p_res_violin, device="png",
       dpi = 300, width = 10, height = 8, units="in", limitsize=FALSE, bg="white")

p.val.res=as.data.frame(p.val.res)
p.val.res=p.val.res[,c(2,3,4,5,6,7,8, 14, 15)]
colnames(p.val.res)[5]="W.statistic"
write.table(p.val.res, file="pvalue_population_nGene_deep.txt",sep="\t", quote=F, row.names=T, col.names=T)

#create nReads violin plot for each cell subsets
for (i in c(1:length(cell1))){
  hoge1 = hoge[hoge$celltype %in%cell1[i],]
  hoge1$orig.ident=factor(hoge1$orig.ident)
  hoge1  = hoge1 %>% group_by(orig.ident)
  p = ggplot(hoge1, aes(x=orig.ident, y=nReads, group=orig.ident, color=orig.ident)) +
    geom_violin(aes(fill=orig.ident), color="black", alpha=1, scale="width") +
    geom_boxplot(aes(fill=orig.ident), width=0.2, color="black", alpha=0.2, outlier.shape=NA) +
    scale_color_manual(values=c(hue_pal()(6))) +
    scale_fill_manual(values=c(hue_pal()(6))) +
    theme_linedraw() +
    theme(legend.position="none") +
    ggtitle(cell1[i]) +
    scale_y_continuous(label=label_number_si(accuracy = 0.1))+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          plot.title = element_text(hjust=0.5,size=12),
          plot.margin = unit(c(1,2,0,1), "mm")) +
    xlab("") + ylab("")
  plotlist[[i]]=p
}

p_res_violin = plot_grid(plotlist = plotlist, ncol = 5) +
  draw_label(expression(paste("Read number (", {log[10]}, ")", sep="")), x=0, y = 0.5, vjust=1, angle=90, size=12)

ggsave(file = "nReads_population.png", plot = p_res_violin, device="png",
       dpi = 300, width = 10, height = 8, units="in", limitsize=FALSE, bg="white")

#wilcox test for nReads
p.val.res=NULL
for (i in c(1:length(cell1))){
  hoge1 = hoge[hoge$celltype %in%cell1[i],]
  res = hoge1 %>% wilcox_test(nReads~orig.ident, exact = TRUE) %>% add_significance() %>% add_xy_position(x="orig.ident")
  res1 = as.data.frame(res)
  res1$celltype = rep(cell1[i], nrow(res1))
  p.val.res = rbind(p.val.res, res1)
}

#adjust multiple comparison
p.val.res$p.adj.BH = p.adjust(p.val.res$p, method="BH")
p.val.res = p.val.res %>% add_significance(p.col="p.adj.BH")
p.val.res = as.data.frame(p.val.res)

p.val.res=as.data.frame(p.val.res)
p.val.res=p.val.res[,c(2,3,4,5,6,7,8, 14, 15)]
colnames(p.val.res)[5]="W.statistic"
write.table(p.val.res, file="pvalue_population_nReads_deep.txt",sep="\t", quote=F, row.names=T, col.names=T)

