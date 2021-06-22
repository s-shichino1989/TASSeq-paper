#perform comparison of nGenes in each cell populations (for Supplementary Figure 3)
library(Seurat)
library(ggplot2)
library(coin)
library(rstatix)
library(dplyr)
library(ggpubr)

#load data

#please copy Seurat v2.3.4 preprocessed datasets into working directory

load("Lung_day00_Rhapsody.rda")
load("Lung_SMARTSeq2.rda")

Lung = MergeSeurat(object1 = Lung_BDRhapsody, object2=Lung_SMARTSeq2, scale.factor = 1000000)

#extract commonly detected cell subsets
cell1 = intersect(as.character(unique(Lung_BDRhapsody@meta.data$celltype)),
                  as.character(unique(Lung_SMARTSeq2@meta.data$celltype)))

#remove doublets
cell1 = cell1[!cell1 %in% c("doublet", "misc")]

#extract cell annotation and detected gene number
hoge = Lung@meta.data[Lung@meta.data$celltype %in% cell1,]
hoge = hoge[,c(1,3, 9)]

p.val.res=NULL
plotlist=list()
cell1 = sort(cell1)
cell1 = cell1[c(4,5,6,7,8,9,10,11,12, 17,18,19,20,1,16,13,14,15,2,3)]


#statistical test of nGene between TAS-Seq and Smart-seq2 for each cell subset by exact wicox_test
for (i in c(1:length(cell1))){
hoge1 = hoge[hoge$celltype %in%cell1[i],c(1,2)]
hoge1$orig.ident=factor(hoge1$orig.ident)
res = hoge1 %>% wilcox_test(nGene~orig.ident) %>% add_significance() %>% add_xy_position(x="orig.ident")
res1 = as.data.frame(res)
p.val.res = rbind(p.val.res, res1)
}

#adjust multiple comparison
p.val.res$p.adj.BH = p.adjust(p.val.res$p, method="BH")
p.val.res = p.val.res %>% add_significance(p.col="p.adj.BH")
p.val.res = as.data.frame(p.val.res)

#create violin plot for each cell subsets
for (i in c(1:length(cell1))){
  hoge1 = hoge[hoge$celltype %in%cell1[i],c(1,2)]
  hoge1$orig.ident=factor(hoge1$orig.ident)
 hoge1  = hoge1 %>% group_by(orig.ident)
 p = ggplot(hoge1, aes(x=orig.ident, y=nGene)) +
  geom_violin(aes(fill=orig.ident)) +
  geom_boxplot(aes(fill=orig.ident), width=0.1, color="black", alpha=0.2, outlier.shape=NA) +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  ylim(c(0, round((max(hoge1$nGene)+2000), -3)))+
  stat_pvalue_manual(p.val.res[i, ,drop=F], tip.length = 0, label = "p.adj.BH.signif",
                     y.position =  round((max(hoge1$nGene)+1000), -3), size = 4) +
  theme_classic()+
  theme(legend.position="none") +
  ggtitle(cell1[i]) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, colour = 1),
        plot.title = element_text(hjust=0.5,size=8)) +
  xlab("") + ylab("")
 plotlist[[i]]=p
}

ggsave(file = "nGene_population.png", plot = plot_grid(plotlist = plotlist, nrow = 5, ncol=4),
       dpi = 300, width = 6, height = 8)

rownames(p.val.res) = cell1
p.val.res=as.data.frame(p.val.res)
p.val.res=p.val.res[,c(2,3,4,5,6,7,8, 13,14)]
colnames(p.val.res)[5]="W.statistic"
write.table(p.val.res, file="pvalue_population_TASseqVSsmartseq2.txt",sep="\t", quote=F, row.names=T, col.names=T)


