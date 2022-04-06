
#Draw supplementary Figure 7-9

library(scales)
library(rDBEC)
library(ggplot2)
library(Seurat)
library(rstatix)
library(dplyr)
library(cowplot)
library(biomaRt)
library(ggExtra)
library(cowplot)
library(ggplotify)

db = useMart("ensembl")
mmu = useDataset("mmusculus_gene_ensembl", mart = db)

dir.create("result_tables")

## ALLene comparison

system("wget -L -O ALLgenes_byCelltype_TAS-Seq.deep.txt https://tus.box.com/shared/static/yf4ou6vqfcbmo65t7ntcyujxs04ql827.txt")
system("wget -L -O ALLgenes_byCelltype_Smart-seq2.txt https://tus.box.com/shared/static/vmans6k929ip441dwjol9ou21431bnw4.txt")
system("wget -L -O ALLgenes_byCelltype_10Xv2-TabulaMuris.txt https://tus.box.com/shared/static/jslqxaef9lavlgkjpf9qdxrrf67y0bai.txt")


files = dir(pattern = "ALLgenes_byCelltype_")

ALLlist = lapply(files, read.table, header=T, sep="\t", quote="")


ALL_table=NULL

for(i in c(1:length(ALLlist))){
  ALL_table = rbind(ALL_table, ALLlist[[i]])
}

ALL_table = ALL_table[ALL_table$within_avg_exp>0.10,]
ALL_table = ALL_table[ALL_table$pct.1>0.1,]

##extract commonly detected celltype
dataset = unique(ALL_table$dataset)

p_res1 = list()
p_res2 = list()
p_res3 = list()
wilcox_res = NULL

comparison = c("10Xv2-TabulaMuris", "TAS-Seq.deep")
celltype = intersect(unique(ALLlist[[1]]$cluster), unique(ALLlist[[3]]$cluster))
celltype = celltype[!celltype %in% c("doublet", "misc")]
                  
##ALL gene comparison analysis, 10Xv2 Tabula Muris vs TAS-Seq deep

for(i in c(1:length(celltype))){
expression_table =  ALL_table[ALL_table$cluster == celltype[i],]

tmp = expression_table[expression_table$dataset %in% comparison,]

TASSeq_ALL1 = tmp[tmp$dataset == comparison[2],"gene"]

TenXv2_TabulaMuris_ALL1 =tmp[tmp$dataset == comparison[1], "gene"]

#analysis of commonly-detected genes
ALL_common = intersect(TASSeq_ALL1, TenXv2_TabulaMuris_ALL1)

TASSeq_ALL = tmp[tmp$dataset == comparison[2],]
rownames(TASSeq_ALL)=TASSeq_ALL$gene
TASSeq_ALL = TASSeq_ALL[ALL_common,c("within_avg_exp", "pct.1")]

TenXv2_TabulaMuris_ALL =tmp[tmp$dataset == comparison[1], ]
rownames(TenXv2_TabulaMuris_ALL)=TenXv2_TabulaMuris_ALL$gene
TenXv2_TabulaMuris_ALL =TenXv2_TabulaMuris_ALL[ALL_common,c("within_avg_exp", "pct.1")]

TASSeq_ALL$diff.pct = TASSeq_ALL$pct.1 - TenXv2_TabulaMuris_ALL$pct.1
TASSeq_ALL$group = rep("common", nrow(TASSeq_ALL))

##add TAS-seq only detected genes
ALL_diff = setdiff(TASSeq_ALL1, TenXv2_TabulaMuris_ALL1)
TASSeq_ALL2 = tmp[tmp$dataset == comparison[2],]
rownames(TASSeq_ALL2)=TASSeq_ALL2$gene
TASSeq_ALL2 = TASSeq_ALL2[ALL_diff,c("within_avg_exp", "pct.1")]
TASSeq_ALL2$diff.pct = TASSeq_ALL2$pct.1
TASSeq_ALL2$group = rep(paste(comparison[2], "only", sep=" "), nrow(TASSeq_ALL2))

TASSeq_ALL = rbind(TASSeq_ALL, TASSeq_ALL2)

##add 10Xv2 only detected genes
ALL_diff_10X = setdiff(TenXv2_TabulaMuris_ALL1, TASSeq_ALL1)

if(length(ALL_diff_10X)>0){
 TenXv2_TabulaMuris_ALL2 = tmp[tmp$dataset == comparison[1],]
 rownames(TenXv2_TabulaMuris_ALL2)=TenXv2_TabulaMuris_ALL2$gene
 TenXv2_TabulaMuris_ALL2 = TenXv2_TabulaMuris_ALL2[ALL_diff_10X,c("within_avg_exp", "pct.1")]
 TenXv2_TabulaMuris_ALL2$diff.pct = -TenXv2_TabulaMuris_ALL2$pct.1
 TenXv2_TabulaMuris_ALL2$group = rep(paste(comparison[1], "only", sep=" "), nrow(TenXv2_TabulaMuris_ALL2))
 TASSeq_ALL = rbind(TASSeq_ALL, TenXv2_TabulaMuris_ALL2)
}

if(length(ALL_diff_10X)>0){
p = ggplot(TASSeq_ALL, aes(within_avg_exp, diff.pct, group=factor(group,levels=c(paste(comparison[1], "only", sep=" "),"common",
                                                                                    paste(comparison[2], "only", sep=" "))), 
                              color=factor(group,levels=c(paste(comparison[1], "only", sep=" "),"common",
                                                          paste(comparison[2], "only", sep=" ")))))+
  geom_point(data=TASSeq_ALL, aes(within_avg_exp, diff.pct, group=factor(group,levels=c(paste(comparison[1], "only", sep=" "),"common",
                                                                                           paste(comparison[2], "only", sep=" ")))),
             size=0.3, shape=20) +
  scale_color_manual(values=c("cyan", "green", "magenta")) +
  scale_fill_manual(values=c("cyan", "green", "magenta")) +
  theme_linedraw()+
  ggtitle(celltype[i])+
  ylab("") + 
  xlab("") + 
  theme(plot.title=element_text(hjust = 0.5), 
        text=element_text(size=12), axis.text = element_text(size=12),legend.position = "none",
        plot.margin = unit(c(2,2,2,2), "mm")) 

xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
  geom_boxplot(data=TASSeq_ALL, aes(y = within_avg_exp, 
                              group = factor(group,
                                             levels=c(paste(comparison[1], "only", sep=" "),
                                                      "common",
                                                      paste(comparison[2], "only", sep=" "))), 
                              color =factor(group,
                                            levels=c(paste(comparison[1], "only", sep=" "),
                                                     "common",
                                                     paste(comparison[2], "only", sep=" "))),
                              fill = factor(group,
                                            levels=c(paste(comparison[1], "only", sep=" "),
                                                     "common",
                                                     paste(comparison[2], "only", sep=" ")))),
               outlier.shape=NA, width=0.7)+

  coord_flip()+
  scale_color_manual(values=c("black", "black", "black")) +
  scale_fill_manual(values=c("cyan", "green", "magenta"))
ybox = axis_canvas(p, axis="y", coord_flip = TRUE) +
  geom_boxplot(data=TASSeq_ALL, aes(y = diff.pct, 
                                       group = factor(group,
                                                      levels=c(paste(comparison[1], "only", sep=" "),
                                                               "common",
                                                               paste(comparison[2], "only", sep=" "))), 
                                       color =factor(group,
                                                     levels=c(paste(comparison[1], "only", sep=" "),
                                                              "common",
                                                              paste(comparison[2], "only", sep=" "))),
                                       fill = factor(group,
                                                     levels=c(paste(comparison[1], "only", sep=" "),
                                                              "common",
                                                              paste(comparison[2], "only", sep=" ")))),
               outlier.shape=NA, width=0.7)+
  scale_color_manual(values=c("black", "black", "black")) +
  scale_fill_manual(values=c("cyan", "green", "magenta"))
p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
p2 = insert_yaxis_grob(p1, ybox, grid::unit(0.5, "in"), position="right")
p2 = as.ggplot(p2)

} else {
  p = ggplot(TASSeq_ALL, aes(within_avg_exp, diff.pct, group=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                                                      paste(comparison[2], "high", sep=" "))), 
                                color=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                            paste(comparison[2], "high", sep=" ")))))+
    geom_point(data=TASSeq_ALL, aes(within_avg_exp, diff.pct, group=factor(group,levels=c("common",
                                                                                             paste(comparison[2], "high", sep=" ")))),
               size=0.3, shape=20) +
    scale_color_manual(values=c("green", "magenta")) +
    scale_fill_manual(values=c("green", "magenta")) +
    theme_linedraw()+
    ggtitle(celltype[i])+
    ylab("") + 
    xlab("") + 
    theme(plot.title=element_text(hjust = 0.5), 
          text=element_text(size=12), axis.text = element_text(size=12),legend.position = "none",
          plot.margin = unit(c(3,3,3,3), "mm"))  
  xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
    geom_boxplot(data=TASSeq_ALL, aes(y =  within_avg_exp, 
                                group = factor(group,
                                               levels=c("common",
                                                        paste(comparison[2], "high", sep=" "))), 
                                color =factor(group,
                                              levels=c("common",
                                                       paste(comparison[2], "high", sep=" "))),
                                fill = factor(group,
                                              levels=c("common",
                                                       paste(comparison[2], "high", sep=" ")))),
                 outlier.shape=NA, width=0.7)+
    coord_flip()+
    scale_color_manual(values=c("black", "black")) +
    scale_fill_manual(values=c("green", "magenta"))
  ybox = axis_canvas(p, axis="y", coord_flip = TRUE) +
    geom_boxplot(data=TASSeq_ALL, aes(y =  diff.pct, 
                                         group = factor(group,
                                                        levels=c("common",
                                                                 paste(comparison[2], "high", sep=" "))), 
                                         color =factor(group,
                                                       levels=c("common",
                                                                paste(comparison[2], "high", sep=" "))),
                                         fill = factor(group,
                                                       levels=c("common",
                                                                paste(comparison[2], "high", sep=" ")))),
                 outlier.shape=NA, width=1)+
    scale_color_manual(values=c("black", "black")) +
    scale_fill_manual(values=c("green", "magenta"))
  p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
  p2 = insert_yaxis_grob(p1, ybox, grid::unit(0.7, "in"), position="right")
  p2 = as.ggplot(p2)
}



p_res1[[i]]=p2

filename = paste("./result_tables/", celltype[i], "_", comparison[2], "vs", comparison[1] ,"_ALLGene_comparison.txt", sep="")
write.table(TASSeq_ALL, file = filename, sep="\t", quote=F, row.names=T, col.names=T)

##wilcoxon test
TASSeq_ALL$gene = rownames(TASSeq_ALL)
TASSeq_ALL$within_avg_exp = round(TASSeq_ALL$within_avg_exp, digits=4)
res4_wilcox3 = TASSeq_ALL %>% wilcox_test(within_avg_exp ~ group, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "group")
colnames(res4_wilcox3)[6]="W_statistic"
res4_wilcox3=as.data.frame(res4_wilcox3)
res4_wilcox3=res4_wilcox3[,c(2:8), drop=F]
res4_wilcox3$variable = rep("within_avg_exp", nrow(res4_wilcox3))
res4_wilcox3$celltype = rep(celltype[i], nrow(res4_wilcox3))
  
##percent of GC analysis
  #highly detected genes in TAS-Seq
  res <- getBM(attributes = c("external_gene_name", "percentage_gene_gc_content", "transcript_length"),
               filters = "external_gene_name",
               values = unique(c(rownames(TASSeq_ALL[TASSeq_ALL$diff.pct >=0.10,]), 
                                  rownames(TASSeq_ALL[TASSeq_ALL$group == "TAS-Seq.deep only",]))),
               mart = mmu)
  res$group = rep(paste(comparison[2], "high", sep=" "), nrow(res))
  res$diff.pct = TASSeq_ALL[res$external_gene_name,3]
  
  res2=NULL
  if(sum(TASSeq_ALL$diff.pct) >0){
    res2 <- getBM(attributes = c("external_gene_name", "percentage_gene_gc_content", "transcript_length"),
                filters = "external_gene_name",
                values = unique(c(rownames(TASSeq_ALL[TASSeq_ALL$diff.pct <= (-0.10),]), 
                                          rownames(TASSeq_ALL[TASSeq_ALL$group == paste(comparison[1], "only", sep=" "),]))),
                mart = mmu)
    if(nrow(res2)>0){
    res2$group = rep(paste(comparison[1], "high", sep=" "), nrow(res2))
    res2$diff.pct = TASSeq_ALL[res2$external_gene_name,3]
    } else {
      res2 =NULL
    }
  } else {
    res2 =NULL
  }
  tmp = TASSeq_ALL[TASSeq_ALL$group == "common",]
  res3 <- getBM(attributes = c("external_gene_name", "percentage_gene_gc_content", "transcript_length"),
               filters = "external_gene_name",
               values = as.vector(rownames(tmp[abs(tmp$diff.pct) < 0.10,])),
               mart = mmu)
  res3$group = rep(paste("common"), nrow(res3))
  res3$diff.pct = TASSeq_ALL[res3$external_gene_name,3]
  
  res4 = rbind(res, res2, res3)
  colnames(res4)[1]="gene"
  
  res4$transcript_length = round(log10(res4$transcript_length), digits=4)
  
  filename = paste("./result_tables/", celltype[i],  "_", comparison[2], "vs", comparison[1] , "ALLGene_Attributes_comparison.txt", sep="")
  write.table(res4, file = filename, sep="\t", quote=F, row.names=T, col.names=T)
  
  
  #wilcox test of % GC and transcript length
  res4_wilcox1 = res4 %>% wilcox_test(percentage_gene_gc_content ~ group, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "group")
  colnames(res4_wilcox1)[6]="W_statistic"
  res4_wilcox1=as.data.frame(res4_wilcox1)
  res4_wilcox1=res4_wilcox1[,c(2:8), drop=F]
  res4_wilcox1$variable = rep("percent_GC", nrow(res4_wilcox1))
  res4_wilcox1$celltype = rep(celltype[i], nrow(res4_wilcox1))
  
  res4_wilcox2 = res4 %>% wilcox_test(transcript_length ~ group, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "group")
  colnames(res4_wilcox2)[6]="W_statistic"
  res4_wilcox2=as.data.frame(res4_wilcox2)
  res4_wilcox2=res4_wilcox2[,c(2:8), drop=F]
  res4_wilcox2$variable = rep("transcript_length", nrow(res4_wilcox2))
  res4_wilcox2$celltype = rep(celltype[i], nrow(res4_wilcox2))
  
  res4_wilcox = rbind(res4_wilcox1, res4_wilcox2)
  res4_wilcox = rbind(res4_wilcox,  res4_wilcox3)
  wilcox_res = rbind(wilcox_res, res4_wilcox)
  
  if(!is.null(res2)){
  p = ggplot(res4, aes(percentage_gene_gc_content, diff.pct,  group=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                                                          paste(comparison[2], "high", sep=" "))), 
                       color=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                   paste(comparison[2], "high", sep=" ")))))+
    geom_point(data=res4, aes(percentage_gene_gc_content, diff.pct, 
                              group=factor(group,
                                           levels=c(paste(comparison[1], "high", sep=" "),
                                                    "common",
                                                    paste(comparison[2], "high", sep=" ")))),
                              size=0.3, shape=20) +
    scale_color_manual(values=c("cyan", "green", "magenta")) +
    scale_fill_manual(values=c("cyan", "green", "magenta")) +
    theme_linedraw()+
    ggtitle(celltype[i])+
    ylab("") + 
    xlab("") + 
    theme(plot.title=element_text(hjust = 0.5), 
          text=element_text(size=12), axis.text = element_text(size=12), legend.position = "none",
          plot.margin = unit(c(3,3,3,3), "mm")) 
  
  xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
    geom_boxplot(data=res4, aes(y = percentage_gene_gc_content, 
                                group = factor(group,
                                               levels=c(paste(comparison[1], "high", sep=" "),
                                                        "common",
                                                        paste(comparison[2], "high", sep=" "))), 
                                color =factor(group,
                                              levels=c(paste(comparison[1], "high", sep=" "),
                                                       "common",
                                                       paste(comparison[2], "high", sep=" "))),
                                fill = factor(group,
                                              levels=c(paste(comparison[1], "high", sep=" "),
                                                       "common",
                                                       paste(comparison[2], "high", sep=" ")))),
                 outlier.shape=NA, width=0.7)+
                   coord_flip()+
                   scale_color_manual(values=c("black", "black", "black")) +
                   scale_fill_manual(values=c("cyan", "green", "magenta"))
  p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
  p1 = as.ggplot(p1)


  p_res2[[i]]=p1
  
  p = ggplot(res4, aes(transcript_length, diff.pct, group=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                                                paste(comparison[2], "high", sep=" "))), 
                                                    color=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                                 paste(comparison[2], "high", sep=" ")))))+
    geom_point(data=res4, aes(transcript_length, diff.pct, 
                              group=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                           paste(comparison[2], "high", sep=" ")))), size=0.3, shape=20) +
    scale_color_manual(values=c("cyan", "green", "magenta")) +
    scale_fill_manual(values=c("cyan", "green","magenta")) +
    theme_linedraw()+
    ggtitle(celltype[i])+
    ylab("") + 
    xlab("") + 
    theme(plot.title=element_text(hjust = 0.5), axis.text = element_text(size=12),
          text=element_text(size=12), legend.position = "none",
          plot.margin = unit(c(3,3,3,3), "mm")) +
    guides(color = guide_legend(override.aes = list(size=4, alpha=1)))
  
  xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
    geom_boxplot(data=res4, aes(y =transcript_length, 
                                group = factor(group,
                                               levels=c(paste(comparison[1], "high", sep=" "),
                                                        "common",
                                                        paste(comparison[2], "high", sep=" "))), 
                                color =factor(group,
                                              levels=c(paste(comparison[1], "high", sep=" "),
                                                       "common",
                                                       paste(comparison[2], "high", sep=" "))),
                                fill = factor(group,
                                              levels=c(paste(comparison[1], "high", sep=" "),
                                                       "common",
                                                       paste(comparison[2], "high", sep=" ")))),
                 outlier.shape=NA, width=0.7)+
    coord_flip()+
    scale_color_manual(values=c("black", "black", "black")) +
    scale_fill_manual(values=c("cyan", "green", "magenta"))
  p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
  p1 = as.ggplot(p1)
  
  p_res3[[i]]=p1
  } else {
    p = ggplot(res4, aes(percentage_gene_gc_content, diff.pct,  
                         group=factor(group,levels=c("common",paste(comparison[2], "high", sep=" "))), 
                         color=factor(group,levels=c("common",
                                                     paste(comparison[2], "high", sep=" ")))))+
      geom_point(data=res4, aes(percentage_gene_gc_content, diff.pct, 
                                group=factor(group,levels=c("common",
                                            paste(comparison[2], "high", sep=" ")))), 
                 size=0.3, shape=20) +
      scale_color_manual(values=c("green","magenta")) +
      scale_fill_manual(values=c("green","magenta")) +
      theme_linedraw()+
      ggtitle(celltype[i])+
      ylab("") + 
      xlab("") + 
      theme(plot.title=element_text(hjust = 0.5), axis.text = element_text(size=12),
            text=element_text(size=12), legend.position = "none",
            plot.margin = unit(c(3,3,3,3), "mm")) 
    
      xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
        geom_boxplot(data=res4, aes(y = percentage_gene_gc_content, 
                                    group = factor(group,
                                                   levels=c(paste(comparison[1], "high", sep=" "),
                                                            "common",
                                                            paste(comparison[2], "high", sep=" "))), 
                                    color =factor(group,
                                                  levels=c(paste(comparison[1], "high", sep=" "),
                                                           "common",
                                                           paste(comparison[2], "high", sep=" "))),
                                    fill = factor(group,
                                                  levels=c(paste(comparison[1], "high", sep=" "),
                                                           "common",
                                                           paste(comparison[2], "high", sep=" ")))),
                     outlier.shape=NA, width=0.7)+
        coord_flip()+
        scale_color_manual(values=c("black", "black")) +
        scale_fill_manual(values=c("green", "magenta"))
      p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
      p1 = as.ggplot(p1)
    
    p_res2[[i]]=p
    
    p = ggplot(res4, aes(transcript_length, diff.pct, group=factor(group,levels=c("common",
                                                                                        paste(comparison[2], "high", sep=" "))), 
                               color=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                           paste(comparison[2], "high", sep=" ")))))+
      geom_point(data=res4, aes(transcript_length, diff.pct, 
                                group=factor(group,levels=c("common",
                                                            paste(comparison[2], "high", sep=" ")))), size=0.3, shape=20) +
      scale_color_manual(values=c("green","magenta")) +
      scale_fill_manual(values=c("green","magenta")) +
      theme_linedraw()+
      ggtitle(celltype[i])+
      ylab("") + 
      xlab("") + 
      theme(plot.title=element_text(hjust = 0.5), axis.text = element_text(size=12),
            text=element_text(size=12), legend.position = "none",
            plot.margin = unit(c(3,3,3,3), "mm"))  +
      guides(color = guide_legend(override.aes = list(size=4, alpha=1)))
    
    xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
      geom_boxplot(data=res4, aes(y =transcript_length, 
                                  group = factor(group,
                                                 levels=c("common",
                                                          paste(comparison[2], "high", sep=" "))), 
                                  color =factor(group,
                                                levels=c("common",
                                                         paste(comparison[2], "high", sep=" "))),
                                  fill = factor(group,
                                                levels=c("common",
                                                         paste(comparison[2], "high", sep=" ")))),
                   outlier.shape=NA, width=0.7)+
      coord_flip()+
      scale_color_manual(values=c("black", "black")) +
      scale_fill_manual(values=c("green", "magenta"))
    p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
    p1 = as.ggplot(p1)
    p_res3[[i]]=p1
  }
  
}

p_res1_1 = plot_grid(plotlist = p_res1, ncol = 5) +
          draw_label("Average gene expression (log)", x=0.5, y = 0, hjust=0.5, vjust=-0.5, angle=0, size=12)+
          draw_label("Difference of the % of positive cells (TAS-Seq - 10Xv2)", x=0, y = 0.5, vjust=1.5, angle=90, size=12)

ggsave(file="TASSeqDeepVS10Xv2TabulaMuris_allGenes_expression.png", plot=p_res1_1, device="png", units="in",
       dpi=300, width=15, height=15, limitsize = FALSE, bg="white")

p_res1_2 = plot_grid(plotlist = p_res2, ncol = 5) +
  draw_label("GC content of genes (%)", x=0.5, y = 0, hjust=0.5, vjust=-0.5, angle=0, size=12)+
  draw_label("Difference of the % of positive cells (TAS-Seq - 10Xv2)", x=0, y = 0.5, vjust=1.5, angle=90, size=12)

ggsave(file="TASSeqDeepVS10Xv2TabulaMuris_allGenes_GCcontent.png", plot=p_res1_2, device="png", units="in",
       dpi=300, width=15, height=15, limitsize = FALSE, bg="white")

p_res1_3 = plot_grid(plotlist = p_res3, ncol = 5) +
  draw_label(expression(paste("Transcript length (", {log[10]}, ")", sep="")), x=0.5, y = 0, hjust=0.5, vjust=-0.5, angle=0, size=12)+
  draw_label("Difference of the % of positive cells (TAS-Seq - 10Xv2)", x=0, y = 0.5, vjust=1.5, angle=90, size=12)

ggsave(file="TASSeqDeepVS10Xv2TabulaMuris_allGenes_transcriptLength.png", plot=p_res1_3, device="png", units="in",
       dpi=300, width=15, height=15, limitsize = FALSE, bg="white")

write.table(wilcox_res, file="wilcox_test_statistics_10Xv2_TAS-Seq.txt", sep="\t", quote=F, row.names=F, col.names=T)

###Smart-seq2 vs TAS-Seq deep

p_res1 = list()
p_res2 = list()
p_res3 = list()
wilcox_res = NULL

comparison = c("Smart-seq2", "TAS-Seq.deep")
celltype = intersect(unique(ALLlist[[2]]$cluster), unique(ALLlist[[3]]$cluster))
celltype = celltype[!celltype %in% c("doublet", "misc")]

##ALL gene comparison analysis, 10Xv2 Tabula Muris vs TAS-Seq deep

for(i in c(1:length(celltype))){
  expression_table =  ALL_table[ALL_table$cluster == celltype[i],]
  
  tmp = expression_table[expression_table$dataset %in% comparison,]
  
  TASSeq_ALL1 = tmp[tmp$dataset == comparison[2],"gene"]
  
  TenXv2_TabulaMuris_ALL1 =tmp[tmp$dataset == comparison[1], "gene"]
  
  #analysis of commonly-detected ALL genes
  ALL_common = intersect(TASSeq_ALL1, TenXv2_TabulaMuris_ALL1)
  
  TASSeq_ALL = tmp[tmp$dataset == comparison[2],]
  rownames(TASSeq_ALL)=TASSeq_ALL$gene
  TASSeq_ALL = TASSeq_ALL[ALL_common,c("within_avg_exp", "pct.1")]
  
  TenXv2_TabulaMuris_ALL =tmp[tmp$dataset == comparison[1], ]
  rownames(TenXv2_TabulaMuris_ALL)=TenXv2_TabulaMuris_ALL$gene
  TenXv2_TabulaMuris_ALL =TenXv2_TabulaMuris_ALL[ALL_common,c("within_avg_exp", "pct.1")]
  
  TASSeq_ALL$diff.pct = TASSeq_ALL$pct.1 - TenXv2_TabulaMuris_ALL$pct.1
  TASSeq_ALL$group = rep("common", nrow(TASSeq_ALL))
  
  ##add TAS-seq only detected genes
  ALL_diff = setdiff(TASSeq_ALL1, TenXv2_TabulaMuris_ALL1)
  TASSeq_ALL2 = tmp[tmp$dataset == comparison[2],]
  rownames(TASSeq_ALL2)=TASSeq_ALL2$gene
  TASSeq_ALL2 = TASSeq_ALL2[ALL_diff,c("within_avg_exp", "pct.1")]
  TASSeq_ALL2$diff.pct = TASSeq_ALL2$pct.1
  TASSeq_ALL2$group = rep(paste(comparison[2], "only", sep=" "), nrow(TASSeq_ALL2))
  
  TASSeq_ALL = rbind(TASSeq_ALL, TASSeq_ALL2)
  
  ##add 10Xv2 only detected genes
  ALL_diff_10X = setdiff(TenXv2_TabulaMuris_ALL1, TASSeq_ALL1)
  
  if(length(ALL_diff_10X)>0){
    TenXv2_TabulaMuris_ALL2 = tmp[tmp$dataset == comparison[1],]
    rownames(TenXv2_TabulaMuris_ALL2)=TenXv2_TabulaMuris_ALL2$gene
    TenXv2_TabulaMuris_ALL2 = TenXv2_TabulaMuris_ALL2[ALL_diff_10X,c("within_avg_exp", "pct.1")]
    TenXv2_TabulaMuris_ALL2$diff.pct = -TenXv2_TabulaMuris_ALL2$pct.1
    TenXv2_TabulaMuris_ALL2$group = rep(paste(comparison[1], "only", sep=" "), nrow(TenXv2_TabulaMuris_ALL2))
    TASSeq_ALL = rbind(TASSeq_ALL, TenXv2_TabulaMuris_ALL2)
  }
  
  if(length(ALL_diff_10X)>0){
    p = ggplot(TASSeq_ALL, aes(within_avg_exp, diff.pct, group=factor(group,levels=c(paste(comparison[1], "only", sep=" "),"common",
                                                                                        paste(comparison[2], "only", sep=" "))), 
                                  color=factor(group,levels=c(paste(comparison[1], "only", sep=" "),"common",
                                                              paste(comparison[2], "only", sep=" ")))))+
      geom_point(data=TASSeq_ALL, aes(within_avg_exp, diff.pct, group=factor(group,levels=c(paste(comparison[1], "only", sep=" "),"common",
                                                                                               paste(comparison[2], "only", sep=" ")))),
                 size=0.3, shape=20) +
      scale_color_manual(values=c("cyan", "green", "magenta")) +
      scale_fill_manual(values=c("cyan", "green", "magenta")) +
      theme_linedraw()+
      ggtitle(celltype[i])+
      ylab("") + 
      xlab("") + 
      theme(plot.title=element_text(hjust = 0.5), 
            text=element_text(size=12), axis.text = element_text(size=12),legend.position = "none",
            plot.margin = unit(c(3,3,3,3), "mm")) 
    
    xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
      geom_boxplot(data=TASSeq_ALL, aes(y = within_avg_exp, 
                                           group = factor(group,
                                                          levels=c(paste(comparison[1], "only", sep=" "),
                                                                   "common",
                                                                   paste(comparison[2], "only", sep=" "))), 
                                           color =factor(group,
                                                         levels=c(paste(comparison[1], "only", sep=" "),
                                                                  "common",
                                                                  paste(comparison[2], "only", sep=" "))),
                                           fill = factor(group,
                                                         levels=c(paste(comparison[1], "only", sep=" "),
                                                                  "common",
                                                                  paste(comparison[2], "only", sep=" ")))),
                   outlier.shape=NA, width=0.7)+
      
      coord_flip()+
      scale_color_manual(values=c("black", "black", "black")) +
      scale_fill_manual(values=c("cyan", "green", "magenta"))
    ybox = axis_canvas(p, axis="y", coord_flip = TRUE) +
      geom_boxplot(data=TASSeq_ALL, aes(y = diff.pct, 
                                           group = factor(group,
                                                          levels=c(paste(comparison[1], "only", sep=" "),
                                                                   "common",
                                                                   paste(comparison[2], "only", sep=" "))), 
                                           color =factor(group,
                                                         levels=c(paste(comparison[1], "only", sep=" "),
                                                                  "common",
                                                                  paste(comparison[2], "only", sep=" "))),
                                           fill = factor(group,
                                                         levels=c(paste(comparison[1], "only", sep=" "),
                                                                  "common",
                                                                  paste(comparison[2], "only", sep=" ")))),
                   outlier.shape=NA, width=1)+
      scale_color_manual(values=c("black", "black", "black")) +
      scale_fill_manual(values=c("cyan", "green", "magenta"))
    p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
    p2 = insert_yaxis_grob(p1, ybox, grid::unit(0.7, "in"), position="right")
    p2 = as.ggplot(p2)
    
  } else {
    p = ggplot(TASSeq_ALL, aes(within_avg_exp, diff.pct, group=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                                                        paste(comparison[2], "high", sep=" "))), 
                                  color=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                              paste(comparison[2], "high", sep=" ")))))+
      geom_point(data=TASSeq_ALL, aes(within_avg_exp, diff.pct, group=factor(group,levels=c("common",
                                                                                               paste(comparison[2], "high", sep=" ")))),
                 size=0.3, shape=20) +
      scale_color_manual(values=c("green", "magenta")) +
      scale_fill_manual(values=c("green", "magenta")) +
      theme_linedraw()+
      ggtitle(celltype[i])+
      ylab("") + 
      xlab("") + 
      theme(plot.title=element_text(hjust = 0.5), 
            text=element_text(size=12), axis.text = element_text(size=12),legend.position = "none",
            plot.margin = unit(c(3,3,3,3), "mm"))  
    xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
      geom_boxplot(data=TASSeq_ALL, aes(y =  within_avg_exp, 
                                           group = factor(group,
                                                          levels=c("common",
                                                                   paste(comparison[2], "high", sep=" "))), 
                                           color =factor(group,
                                                         levels=c("common",
                                                                  paste(comparison[2], "high", sep=" "))),
                                           fill = factor(group,
                                                         levels=c("common",
                                                                  paste(comparison[2], "high", sep=" ")))),
                   outlier.shape=NA, width=0.7)+
      coord_flip()+
      scale_color_manual(values=c("black", "black")) +
      scale_fill_manual(values=c("green", "magenta"))
    ybox = axis_canvas(p, axis="y", coord_flip = TRUE) +
      geom_boxplot(data=TASSeq_ALL, aes(y =  diff.pct, 
                                           group = factor(group,
                                                          levels=c("common",
                                                                   paste(comparison[2], "high", sep=" "))), 
                                           color =factor(group,
                                                         levels=c("common",
                                                                  paste(comparison[2], "high", sep=" "))),
                                           fill = factor(group,
                                                         levels=c("common",
                                                                  paste(comparison[2], "high", sep=" ")))),
                   outlier.shape=NA, width=1)+
      scale_color_manual(values=c("black", "black")) +
      scale_fill_manual(values=c("green", "magenta"))
    p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
    p2 = insert_xaxis_grob(p1, ybox, grid::unit(0.7, "in"), position="right")
    p2 = as.ggplot(p2)
  }
  
  p_res1[[i]]=p2
  
  
  filename = paste("./result_tables/", celltype[i], "_", comparison[2], "vs", comparison[1] ,"_ALLGene_comparison.txt", sep="")
  write.table(TASSeq_ALL, file = filename, sep="\t", quote=F, row.names=T, col.names=T)
  
  ##wilcox test
  TASSeq_ALL$gene = rownames(TASSeq_ALL)
  TASSeq_ALL$within_avg_exp = round(TASSeq_ALL$within_avg_exp, digits=4)
  res4_wilcox3 = TASSeq_ALL %>% wilcox_test(within_avg_exp ~ group, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "group")
  colnames(res4_wilcox3)[6]="W_statistic"
  res4_wilcox3=as.data.frame(res4_wilcox3)
  res4_wilcox3=res4_wilcox3[,c(2:8), drop=F]
  res4_wilcox3$variable = rep("within_avg_exp", nrow(res4_wilcox3))
  res4_wilcox3$celltype = rep(celltype[i], nrow(res4_wilcox3))
  
  ##percent of GC analysis
  res <- getBM(attributes = c("external_gene_name", "percentage_gene_gc_content", "transcript_length"),
               filters = "external_gene_name",
               values = unique(c(rownames(TASSeq_ALL[TASSeq_ALL$diff.pct >=0.10,]), 
                                         rownames(TASSeq_ALL[TASSeq_ALL$group == "TAS-Seq.deep only",]))),
               mart = mmu)
  res$group = rep(paste(comparison[2], "high", sep=" "), nrow(res))
  res$diff.pct = TASSeq_ALL[res$external_gene_name,3]
  
  res2=NULL
  if(sum(TASSeq_ALL$diff.pct) >0){
    res2 <- getBM(attributes = c("external_gene_name", "percentage_gene_gc_content", "transcript_length"),
                  filters = "external_gene_name",
                  values = unique(c(rownames(TASSeq_ALL[TASSeq_ALL$diff.pct <= (-0.10),]), 
                                            rownames(TASSeq_ALL[TASSeq_ALL$group == paste(comparison[1], "only", sep=" "),]))),
                  mart = mmu)
    if(nrow(res2)>0){
      res2$group = rep(paste(comparison[1], "high", sep=" "), nrow(res2))
      res2$diff.pct = TASSeq_ALL[res2$external_gene_name,3]
    } else {
      res2 =NULL
    }
  } else {
    res2 =NULL
  }
  tmp = TASSeq_ALL[TASSeq_ALL$group == "common",]
  res3 <- getBM(attributes = c("external_gene_name", "percentage_gene_gc_content", "transcript_length"),
                filters = "external_gene_name",
                values = as.vector(rownames(tmp[abs(tmp$diff.pct) < 0.10,])),
                mart = mmu)
  res3$group = rep(paste("common"), nrow(res3))
  res3$diff.pct = TASSeq_ALL[res3$external_gene_name,3]
  
  
  res4 = rbind(res, res2, res3)
  colnames(res4)[1]="gene"
  
  filename = paste("./result_tables/", celltype[i],  "_", comparison[2], "vs", comparison[1] , "ALLGene_Attributes_comparison.txt", sep="")
  write.table(res4, file = filename, sep="\t", quote=F, row.names=T, col.names=T)
  
  res4$transcript_length = round(log10(res4$transcript_length), digits=4)
  
  #wilcox test of % GC and transcript length
  res4_wilcox1 = res4 %>% wilcox_test(percentage_gene_gc_content ~ group, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "group")
  colnames(res4_wilcox1)[6]="W_statistic"
  res4_wilcox1=as.data.frame(res4_wilcox1)
  res4_wilcox1=res4_wilcox1[,c(2:8), drop=F]
  res4_wilcox1$variable = rep("percent_GC", nrow(res4_wilcox1))
  res4_wilcox1$celltype = rep(celltype[i], nrow(res4_wilcox1))
  
  res4_wilcox2 = res4 %>% wilcox_test(transcript_length ~ group, exact = TRUE) %>% add_significance() %>% add_xy_position(x = "group")
  colnames(res4_wilcox2)[6]="W_statistic"
  res4_wilcox2=as.data.frame(res4_wilcox2)
  res4_wilcox2=res4_wilcox2[,c(2:8), drop=F]
  res4_wilcox2$variable = rep("transcript_length", nrow(res4_wilcox2))
  res4_wilcox2$celltype = rep(celltype[i], nrow(res4_wilcox2))
  
  res4_wilcox = rbind(res4_wilcox1, res4_wilcox2)
  res4_wilcox = rbind(res4_wilcox,  res4_wilcox3)
  wilcox_res = rbind(wilcox_res, res4_wilcox)
  
  if(!is.null(res2)){
    p = ggplot(res4, aes(percentage_gene_gc_content, diff.pct,  group=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                                                            paste(comparison[2], "high", sep=" "))), 
                         color=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                     paste(comparison[2], "high", sep=" ")))))+
      geom_point(data=res4, aes(percentage_gene_gc_content, diff.pct, 
                                group=factor(group,
                                             levels=c(paste(comparison[1], "high", sep=" "),
                                                      "common",
                                                      paste(comparison[2], "high", sep=" ")))),
                 size=0.3, shape=20) +
      scale_color_manual(values=c("cyan", "green", "magenta")) +
      scale_fill_manual(values=c("cyan", "green", "magenta")) +
      theme_linedraw()+
      ggtitle(celltype[i])+
      ylab("") + 
      xlab("") + 
      theme(plot.title=element_text(hjust = 0.5), 
            text=element_text(size=12), axis.text = element_text(size=12), legend.position = "none",
            plot.margin = unit(c(3,3,3,3), "mm")) 
    
    xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
      geom_boxplot(data=res4, aes(y = percentage_gene_gc_content, 
                                  group = factor(group,
                                                 levels=c(paste(comparison[1], "high", sep=" "),
                                                          "common",
                                                          paste(comparison[2], "high", sep=" "))), 
                                  color =factor(group,
                                                levels=c(paste(comparison[1], "high", sep=" "),
                                                         "common",
                                                         paste(comparison[2], "high", sep=" "))),
                                  fill = factor(group,
                                                levels=c(paste(comparison[1], "high", sep=" "),
                                                         "common",
                                                         paste(comparison[2], "high", sep=" ")))),
                   outlier.shape=NA, width=0.7)+
      coord_flip()+
      scale_color_manual(values=c("black", "black", "black")) +
      scale_fill_manual(values=c("cyan", "green", "magenta"))
    p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
    p1 = as.ggplot(p1)
    
    
    p_res2[[i]]=p1
    
    p = ggplot(res4, aes(transcript_length, diff.pct, group=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                                                  paste(comparison[2], "high", sep=" "))), 
                         color=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                     paste(comparison[2], "high", sep=" ")))))+
      geom_point(data=res4, aes(transcript_length, diff.pct, 
                                group=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                            paste(comparison[2], "high", sep=" ")))), size=0.3, shape=20) +
      scale_color_manual(values=c("cyan", "green", "magenta")) +
      scale_fill_manual(values=c("cyan", "green","magenta")) +
      theme_linedraw()+
      ggtitle(celltype[i])+
      ylab("") + 
      xlab("") + 
      theme(plot.title=element_text(hjust = 0.5), axis.text = element_text(size=12),
            text=element_text(size=12), legend.position = "none",
            plot.margin = unit(c(3,3,3,3), "mm"))  +
      guides(color = guide_legend(override.aes = list(size=4, alpha=1)))
    
    xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
      geom_boxplot(data=res4, aes(y =transcript_length, 
                                  group = factor(group,
                                                 levels=c(paste(comparison[1], "high", sep=" "),
                                                          "common",
                                                          paste(comparison[2], "high", sep=" "))), 
                                  color =factor(group,
                                                levels=c(paste(comparison[1], "high", sep=" "),
                                                         "common",
                                                         paste(comparison[2], "high", sep=" "))),
                                  fill = factor(group,
                                                levels=c(paste(comparison[1], "high", sep=" "),
                                                         "common",
                                                         paste(comparison[2], "high", sep=" ")))),
                   outlier.shape=NA, width=0.7)+
      coord_flip()+
      scale_color_manual(values=c("black", "black", "black")) +
      scale_fill_manual(values=c("cyan", "green", "magenta"))
    p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
    p1 = as.ggplot(p1)
    
    p_res3[[i]]=p1
  } else {
    p = ggplot(res4, aes(percentage_gene_gc_content, diff.pct,  
                         group=factor(group,levels=c("common",paste(comparison[2], "high", sep=" "))), 
                         color=factor(group,levels=c("common",
                                                     paste(comparison[2], "high", sep=" ")))))+
      geom_point(data=res4, aes(percentage_gene_gc_content, diff.pct, 
                                group=factor(group,levels=c("common",
                                                            paste(comparison[2], "high", sep=" ")))), 
                 size=0.3, shape=20) +
      scale_color_manual(values=c("green","magenta")) +
      scale_fill_manual(values=c("green","magenta")) +
      theme_linedraw()+
      ggtitle(celltype[i])+
      ylab("") + 
      xlab("") + 
      theme(plot.title=element_text(hjust = 0.5), axis.text = element_text(size=12),
            text=element_text(size=12), legend.position = "none",
            plot.margin = unit(c(3,3,3,3), "mm")) 
    
    xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
      geom_boxplot(data=res4, aes(y = percentage_gene_gc_content, 
                                  group = factor(group,
                                                 levels=c(paste(comparison[1], "high", sep=" "),
                                                          "common",
                                                          paste(comparison[2], "high", sep=" "))), 
                                  color =factor(group,
                                                levels=c(paste(comparison[1], "high", sep=" "),
                                                         "common",
                                                         paste(comparison[2], "high", sep=" "))),
                                  fill = factor(group,
                                                levels=c(paste(comparison[1], "high", sep=" "),
                                                         "common",
                                                         paste(comparison[2], "high", sep=" ")))),
                   outlier.shape=NA, width=0.7)+
      coord_flip()+
      scale_color_manual(values=c("black", "black")) +
      scale_fill_manual(values=c("green", "magenta"))
    p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
    p1 = as.ggplot(p1)
    
    p_res2[[i]]=p
    
    p = ggplot(res4, aes(transcript_length, diff.pct, group=factor(group,levels=c("common",
                                                                                  paste(comparison[2], "high", sep=" "))), 
                         color=factor(group,levels=c(paste(comparison[1], "high", sep=" "),"common",
                                                     paste(comparison[2], "high", sep=" ")))))+
      geom_point(data=res4, aes(transcript_length, diff.pct, 
                                group=factor(group,levels=c("common",
                                                            paste(comparison[2], "high", sep=" ")))), size=0.3, shape=20) +
      scale_color_manual(values=c("green","magenta")) +
      scale_fill_manual(values=c("green","magenta")) +
      theme_linedraw()+
      ggtitle(celltype[i])+
      ylab("") + 
      xlab("") + 
      theme(plot.title=element_text(hjust = 0.5), axis.text = element_text(size=12),
            text=element_text(size=12), legend.position = "none",
            plot.margin = unit(c(3,3,3,3), "mm"))  +
      guides(color = guide_legend(override.aes = list(size=4, alpha=1)))
    
    xbox = axis_canvas(p, axis="x", coord_flip = TRUE) +
      geom_boxplot(data=res4, aes(y =transcript_length, 
                                  group = factor(group,
                                                 levels=c("common",
                                                          paste(comparison[2], "high", sep=" "))), 
                                  color =factor(group,
                                                levels=c("common",
                                                         paste(comparison[2], "high", sep=" "))),
                                  fill = factor(group,
                                                levels=c("common",
                                                         paste(comparison[2], "high", sep=" ")))),
                   outlier.shape=NA, width=0.7)+
      coord_flip()+
      scale_color_manual(values=c("black", "black")) +
      scale_fill_manual(values=c("green", "magenta"))
    p1 = insert_xaxis_grob(p, xbox, grid::unit(0.5, "in"), position="top")
    p1 = as.ggplot(p1)
    p_res3[[i]]=p1
  }
  
}

p_res1_1 = plot_grid(plotlist = p_res1, ncol = 5) +
  draw_label("Average gene expression (log)", x=0.5, y = 0, hjust=0.5, vjust=-0.5, angle=0, size=12)+
  draw_label("Difference of the % of positive cells (TAS-Seq - Smart-seq2)", x=0, y = 0.5, vjust=1.5, angle=90, size=12)

ggsave(file="TASSeqDeepVSsmartseq2_allGenes_expression.png", plot=p_res1_1, device="png", units="in",
       dpi=300, width=18, height=15, limitsize = FALSE, bg="white")

p_res1_2 = plot_grid(plotlist = p_res2, ncol = 5) +
  draw_label("GC content of genes (%)", x=0.5, y = 0, hjust=0.5, vjust=-0.5, angle=0, size=12)+
  draw_label("Difference of the % of positive cells (TAS-Seq - Smart-seq2)", x=0, y = 0.5, vjust=1.5, angle=90, size=12)

ggsave(file="TASSeqDeepVSsmartseq2_allGenes_GCcontent.png", plot=p_res1_2, device="png", units="in",
       dpi=300, width=15, height=15, limitsize = FALSE, bg="white")

p_res1_3 = plot_grid(plotlist = p_res3, ncol = 5) +
  draw_label(expression(paste("Transcript length (", {log[10]}, ")", sep="")), x=0.5, y = 0, hjust=0.5, vjust=-0.5, angle=0, size=12)+
  draw_label("Difference of the % of positive cells (TAS-Seq - Smart-seq2)", x=0, y = 0.5, vjust=1.5, angle=90, size=12)

ggsave(file="TASSeqDeepVSsmartseq2_allGenes_transcriptLength.png", plot=p_res1_3, device="png", units="in",
       dpi=300, width=15, height=15, limitsize = FALSE, bg="white")


write.table(wilcox_res, file="wilcox_test_statistics_SmartSeq2_TAS-Seq.txt", sep="\t", quote=F, row.names=F, col.names=T)


