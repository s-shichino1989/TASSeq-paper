##draw Figure 6

library(Seurat)
library(ggplot2)
library(cowplot)
library(org.Mm.eg.db) # you can use "org.Hs.eg.db" for human samples
library(AnnotationDbi)
library(GO.db)

dir.name="./"
sample.name="test"


##download annotated Seurat objects for mouse lung data
system("wget -L -O UTlung_deep_Seurat_annot.qs https://tus.box.com/shared/static/ohoohhgf3t12tcmupueajbbapu4xdxm6.qs")
system("wget -L -O UTlung_SmartSeq2_Seurat_annot.qs https://tus.box.com/shared/static/ej8tff4k7szdfqsw6niljnzrronnx910.qs")
system("wget -L -O UTlung_10Xv2_TabulaMuris_Seurat_annot.qs https://tus.box.com/shared/static/fqiabrfwxd7y15jgcwlovelhoz2tgz65.qs")


fnames = dir(pattern="UTlung_")
seu_list = lapply(fnames, qread, nthreads=24)


fnames = dir(pattern="UTlung_")
seu_list = lapply(fnames, qread, nthreads=24)

#rename Smart-seq2 orig.ident
seu_list[[3]]@meta.data$orig.ident = rep("Smart-seq2", nrow(seu_list[[5]]@meta.data))

#rename 10X v2 Tabula Muris data
seu_list[[1]]@meta.data$orig.ident = rep("10Xv2_TabulaMuris", nrow(seu_list[[2]]@meta.data))


#remove doublets
for(i in c(1:3)){
  doublet_key = seu_list[[i]]@meta.data$celltype %in% c("doublet", "misc")
  seu_list[[i]]=subset(seu_list[[i]], cells = rownames(seu_list[[i]]@meta.data[!doublet_key,]))
}

#merge data
for(i in c(1:3)){
  if(i>1){
    seu = merge(x=seu, y=seu_list[[i]])
  } else {
    seu = seu_list[[1]]
  }
}

idents = c("10Xv2-TabulaMuris","TAS-Seq.deep", "Smart-seq2")

#plot interleukin genes
for(i in c(1:3)){
interleukins = grep(pattern="Il[1-9].*", x = rownames(seu_list[[i]]@assays$RNA@counts), value=TRUE)
interleukins = interleukins[!interleukins %in% c(grep(pattern="r", x = interleukins, value=TRUE), "IL18bp", "Il1bos", "Il4i1", "Il6st")]
interleukins=sort(interleukins)

p = DotPlot(seu_list[[i]], cols = c("grey", "red"), features=interleukins, scale = FALSE, dot.min = 0.05,
            group.by="celltype") +
   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), ) +
   xlab("") + ylab("") +
   ggtitle(idents[i])+
   theme(plot.title=element_text(hjust = 0.5),
        text=element_text(size=10))

ggsave(filename=paste0(idents[i], "_IL_dotplot.png"), plot = p, device = "png", width = 7, height=8, units = "in", dpi = 300,
       limitsize = FALSE, bg = "white")
}


#extract growth factor genes
data2 = rownames(seu@assays$RNA@counts)
res = select(org.Mm.eg.db, keys = data2, keytype = "SYMBOL",columns = c("SYMBOL", "GOALL"))
data3 = res[res$ONTOLOGYALL %in% "MF",]
key1 = as.character(as.vector(data3[,2]))
key2 = Term(key1)
data4 = transform(data3,Target=key2)

#specify GO term, you can search by using QuickGO
annotation = c("growth factor activity")
data5 = data4[as.character(data4[,5]) %in% annotation,]
annotation_res = unique(data5[,1])

#remove interleukin genes
annotation_res = annotation_res[!annotation_res %in% interleukins]
annotation_res = sort(annotation_res)

#plot GF genes
for(i in c(1:3)){
  GF_genes = annotation_res[annotation_res %in% rownames(seu_list[[i]]@assays$RNA@counts)]

  p = DotPlot(seu_list[[i]], cols = c("grey", "red"), features=GF_genes, scale = FALSE, dot.min = 0.05,
              group.by="celltype") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), ) +
    xlab("") + ylab("") +
    ggtitle(idents[i])+
    theme(plot.title=element_text(hjust = 0.5),
          text=element_text(size=10))

  ggsave(filename=paste0(idents[i], "_GF_dotplot.png"), plot = p, device = "png", width = 20, height=8, units = "in", dpi = 300,
         limitsize = FALSE, bg = "white")
}

