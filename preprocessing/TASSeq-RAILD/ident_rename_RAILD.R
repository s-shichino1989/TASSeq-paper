
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

