
library(Seurat)
library(dplyr)
library(rDBEC)
library(qs)

#rename idents of deep sequencing data

#rename DC_IM object
load("./cDC1_cDC2_gate.rda")
new.ident = c("DC_cDC1", "DC_cDC2")
names(new.ident)=levels(mBC@active.ident)

mBC = RenameIdents(mBC, new.ident)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

DC_Metadata = mBC@meta.data

#rename entire data

mBC = qread("./Seurat/day00UT_results_Seurat.qs", nthreads=12)
mBC = FindClusters(mBC, resolution=3.5, random.seed=42)
new.ident = read.table("./TASSeq_dataset2_annot.txt", row.names=1, header = T, sep="\t", quote="")

tmp = new.ident[,1]
names(tmp)=as.numeric(rownames(new.ident))

mBC = RenameIdents(mBC, tmp)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

#add DC/IM annotation
hoge = mBC@meta.data

hoge = hoge[!hoge$RNA_snn_res.3.5 %in% "22",]
hoge = hoge[,c(1,24)]
hoge1 = DC_Metadata[,c(1,9)]

hoge2 = rbind(hoge, hoge1)

hoge3 = hoge2$celltype
names(hoge3)=rownames(hoge2)
mBC = AddMetaData(mBC, metadata = hoge3, col.name="celltype")

#change active.ident
hoge3 = factor(hoge3, levels = sort(unique(hoge3)))
mBC@active.ident = hoge3

#save renamed object
qsave(mBC, file="./Seurat/UTlung_TASSeq_dataset2_Seurat_annot.qs", nthreads=12)

#export cell subset number
tmp = table(mBC@meta.data$celltype, mBC@meta.data$orig.ident)
tmp = as.data.frame(tmp)
tmp1 = tmp %>% tidyr::spread(key = Var2, value=Freq) 
colnames(tmp1)[1]="celltype"
hoge = c("doublet", "not-detected")
tmp1 = tmp1[,!colnames(tmp1) %in% hoge]

write.table(tmp1, "./Seurat/Seurat_tables/UTlung_TASSeq_dataset2_results_subsetNumber_annot.txt", row.names=F, sep="\t", quote=F)


