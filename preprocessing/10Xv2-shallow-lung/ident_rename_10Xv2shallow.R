
library(Seurat)
library(dplyr)
library(rDBEC)
library(qs)

#rename idents of deep sequencing data

#rename Tcell object
load("./Tcell_gate.rda")
new.ident = c("Tcell_CD4T", "Tcell_CD8T", "Tcell_CD8T", "Tcell_CD8T")
names(new.ident)=levels(mBC@active.ident)

mBC = RenameIdents(mBC, new.ident)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

Tcell_Metadata = mBC@meta.data

#rename gdTcell object
load("./gdTcell_ILC2_gate.rda")
new.ident = c("ILC2", "Tcell_proliferated", "misc", "Tcell_gdT", "misc")
names(new.ident)=levels(mBC@active.ident)

mBC = RenameIdents(mBC, new.ident)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

gdTcell_Metadata = mBC@meta.data


#rename entire data

load("./Seurat/UTlung_10Xv2_shallow_Seurat.rda")

new.ident = read.table("./10Xv2_shallow_annot.txt", row.names=1, header = T, sep="\t", quote="")

tmp = new.ident[,1]
names(tmp)=as.numeric(rownames(new.ident))

mBC = RenameIdents(mBC, tmp)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

#add Tcell annotation
hoge = mBC@meta.data

hoge = hoge[!hoge$seurat_clusters %in% c("4", "20", "8"),]
hoge = hoge[,c(1,23)]
hoge1 = Tcell_Metadata[,c(1,ncol(Tcell_Metadata))]
hoge1_1 = gdTcell_Metadata[,c(1,ncol(gdTcell_Metadata))]

hoge2 = rbind(hoge, hoge1, hoge1_1)

hoge3 = hoge2$celltype
names(hoge3)=rownames(hoge2)
mBC = AddMetaData(mBC, metadata = hoge3, col.name="celltype")

#change active.ident
hoge3 = factor(hoge3, levels = sort(unique(hoge3)))
mBC@active.ident = hoge3

#save renamed object
qsave(mBC, file="./Seurat/UTlung_10Xv2_shallow_Seurat_annot.qs", nthreads=12)


#create gene expression percentage table

norm_expression = mBC@assays$RNA@data
idents = sort(unique(as.character(mBC@active.ident)))
res=NULL

for(i in c(1:length(idents))){
  tmp_matrix = norm_expression[,colnames(norm_expression) %in% rownames(mBC@meta.data[mBC@meta.data$celltype == idents[i],])]
  tmp_matrix = tmp_matrix[rowSums(tmp_matrix)>0,]
  within_avg_exp = apply(tmp_matrix, MARGIN=1, FUN = function(x)(log(x = mean(x = expm1(x = x))+1)))
  pct.1 <- round(
    x = Matrix::rowSums(x = tmp_matrix > 0) /ncol(tmp_matrix),
    digits = 3)
  gene = rownames(tmp_matrix)
  cluster = rep(idents[i], nrow(tmp_matrix))
  tmp = data.frame(gene, cluster,  within_avg_exp, pct.1)
  res = rbind(res, tmp)
}

res$dataset = rep("10Xv2-GSM3926540", nrow(res))

#filter out genes less than 5% of cells expressed 

write.table(res, file="ALLgenes_byCelltype_10Xv2-GSM3926540.txt", sep="\t", quote=F, row.names=F, col.names=T)

##export metadata
hoge = mBC@meta.data
hoge = hoge[,c(1:3, ncol(hoge))]
write.table(hoge, file="10Xv2_shallow_metadata.txt", sep="\t", quote=F, row.names=T, col.names=T)



