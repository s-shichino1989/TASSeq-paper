
library(Seurat)
library(dplyr)
library(rDBEC)
library(qs)

#rename idents of deep sequencing data

#rename Tcell object
load("./Tcell_gate.rda")
new.ident = c("Tcell_Treg", "Tcell_CD8T","Tcell_CD4T", "Tcell_CD8T", "Tcell_CD8T", "Tcell_CD4T","Tcell_CD4T","Tcell_CD8T","Tcell_CD8T")
names(new.ident)=levels(mBC@active.ident)

mBC = RenameIdents(mBC, new.ident)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

Tcell_Metadata = mBC@meta.data

#rename DC object
load("./DC_gate2.rda")
new.ident = c("DC_cDC1", "DC_cDC2","DC_pDC")
names(new.ident)=levels(mBC@active.ident)

mBC = RenameIdents(mBC, new.ident)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

DC_Metadata = mBC@meta.data

#rename IM/Mo object
load("DC_IM_gate1.rda")
new.ident = c("Mo_Ly6Chi","DC", "Mo_Ly6Chi", "Mac_IM")
names(new.ident)=levels(mBC@active.ident)

mBC = RenameIdents(mBC, new.ident)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

mBC = subset(mBC, idents=c("Mo_Ly6Chi", "Mac_IM"))
IM_Metadata = mBC@meta.data

#rename entire data
mBC = qread("./Seurat/UTlung_10Xv2_TabulaMuris_Seurat.qs")

new.ident = read.table("./10Xv2_TabulaMuris_annot.txt", row.names=1, header = T, sep="\t", quote="")

tmp = new.ident[,1]
names(tmp)=as.numeric(rownames(new.ident))

mBC = RenameIdents(mBC, tmp)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

#add Tcell annotation
hoge = mBC@meta.data

hoge = hoge[!hoge$seurat_clusters %in% c("1", "8"),]
hoge = hoge[,c(1,ncol(hoge))]
hoge1 = Tcell_Metadata[,c(1,ncol(Tcell_Metadata))]
hoge1_1 = IM_Metadata[,c(1,ncol(IM_Metadata))]
hoge1_2 = DC_Metadata[,c(1,ncol(DC_Metadata))]

hoge2 = rbind(hoge, hoge1)
hoge2 = rbind(hoge2, hoge1_1)
hoge2 = rbind(hoge2, hoge1_2)

hoge3 = hoge2$celltype
names(hoge3)=rownames(hoge2)
mBC = AddMetaData(mBC, metadata = hoge3, col.name="celltype")

#change active.ident
hoge3 = factor(hoge3, levels = sort(unique(hoge3)))
mBC@active.ident = hoge3

#save renamed object
qsave(mBC, file="./Seurat/UTlung_10Xv2_TabulaMuris_Seurat_annot.qs", nthreads=12)

#perform DE analysis by celltype
mBC.markers = rDBEC::FindMarkers_parallel_lite(mBC, AllcellsIdent=mBC@active.ident,
                                               test.use="wilcox", only.pos=TRUE, min.pct=0.1,
                                               features.use = NULL, nthreads = 16,
                                               adj.p.val.threshold=0.05)

mBC.markers = mBC.markers %>% arrange(desc(avg_logFC))  %>% arrange(cluster)
mBC.markers = as.data.frame(mBC.markers)

quiet(gc())
mBC.markers = mBC.markers[,c(10, 9, 8, 2, 3,4, 6,7,1)]
tmp1 = exp(mBC.markers$avg_logFC)
tmp1 = log2(tmp1)
mBC.markers$avg_logFC = tmp1

#add expression SD info to marker gene table

idents = unique(mBC.markers$cluster)
norm_expression = mBC@assays$RNA@data
res=NULL

for(i in c(1:length(idents))){
  tmp_matrix = norm_expression[mBC.markers[mBC.markers$cluster ==idents[i],"gene"],
                               colnames(norm_expression) %in% rownames(mBC@meta.data[mBC@meta.data$celltype == idents[i],])]
  sd_expression_log = apply(tmp_matrix, MARGIN=1, FUN = function(x)(log(x = sd(x = expm1(x = x)))))
  sd_expression = apply(tmp_matrix, MARGIN=1, FUN = function(x)(log(x = sd(x = expm1(x = x)))))
  tmp = data.frame(sd_expression_log, sd_expression)
  res = rbind(res, tmp)
}

mBC.markers = cbind(mBC.markers, res)

mBC.markers$dataset = rep("10Xv2-TabulaMuris", nrow(mBC.markers))

write.table(mBC.markers, file="markers_byCelltype_10Xv2-TabulaMuris.txt", sep="\t", quote=F, row.names=F, col.names=T)


#create gene expression percentage table
mBC = qread("./Seurat/UTlung_10Xv2_TabulaMuris_Seurat_annot.qs")
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

res$dataset = rep("10Xv2-TabulaMuris", nrow(res))

#filter out genes less than 5% of cells expressed 

write.table(res, file="ALLgenes_byCelltype_10Xv2-TabulaMuris.txt", sep="\t", quote=F, row.names=F, col.names=T)






