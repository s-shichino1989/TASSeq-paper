
library(Seurat)
library(dplyr)
library(rDBEC)
library(qs)

#rename idents of deep sequencing data

#rename Tcell object
load("./DC_IM_gate.rda")
new.ident = c("DC_cDC1", "DC_cDC2", "DC_cDC2", "DC_cDC2",  "DC_cDC2", "DC_pDC", "Mac_IM", "DC_cDC3")
names(new.ident)=levels(mBC@active.ident)

mBC = RenameIdents(mBC, new.ident)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

DC_IM_Metadata = mBC@meta.data

#rename entire data
mBC = qread("./Seurat/UTlung_shallow_Seurat.qs")

new.ident = read.table("./shallow_annot.txt", row.names=1, header = T, sep="\t", quote="")

tmp = new.ident[,1]
names(tmp)=as.numeric(rownames(new.ident))

mBC = RenameIdents(mBC, tmp)

hoge = mBC@active.ident
hoge = as.character(hoge)
names(hoge)=rownames(mBC@meta.data)

mBC = AddMetaData(mBC, metadata = hoge, col.name="celltype")

#add DC IM cell annotation
hoge = mBC@meta.data

hoge = hoge[!hoge$seurat_clusters %in% c("24"),]
hoge = hoge[,c(1,ncol(hoge))]
hoge1 = DC_IM_Metadata [,c(1,ncol(DC_IM_Metadata ))]

hoge2 = rbind(hoge, hoge1)

hoge3 = hoge2$celltype
names(hoge3)=rownames(hoge2)
mBC = AddMetaData(mBC, metadata = hoge3, col.name="celltype")

#change active.ident
hoge3 = factor(hoge3, levels = sort(unique(hoge3)))
mBC@active.ident = hoge3

#save renamed object
qsave(mBC, file="./Seurat/UTlung_shallow_Seurat_annot.qs", nthreads=12)

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

mBC.markers$dataset = rep("TAS-Seq.shallow", nrow(mBC.markers))

write.table(mBC.markers, file="markers_byCelltype_TAS-Seq-shallow.txt", sep="\t", quote=F, row.names=F, col.names=T)


#create gene expression percentage table
mBC = qread("./Seurat/UTlung_shallow_Seurat_annot.qs")
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

res$dataset = rep("TAS-Seq.shallow", nrow(res))

#filter out genes less than 5% of cells expressed 

write.table(res, file="ALLgenes_byCelltype_TAS-Seq_shallow.txt", sep="\t", quote=F, row.names=F, col.names=T)

##export metadata
hoge = mBC@meta.data
hoge = hoge[,c(1:3, ncol(hoge))]
write.table(hoge, file="TASSeq_shallow_metadata.txt", sep="\t", quote=F, row.names=T, col.names=T)





