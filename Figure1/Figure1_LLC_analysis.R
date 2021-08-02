#re-analysis of murine Lewis Lung carcinoma TAS-Seq data
options(stringsAsFactors = FALSE)

#please install Seurat v2.3.4. and rDBEC package before use (could be downloaded from this repository)
suppressWarnings(suppressMessages(source("library_source_Seurat.R")))

tableread_fast_hashtag = function(x, sep="\t", header=TRUE){
  tmp = data.table::fread(x, header=header, sep=sep, quote="")
  tmp = as.data.frame(tmp)
  return(tmp)
}


#Download data
fileUrl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167192/suppl/GSE167192%5FProcessedData%5FFinal%5FMatrix%5FHashtag%2Etxt%2Egz"
download.file(fileUrl, destfile = "./matrix_inflection_LLC.txt.gz", method = "wget")

#load gene expression matrix
tablelist = lapply("matrix_inflection_LLC.txt.gz", 
                   tableread_fast_sparse)
names(tablelist)="LLC"

#remove hashtag pre-annotation
colnames(tablelist[[1]])=gsub(".*_", "", colnames(tablelist[[1]]))

#load hashtag-expression matrix
hashtag_data=lapply("Hashtag_top1M_LLCTag1.txt.gz", tableread_fast_hashtag)
hashtag_data[[1]][,1] = paste(names(tablelist)[1], hashtag_data[[1]][,1], sep='_')
quiet(gc())

##create Seurat object (Seurat v2 workflow) and annotate by hashtag
colnames(tablelist[[1]]) = paste(names(tablelist)[1], colnames(tablelist[[1]]), sep='_')
seu = lapply(tablelist, CreateSeuratObject, min.cells = 5, min.genes = 500)

#Demultiplexing and generate centered log heatmap (Figure 1f)
seu = Demultiplex_DNAtags(hashtag_data[[1]], seu[[1]],
                          scale.factor = max(colSums(hashtag_data[[1]][,2:ncol(hashtag_data[[1]])])),
                          nTags=14, nCells_in_cartridge=nrow(seu[[1]]@meta.data), 
                          sample.name="LLC",
                          dir.name="./")

quiet(file.remove("LLC_demultiplex.txt"))

#export the number of cells which were assigned toeach Hashtag
tmp = table(seu[[1]]@meta.data$TagIDs)
tmp = as.data.frame(tmp)
rownames(tmp)=tmp[,1]
logs_Ncells = as.data.frame(t(tmp))
logs_Ncells = logs_Ncells[2,,drop=F]
rownames(logs_Ncells)="cell number"
write.table(logs_Ncells, "cellnumber.txt", sep="\t", quote=F)

#re-generate Seurat object for demultiplexed sample
raw.data = seu[[1]]@raw.data
tmp = paste(seu[[1]]@meta.data$TagIDs, rownames(seu[[1]]@meta.data), sep="_")
colnames(raw.data)=as.character(tmp)
colnames(raw.data)=gsub("not_detected", "not-detected", colnames(raw.data))
seu=NULL
seu = CreateSeuratObject(raw.data, min.cells=5, min.genes=500)
colnames(seu@meta.data)[2]="nReads"
save.pigz(seu, file="LLC_demultiplex.rda", n.cores=16)

#Ridgeline plot (Figure 1f)
p = RidgePlot(seu, features.plot="nGene", size.title.use = 0)
p = p + theme(title = element_blank(), axis.title.y = element_blank())

ggsave(file="genes_Ridge.png", plot=p, device="png", units="in",
       dpi=300, width=3.5, height=3.5, limitsize = FALSE)







