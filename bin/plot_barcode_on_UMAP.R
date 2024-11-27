# Fun fact: Jamie helped Andy write this
plot_barcode <- function(input_barcode, seurat_object, color) {
  barcode <- unlist(strsplit(input_barcode,split="_"))
  # identify the barcode sites in the input barcode that have the NA wildcard
  list.of.sites <- c("site1","site2","site3","site4","site5","site6","site7","site8","site9","site10")
  exclude <- barcode == "NA"
  # EXCLUDE all sites that are NA from our search
  good.sites <- list.of.sites[!exclude]
  edits <- barcode[!exclude]
  meta = seurat_object@meta.data
  Cells.Sites <- subset(meta, select = list.of.sites)
  # move cells.sites into a temporary data.frame
  q <- Cells.Sites
  for (i in 1:length(good.sites)){
    q <- q %>% filter(.data[[good.sites[[i]]]] == edits[[i]])
  }
  # copy the cell barcode into a column to grab (new cell.id column)
  q$cell.id <- row.names(q)
  # convert cell.id column into a vector (need vector of cells for cells.highlight)
  cells.to.highlight <- q$cell.id
  #print(length(cells.to.highlight))
  title_barcode <- paste(paste(barcode[[1]], barcode[[2]], barcode[[3]], barcode[[4]], barcode[[5]], barcode[[6]], barcode[[7]], barcode[[8]], barcode[[9]], barcode[[10]], sep = "_"))
  UMAP_coordinates <- as.data.frame(seurat_object@reductions[["umap"]]@cell.embeddings)
  
  DimPlot(seurat_object, reduction = 'umap', cells.highlight = cells.to.highlight, sizes.highlight = 1.5, pt.size = 0.5) +  scale_color_manual(labels = c("barcode absent", "barcode present"), values = c('gray', color)) + NoLegend() + labs(title = title_barcode) + theme(plot.title = element_text(hjust = 0.5, size = 10)) + annotate("text", y=min(UMAP_coordinates$UMAP_2)+0.5, x=min(UMAP_coordinates$UMAP_1)+1.5, label= paste(length(cells.to.highlight), " cells", sep = ""), colour = color, size = 6)
}