seuratV3ToURD <- function(seurat.object, reduction.copy = c("umap", "tsne"), assay.copy="RNA") {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    
    # Copy over data
    ds@logupx.data <- as(as.matrix(seurat.object@assays[[assay.copy]]@data), "dgCMatrix")
    if(!any(dim(seurat.object@assays[[assay.copy]]@counts) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays[[assay.copy]]@counts[rownames(seurat.object@assays[[assay.copy]]@data), colnames(seurat.object@assays[[assay.copy]]@data)]), "dgCMatrix")
    
    # Copy over metadata
    ## TO DO - grab kmeans clustering info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) { 
      get.data <- as.data.frame(seurat.object@assays[[assay.copy]]@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) { 
      get.data <- as.data.frame(seurat.object@meta.data) 
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    
    # Copy over var.genes
    if(length(seurat.object@assays[[assay.copy]]@var.features > 0)) ds@var.genes <- seurat.object@assays[[assay.copy]]@var.features
    
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    } else if (.hasSlot(seurat.object, "reductions")) {
      reduction.copy <- reduction.copy[which(reduction.copy %in% names(seurat.object@reductions))]
      if (length(reduction.copy) > 0) {
        reduction.copy.present <- names(which(sapply(reduction.copy, function(r) !any(dim(seurat.object@reductions[[r]]) == 0))))
        if (length(reduction.copy.present) > 0)  {
          reduction.copy.present <- reduction.copy.present[1]
          ds@tsne.y <- as.data.frame(seurat.object@reductions[[reduction.copy.present]]@cell.embeddings)
          colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
        }
      }
    }
    
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("pca" %in% names(seurat.object@reductions)) && !any(dim(Loadings(seurat.object, reduction = "pca")) == 0)) {
        ds@pca.load <- as.data.frame(Loadings(seurat.object, reduction = "pca"))
        ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
        ds@pca.sdev <- seurat.object@reductions$pca@stdev
        ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
}