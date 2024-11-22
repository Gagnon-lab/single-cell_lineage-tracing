# single-cell_lineage-tracing
This is a pipeline for connecting lineage barcodes to cells from single-cell RNA-sequencing (transcriptome).
### Author: Andy Sposato 
with generous assistance from James Gagnon (`bin/plot_barcode_on_UMAP.R`) and Jeff Farrell (`bin/SeuratV3ToURD.R`)

### Workflow: 
1) generate Seurat object from CellRanger outputs (R)
2) score cells according to differentiation with URD (R)
3) find the most abundant lineage barcode per cell (Python)
4) add best lineage barcodes to transcriptome (R)
5) visualize barcodes (R & Python)

### Steps 

In `scoring-cells-by-differentiation.Rmd`, I generate Seurat objects, run URD on URD objects to score cells according to differentiation, and add the differentiation score to the Seurat objects.

In `find-best-lineage-barcode.ipynb`, I use the `.allReadCounts` and `.stats` from 10X_GESTALT outputs to generate a table that contains 10X barcodes (cell.id) and lineage barcodes (barcode). 

In `add-best-lineage-barcode.Rmd`,  I add the lineage barcode information from `sample_name_bestlineagebarcode.tsv` to the corresponding Seurat object. 

In `visualize-lineage-barcodes.Rmd`, I plot individual lineage barcodes on a UMAP by highlighting the cells containing that lineage barcode. I also use a combination of R and python to wrangle the data so we can examine normalized lineage barcode abundance across differentiation in the form of a line graph. 

### Directory Structure 

```
single-cell_lineage-tracing
  ├── bin
    ├── plot_barcode_on_UMAP.R
    ├── SeuratV3ToURD.R
  ├── data
  	├── 10X_GESTALT_OUTPUT
  		├── unedited1
			├── unedited1.allReadCounts
			├── unedited1.stats
		.
		.
		.
   	├── 18827X1
    		├── filtered_feature_bc_matrix
        		├── barcodes.tsv.gz
        		├── features.tsv.gz
        		├── matrix.mtx.gz
    	.
    	.
    	.
	├── tables
        	├── intermediate tab separated files
    	├── objects (notice this folder is missing from this repository!)
		├── several Seurat objects
  ├── scoring-cells-by-differentiation.Rmd
  ├── find-best-lineage-barcode.ipynb
  ├── add-best-lineage-barcodes.Rmd
  ├── visualize-lineage-barcodes.Rmd
  ├── single-cell_lineage-tracing.Rproj

```

`bin` holds functions that will be using in some of the Rnotebooks. 

`data` holds CellRanger output files in `18827X1` through `18827X4`. `data` also holds GESTALT output files in `10X_GESTALT_OUTPUT`. 
These `data` folders should be fully populated before you even start. You can see what they should contain in this github. 

`objects` and `tables` in `data` contain intermediate files that will be populated as you work through this pipeline. 

To get started with premade Seurat objects, you can the 20 mo. and 22 mo. objects from Andy's Google Drive: 

20 mo. [DOWNLOAD](https://drive.google.com/file/d/18jMrmVg-Rs8qGvF6EO_MvFwRkEab6W47/view?usp=sharing) (edited)

22 mo. [DOWNLOAD](https://drive.google.com/file/d/1rp9jZmYslALasu2Gh-w5HVno8Z-Lo8jb/view?usp=sharing) (unedited)

Be sure to store them in a folder called `objects` within `data`. 

### Installing a copy of this repository
Maybe you're still learning and don't have data ready to analyze. You can practice working with this pipeline by cloning this repository. 
Open up a terminal and clone the repository onto your local computer: 
```
git clone https://github.com/asposato/single-cell_lineage-tracing.git
cd single-cell_lineage-tracing
```
!!! Remember to make a folder for `objects`!!!


