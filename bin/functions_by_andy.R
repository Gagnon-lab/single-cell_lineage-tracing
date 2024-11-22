# function for adding characters to cell IDs from lineage barcode data to match Seurat object cell IDs
# in un1_data (lineage barcode dataframe), the cell-ids are just a string of 16 letters, but in un_meta_1 (transcriptome dataframe), they have a -1_1 at the end. This is because we integrated two channels of 10X to make the unedited object. 
# So we'll add -1_1 to the cell-ids in the un1_data
# un2_data will need the suffix -1_2
# This a function we can use each time
# here x is the cell-id and y is the suffix we're going to add 
add_chars <- function(x, y){
  # remove the tails on cell.ids 
  paste0(substr(x, 1, 16), y, sep = "")
}

# function for normalizing barcode abundance according to how many cells have a lineage barcode within a differentiation range
normalize <- function(num,dividend) {
  num <- (num/dividend)*100
}