#Required packages
library(data.table)
library(tidyverse)
library(stringr)
library(reutils)
library(readxl)
library(ggplot2)
library(ComplexHeatmap)
library(grid)

setwd("Input_the_path_to_your_working_directory")

getwd()

#####################################################

expression_data_encode_table <- data.table(Gene_Symbol=character(),
                              RPKM_Expression_Value = character(),
                              stage = character())
selected_columns_encode <- c("Gene_Symbol", 
                      "RPKM_Expression_Value","stage")
deve_stages <- c('e00-02','e02-04','embryo 04-06hr','embryo 06-08hr','embryo 08-10hr','embryo 10-12hr',
                 'embryo 12-14hr','embryo 14-16hr','embryo 16-18hr','embryo 18-20hr','embryo 20-22hr',
                 'embryo 22-24hr','L1','L2','L3-12hr','larva L3 puffstage 1-2','larva L3 puffstage 3-6','larva L3 puffstage 7-9',
                 'white prepupa','prepupa 12hr','pupa 1d','pupa 2d','pupa 3d','pupa 4d','adult male 01day','adult male 05day',
                 'adult male 30day','adult female 01day','adult female 05day','adult female 30day')
for (i in c(1:8)){
  #i <- 1
  cg_expr_encode_data <- read_excel("Enter_the_full_path_to_your_file", sheet = i)
  setDT(cg_expr_encode_data)
  cg_expr_encode_data[, 'stage' := as.list(deve_stages)]
  colnames(cg_expr_encode_data)[c(3, 7)] <- selected_columns_encode[1:2]
  filtered_data_encode <- cg_expr_encode_data[, ..selected_columns_encode]
  expression_data_encode_table <- rbind(expression_data_encode_table, filtered_data_encode)
}
expression_data_encode_table <- as.data.table(expression_data_encode_table)
expression_data_encode_table[, RPKM_Expression_Value := as.numeric(RPKM_Expression_Value)]
expression_data_encode_table[, stage := unlist(stage)]


wide_data_encode <- dcast(expression_data_encode_table, 
                   Gene_Symbol ~ stage, 
                   value.var = "RPKM_Expression_Value", 
                   fun.aggregate = mean) 
wide_data_encode <- wide_data_encode[, ..deve_stages ]
Gene_Symbol_encode <- wide_data_encode$Gene_Symbol
colnames(wide_data_encode)
rownames(wide_data_encode) <- Gene_Symbol_encode

###########################################################################
#Normalization


normalizerow <- function(row){
  row <- as.numeric(row)
  if (max(row) == min(row)) {
    return(rep(0, length(row)))
  } else {
    return((row - min (row)) / (max(row) - min(row)))
  }
}

normalized_encode <- t(apply(wide_data_encode,1,normalizerow))
colnames(normalized_encode) <- colnames(wide_data_encode)
rownames(normalized_encode) <- Gene_Symbol_encode
normalized_encode <- as.matrix(normalized_encode)
colnames(normalized_encode)


###########################################################################
#Visualization

heatmap_encode <- Heatmap(
  normalized_encode, 
  name = 'Expression',  
  col = hcl.colors(10, palette = 'Spectral', alpha = NULL, rev = FALSE, fixup = TRUE),
  cluster_rows = FALSE, 
  cluster_columns = FALSE,   
  show_row_names = TRUE, 
  show_column_names = TRUE,   
  row_names_gp = gpar(fontsize = 12), 
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45
)


dev.new(width = 10, height = 8)

draw(heatmap_encode)
dev.off()
tiff("heatmap.tiff", width = 4000, height = 3000, res = 300)
