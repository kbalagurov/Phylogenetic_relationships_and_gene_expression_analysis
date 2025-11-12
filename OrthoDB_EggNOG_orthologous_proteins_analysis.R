#Required packages
library(data.table)
library(tidyverse)
library(stringr)
library(rentrez)
library(taxize)
library(ape)
library(ggtree)
library(Biostrings)



setwd("Input the path to your working directory")

getwd()

####################################################################################
#Data from EggNOG  (File name format: GeneID.fa/.txt (CG1602.fa))
#Data from OrthoDB (File name format: FlybaseID_Metazoa.fa/.txt (FBgn0025185_Metazoa.fa))

file_list <- list.files(path="Path_to_directory")
flybase_to_CG_ortho <- fread("Path_to_file") #Table of GeneID and FlybaseID (CG1602 FBgn0025185)
flybase_to_CG_ortho <- as.data.table(flybase_to_CG_ortho)

Nog_list <- file_list[!str_detect(file_list, "[:digit:](?=_Metazoa.txt$)|[:digit:](?=_Metazoa.fa$)")]
Nog_table_list <- Nog_list[str_detect(Nog_list, "\\d.txt")]
fasta_list_NOG <- Nog_list[str_detect(Nog_list, "\\d.fa")]

Ortho_table_list <- file_list[str_detect(file_list, "[:digit:](?=_Metazoa.txt$)")]
fasta_list_ortho <- file_list[str_detect(file_list, "[:digit:](?=_Metazoa.fa$)")]


tmp_tab_nog <- data.table(ID_1=character(),
                          ID_2= character(),
                          organism_name=character(),
                          organism_taxid= character(),
                          other= character(),
                          ID_sp = character(), 
                          interest_gene = character(),
                          XP_NCBI = character(),
                          gene_id = character())


for (i in 1:length(Nog_table_list)) {
  tmp_nog <- fread(Nog_table_list[i])
  setnames(tmp_nog, colnames(tmp_nog), c("ID_1", "ID_2", "organism_name", "organism_taxid", "other"))
  tmp_nog[, ID_sp := gsub('(^\\w+\\s)', "", tmp_nog[, organism_name])]
  tmp_nog[, ID_sp := str_sub(tmp_nog[, ID_sp], 1, 3)]
  tmp_nog[, ID_sp := str_c(str_sub(tmp_nog[, organism_name], 1, 1),".", tmp_nog[,ID_sp])]
  
  tmp_nog[, XP_NCBI := str_extract(tmp_nog[, other], "XP_\\d+\\.1" )]
  
  a <- gsub("(.txt)$", "", Nog_table_list[i])
  tmp_nog[, interest_gene := a]
  tmp_nog[, gene_id := a] 
  tmp_tab_nog <- rbind(tmp_tab_nog, tmp_nog, fill = TRUE)
  
}
id_list <- tmp_tab_nog$XP_NCBI

for (i in 1:length(id_list)) {
  #i <- 1
  if (!is.na(tmp_tab_nog$XP_NCBI[i])) {
    x <- entrez_fetch(tmp_tab_nog$XP_NCBI[i],db = "protein", rettype = "gp", retmode = "xml")
    id_of_gene <- str_extract_all(x, "GeneID:\\d+")[[1]]
    
    num_ids <- length(id_of_gene)
    
    if (num_ids > 1) {
      tmp_tab_nog[i, gene_id := ID_2]
    } else if (num_ids == 1) {
      tmp_tab_nog[i, gene_id := str_extract(id_of_gene[1], "\\d+")]
    } else {
      tmp_tab_nog[i, gene_id := NA]
    }
  }
  
  if (is.na(tmp_tab_nog[i, XP_NCBI])) {
    tmp_tab_nog[i, gene_id := ID_2]
  }
}

tmp_tab_nog[, XP_NCBI := ifelse(!grepl("^XP_", XP_NCBI), ID_1, XP_NCBI)]


################################################################################


tmp_tab_ortho <- data.table(interest_gene=character(),
                            organism_name= character(),
                            organism_taxid= character(),
                            pub_gene_id= character(),
                            corrected_pubgene = character(),
                            description= character(),
                            ID_sp= character())

for (i in 1:length(Ortho_table_list)) {
  tmp_orthodb <- fread(Ortho_table_list[i])
  corrected_pubgene <- str_replace_all(tmp_orthodb[, pub_gene_id], "\\\\", "_")
  corrected_pubgene <- str_replace_all(corrected_pubgene, "^\\w+\\d+;", "")
  tmp_orthodb <- cbind(tmp_orthodb,corrected_pubgene )
  tmp_a <- gsub("(_Metazoa.txt)$", "", Ortho_table_list[i])
  tmp_orthodb[, interest_gene := tmp_a]
  tmp_orthodb[, ID_sp := gsub('(^\\w+\\s)', "", tmp_orthodb[, organism_name])]
  tmp_orthodb[, ID_sp := str_sub(tmp_orthodb[, ID_sp], 1, 3)]
  tmp_orthodb[, ID_sp := str_c(str_sub(tmp_orthodb[, organism_name], 1, 1),".", tmp_orthodb[,ID_sp])]
  tmp_tab_ortho <- rbind(tmp_tab_ortho, tmp_orthodb[, .(interest_gene,corrected_pubgene, organism_name, organism_taxid,
                                                        pub_gene_id, description, ID_sp)])
}  
  for (j in 1:length(tmp_tab_ortho$corrected_pubgene)) {
    query_term <- paste(tmp_tab_ortho[j, 'organism_name'], "AND", tmp_tab_ortho[j, 'corrected_pubgene'])
    ncbi_from_pubog <- entrez_search(db = "gene", term = query_term)
    num_ids <- length(ncbi_from_pubog$ids)
    if (num_ids == 0) {
      tmp_tab_ortho[j, 'corrected_pubgene'] <- tmp_tab_ortho[j, 'pub_gene_id']
      next
    } else if (num_ids > 1) {
      tmp_tab_ortho[j, 'corrected_pubgene'] <- tmp_tab_ortho[j, 'pub_gene_id']
    } else if (num_ids == 1) {
      tmp_tab_ortho[j, 'corrected_pubgene'] <- ncbi_from_pubog$ids[1]
    } else {
      tmp_tab_ortho[j, 'corrected_pubgene'] <- tmp_tab_ortho[j, 'pub_gene_id']
    }
  }

tmp_tab_ortho[, pub_gene_id := str_replace_all(pub_gene_id, '\\\\','')]
  
colnames(tmp_tab_ortho)[5] <- "gene_id"


for (i in 1:nrow(flybase_to_CG_ortho)) {
  #i <- 1
  tmp_tab_ortho[interest_gene==flybase_to_CG_ortho$FlybaseID[i], interest_gene:= flybase_to_CG_ortho$GeneID[i]]
}

#Combine data from OrthoDB and EggNOG

colnames(tmp_tab_nog)[1] <- "pub_gene_id"
c <- colnames(tmp_tab_nog)[c(7,9,6, 3, 4,1 )]
ct <- colnames(tmp_tab_ortho)[c(1,5,7,2,3,4 )]
tmp_tab_ortho1 <- tmp_tab_ortho[, ..ct]
tmp_tab_nog1 <- tmp_tab_nog[, ..c]

tab_NOGortho <- data.table(interest_gene=character(),
                           gene_id= character(),
                           ID_sp= character(),
                           organism_name= character(),
                           organism_taxid = character())
tab_NOGortho <- rbind(tmp_tab_ortho1, tmp_tab_nog1)


fwrite(tab_NOGortho , file = "combined_data_OrthoNOG.txt", append = TRUE)

###################################################
#Rename_fasta

#OrthoDB

general_fasta_ortho <- data.table(name_and_seq=character())

for (j in 1:length(fasta_list_ortho)) {
  #j <- 1
  str <- read_lines(fasta_list_ortho[j], skip_empty_rows = TRUE)
  str <- str_replace_all(str, '\\\\','')
  
  
  b <- fread(Ortho_table_list[j])
  b[, ID := gsub('(^\\w+\\s)', "", b[, organism_name])]
  b[, ID := str_sub(b[, ID], 1, 3)]
  b[, ID := str_c(str_sub(b[, organism_name], 1, 1),".", b[,ID])]
  b[, pub_gene_id := str_replace_all(pub_gene_id, '\\\\','')]
  tmp_a <- gsub("(_Metazoa.txt)$", "", Ortho_table_list[j])
  tmp_b <- flybase_to_CG_ortho[FlybaseID==tmp_a,][1,1]
  b[, interest_gene := tmp_b]
  b[tab_NOGortho, gene_id := i.gene_id, on = "pub_gene_id"]
  tmp_fasta <-  data.table(name_and_seq=character())
  for (i in 1:(length(str))) {
    #i <- 1
    
    if ((i %% 2) != 0) {
      idid <- str_extract(str[i], ('"pub_gene_id":.*,'))
      tmp_idseq <- data.table(pub_gene_id = idid,  pr_sequence = str[i+1])
      tmp_idseq[, pub_gene_id := str_extract(pub_gene_id, (":.*"))]
      tmp_idseq[, pub_gene_id := str_replace_all(pub_gene_id, '"','')]
      tmp_idseq[, pub_gene_id := str_replace_all(pub_gene_id, ',','')]
      tmp_idseq[, pub_gene_id := str_replace_all(pub_gene_id, ':','')]
      tmp_idseq[, pub_gene_id := str_replace_all(pub_gene_id, '\\\\','')]
      tmp_string <- b[pub_gene_id == tmp_idseq[,pub_gene_id], ]
      tmp_idseq[, ':=' ( ID = tmp_string[1,8], gene_id = tmp_string[1,10], interest_gene = tmp_string[1,9])]
      new_name <- paste0(
        '>', tmp_idseq$ID, '_', tmp_idseq$interest_gene,
        '_', tmp_idseq$gene_id, '_', tmp_idseq$pub_gene_id,
        '\n', tmp_idseq$pr_sequence
      )
      
      
      if ((str_detect(str_extract(new_name, ('>.*\\n')), 'NA')) == TRUE ) { #b| (length(grep(sp_or_sp, new_name))==0)) {
        next} else {
          new_name <- strsplit(new_name, "\n")
          tmp_fasta <- rbind(tmp_fasta,new_name) 
        }
    } 
  }
  
  general_fasta_ortho <- rbind(general_fasta_ortho, tmp_fasta) 
  
}

#EggNOG

general_fasta_NOG <- data.table(name_and_seq=character())
for (j in 1:length(fasta_list_NOG)) {
  #j <- 5
  str <- read_lines(fasta_list_NOG[j], skip_empty_rows = TRUE)
  str <- str_replace_all(str, '\\\\','')
  
  
  b <- fread(Nog_table_list[j])
  setnames(b, colnames(b), c("pub_gene_id", "ID_2", "organism_name", "organism_taxid", "other"))
  #b <- tab
  b[, ID := gsub('(^\\w+\\s)', "", b[, organism_name])]
  b[, ID := str_sub(b[, ID], 1, 3)]
  b[, ID := str_c(str_sub(b[, organism_name], 1, 1),".", b[,ID])]
  b[, pub_gene_id := str_replace_all(pub_gene_id, '\\\\','')]
  tmp_a <- gsub("(.txt)$", "", Nog_table_list[j])
  tmp_b <- flybase_to_CG_ortho[GeneID==tmp_a,][1,1]
  b[, interest_gene := tmp_b]
  b[tab_NOGortho, gene_id := i.gene_id, on = "pub_gene_id"]
  #b[tmp_tab_ortho, gene_id := i.gene_id, on = "pub_gene_id"]
  tmp_fasta <-  data.table(name_and_seq=character())
  for (i in 1:(length(str))) {
    #i <- 1
    
    if ((i %% 2) != 0) {
      idid <- str_extract(str[i],  "(?<=>)[^>]+")
      tmp_idseq <- data.table(pub_gene_id = idid,  pr_sequence = str[i+1])
      tmp_string <- b[pub_gene_id == tmp_idseq[,pub_gene_id], ]
      if (nrow(tmp_string) == 0) {
        print(paste(tmp_idseq[,pub_gene_id], "is not in fastalist"))
        next
      }
      tmp_idseq[, ':=' ( ID = tmp_string[1,6], gene_id = tmp_string[1,8], interest_gene = tmp_string[1,7])]
      new_name <- paste0('>',tmp_idseq$ID,'_',tmp_idseq$interest_gene,
                         '_',tmp_idseq$gene_id,'_',tmp_idseq$pub_gene_id,'\n', tmp_idseq$pr_sequence )
      
      if ((str_detect(str_extract(new_name, ('>.*\\n')), 'NA')) == TRUE ) { #b| (length(grep(sp_or_sp, new_name))==0)) {
        next} else {
          new_name <- strsplit(new_name, "\n")
          tmp_fasta <- rbind(tmp_fasta,new_name) 
        }
    } 
  }
  
  general_fasta_NOG <- rbind(general_fasta_NOG, tmp_fasta) 
  
}

general_fasta <- rbind(general_fasta_ortho, general_fasta_NOG) 
fwrite( general_fasta , file = "general_fasta.fa", append = TRUE)

#Manually process general_fasta.fa if needed. Save the processed file as general_fasta_processed.fa

################################################################################
#Specify which FASTA file to use: processed or unprocessed

#fasta_lines <- readLines("general_fasta.fa")
#fasta_lines <- readLines("general_fasta_processed.fa")

header_lines <- fasta_lines[grepl("^>", fasta_lines)]
fasta_ids <- sub("^>", "", header_lines)
fasta_ids_cleaned <- sub("/.*", "", fasta_ids)
fasta_list_cleaned <- data.table(name = fasta_ids_cleaned)

tab_NOGortho_upd <- copy(tab_NOGortho)
tab_NOGortho_upd[, fasta_id := paste0(ID_sp, "_", interest_gene, "_", gene_id, "_", pub_gene_id)]
tab_NOGortho_upd <- tab_NOGortho_upd[fasta_id %in% fasta_list_cleaned$name]

#################################################################################
#Analysis

length(unique(tab_NOGortho_upd$organism_name))

#List of species
species <- (unique(tab_NOGortho_upd$organism_name))
#taxize check specie name
specie_name_check <- gnr_resolve(species)
# retrieve class and order
get_class_and_order <- function(specie_name) {
  tryCatch({
    classification <- taxize::classification(specie_name, db = "ncbi")
    order <- classification[[1]]$name[classification[[1]]$rank == "order"]
    class <- classification[[1]]$name[classification[[1]]$rank == "class"]
    family <- classification[[1]]$name[classification[[1]]$rank == "family"]
    if (length(order) == 0) order <- NA
    if (length(class) == 0) class <- NA
    if (length(family) == 0) family <- NA
    return(list(Order = order, Class = class, Family = family))
  }, error = function(e) {
    return(list(Order = NA, Class = NA, Family = NA)) 
  })
}

classes_and_orders_ncbi <- lapply(species, get_class_and_order)

specie_and_order_ncbi <- data.table(
  Species = species,
  Order = sapply(classes_and_orders_ncbi, `[[`, "Order"),
  Class = sapply(classes_and_orders_ncbi, `[[`, "Class"),
  Family = sapply(classes_and_orders_ncbi, `[[`, "Family")
)

unique_count_orders <- specie_and_order_ncbi[, .(
  Count = .N,      
  Order = unique(Order)  
), by = Family] 

#################################################################################
#Visualization

species_list_simple <- c("Drosophila virilis","Drosophila melanogaster", "Ceratitis capitata","Musca domestica"
                         , "Homo sapiens",
                  "Lucilia sericata", 'Contarinia nasturtii', 'Culex pipiens','Danaus plexippus')

tax_info_simple <- classification(species_list_simple, db = "ncbi")

tree_from_taxonomy <- class2tree(tax_info_simple)
phylo_tree_simple <- tree_from_taxonomy$phylo
plot(phylo_tree_simple, cex = 0.8)
title("Taxonomy-Based Phylogenetic Tree")

ggtree_obj <- ggtree(phylo_tree_simple) + 
  theme_tree2() + 
  geom_tiplab(size = 3, angle = 0, hjust = 0) 
ggsave("phylogenetic_tree.svg", ggtree_obj, width = 20, height = 10, dpi = 300)
plot(ggtree_obj, cex =0.8)
#################################################################################
#Number of each protein in different species

filtered_data1 <- tab_NOGortho_upd[organism_name %in% species_list_simple]
selected_columns1<- c("interest_gene", "organism_name")

filtered_data1 <- filtered_data1[, ..selected_columns1]
unique_count_filtered_data1 <- filtered_data1[, .(
  Protein_Count = .N  
), by = .(organism_name, interest_gene)]


protein_matrix <- dcast(unique_count_filtered_data1, 
                        organism_name ~ interest_gene, value.var = "Protein_Count", fill = 0)
protein_matrix$Total_Proteins <- rowSums(protein_matrix[, -1], na.rm = TRUE)


desired_column_order <- c("organism_name", "CG1602", "CG1603","CG1605", "CG2129", 
                          "CG10959", "CG18262", 'CG8643', 'CG8944', 'Total_Proteins') 
setcolorder(protein_matrix, desired_column_order)


desired_row_order <- c("Drosophila virilis","Drosophila melanogaster", "Ceratitis capitata","Musca domestica", "Lucilia sericata",
                       'Culex pipiens', 'Contarinia nasturtii', 'Danaus plexippus', 'Homo sapiens')

protein_matrix[, organism_name := factor(organism_name, levels = desired_row_order)]
setorder(protein_matrix, organism_name)

protein_long <- melt(protein_matrix, id.vars = "organism_name", 
                     variable.name = "Protein", value.name = "Count")
plot <- ggplot(protein_long, aes(x = Protein, y = organism_name, fill = Count)) +
  geom_tile(color = "white") +  
  geom_text(aes(label = Count), size = 3) +  # Add counts inside the tiles
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_minimal() +  
  labs(title = "Protein Count in Species",
       x = "Protein",
       y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Helvetica"),  
        axis.text.y = element_text(size = 10, family = "Helvetica"),
        axis.title = element_text(size = 12),family = "Helvetica",
        plot.title = element_text(family = "Helvetica")) +
  scale_y_discrete(limits = rev(levels(factor(protein_long$organism_name))))


ggsave("protein_count_heatmap.svg", plot = plot, width = 10, height = 7, dpi = 300)


########################################################################################



