# Heatmap

#'
#'This script will make a heatmap of which genes have which GO term.
#'Customization for ....
#'

## Libraries & Functions ----

library(ComplexHeatmap)
library(circlize)       # add for color functions in continuous annotations

#

## Data ----

#Import subset of 25 significant genes from an analysis that have protein ubiquination or dna repair annotations
#Dataframe is a binary table indicating if the gene has that term or not
gene_terms <- read.delim("Example_Data_Files/example_binary_gene_terms.txt")
#Import uses "." instead of space character in column names
colnames(gene_terms) <- gsub("\\.", " ", colnames(gene_terms))

#Import which genes were present in which screen for annotations
gene_screens_neg <- read.delim("Example_Data_Files/example_gene_hits_negative.txt")
gene_screens_pos <- read.delim("Example_Data_Files/example_gene_hits_positive.txt")

#
## Basic Overlap Plot (gene/GO) ----

#Start with removing the DNA repair only genes
gene_terms_sub <- gene_terms[,-7]
#If needed, subset table for only rows/cols with minimum number of annotations
gene_terms_sub <- gene_terms_sub[apply(gene_terms_sub, 1, sum)>=1,]
gene_terms_sub <- gene_terms_sub[,apply(gene_terms_sub, 2, sum)>=1]

hmapInfo1 <- Heatmap(as.matrix(t(gene_terms_sub)),
                     name = 'Enriched terms',
                     col = c('0' = 'white', '1' = 'darkblue'),
                     rect_gp = gpar(col = 'grey85'),
                     
                     cluster_rows = TRUE,
                     show_row_dend = TRUE,
                     row_title = 'GO Terms',
                     row_title_side = 'left',
                     row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                     row_title_rot = 90,
                     show_row_names = TRUE,
                     row_names_rot = 0,
                     row_names_gp = gpar(fontsize = 7, just = 'left'),
                     row_names_side = 'right',
                     row_names_max_width = max_text_width(
                       colnames(gene_terms_sub) #, 
                       #gp = gpar(fontsize = 12)
                     ),
                     
                     cluster_columns = TRUE,
                     show_column_dend = TRUE,
                     column_title = 'Genes',
                     column_title_side = 'bottom',
                     column_title_gp = gpar(fontsize = 10, fontface = 'bold'),
                     column_title_rot = 0,
                     show_column_names = TRUE,
                     column_names_gp = gpar(fontsize = 10),
                     
                     show_heatmap_legend = FALSE
)

draw(hmapInfo1)

#
## Basic Overlap Plot (gene/Screen) FINGERPRINT PLOT ----

#Using "-1" for sensitization and "1" for resistance combine screen hit tables
all_genes <- unique(c(rownames(gene_screens_neg), 
                      rownames(gene_screens_pos)))
gene_hits_combine <- matrix(0, nrow = length(all_genes), ncol = ncol(gene_screens_neg))
rownames(gene_hits_combine) <- all_genes
colnames(gene_hits_combine) <- colnames(gene_screens_neg)
gene_hits_combine <- as.data.frame(gene_hits_combine)

for (row in 1:nrow(gene_hits_combine)) {
  temp_gene <- rownames(gene_hits_combine)[row]
  if (temp_gene %in% rownames(gene_screens_neg)) {
    if (temp_gene %in% rownames(gene_screens_pos)) {
      gene_hits_combine[row,] <- -1*gene_screens_neg[which(rownames(gene_screens_neg) == temp_gene),] + gene_screens_pos[which(rownames(gene_screens_pos) == temp_gene),]
    } else {
      gene_hits_combine[row,] <- -1*gene_screens_neg[which(rownames(gene_screens_neg) == temp_gene),]
    }
  } else {
    gene_hits_combine[row,] <- gene_screens_pos[which(rownames(gene_screens_pos) == temp_gene),]
  }
}

#If needed, subset table for only rows/cols with minimum number of annotations
#e.g. only want to see the genes in more than 1 screen
gene_hits_sub <- gene_hits_combine[apply(gene_hits_combine, 1, sum)>=1,]
#If a screen doesn't have genes after other filters, can remove here
gene_hits_sub <- gene_hits_combine[,apply(gene_hits_combine, 2, sum)>=1]

#Heatmap requires matrix input, so here choose if you want to use gene_hits_combine, or gene_hits_sub for the conversion
gene_hits_mat <- as.matrix(gene_hits_combine)
#OR
#gene_hits_mat <- as.matrix(gene_hits_sub)

#Use unique symbols to define order of gene groupings
row_order_split <- c("a", "a", letters[2:16], "r", "r", "s", "t", "t", "t", "u", "q")

#PLOT
hmapInfo2 <- Heatmap(gene_hits_mat,
                     name = 'Fingerprint Plot',
                     col = c('-1' = '#994F00', '0' = 'gray95', '1' = '#006CD1'),
                     rect_gp = gpar(col = 'gray80'),
                     
              ## ROWS  (Genes) ##
                     cluster_rows = FALSE,
                     #show_row_dend = FALSE,
                  #If including Row Title
                     #row_title = 'Genes',
                     #row_title_side = 'left',
                     #row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                     #row_title_rot = 90,
                  #If splitting rows
                     row_split = row_order_split,
                     row_title = NULL,
                     row_gap = unit(2, "mm"),
                  #Row names
                     show_row_names = TRUE,
                     row_names_rot = 0,
                     row_names_gp = gpar(fontsize = 7, just = 'right'),
                     row_names_side = 'left',
                     #row_names_max_width = max_text_width(
                      # colnames(gene_hits_mat) #, 
                       #gp = gpar(fontsize = 12)
                     #),
                     
              ## COLUMNS (Screens) ##
                     cluster_columns = FALSE,
                     #show_column_dend = TRUE,
                     #column_title = 'Screens',
                     #column_title_side = 'bottom',
                     #column_title_gp = gpar(fontsize = 10, fontface = 'bold'),
                     #column_title_rot = 0,
                     show_column_names = TRUE,
                     column_names_gp = gpar(fontsize = 10),
                     column_names_side = 'top',
                   #Column widths are not individually controllable, but plot dimension in
                     width = unit(2, "cm"),
                     
              ## LEGEND ##
                     #show_heatmap_legend = FALSE,
                     heatmap_legend_param = list(
                       at = c(-1, 0, 1),
                       labels = c("Sensitizing", "", "Resistance"),
                       title = NULL,
                       #title_position = "lefttop-rot",
                       legend_height = unit(4, "cm"))
              
)

draw(hmapInfo2)
