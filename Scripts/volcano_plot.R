# Volcano Plot

#'
#'This script will make a volcano plot of LFC data from MAGeCK scores.
#'Customization for threshold of identification for sensitization and resistance genes.
#'

## Libraries & Functions ----

library(ggplot2)
library(scales)
library(ggrepel)

#Lab imac was unable to update the scales project to include a required function,
# compose_trans()
# so the source code has been copied and imported.
source("Functions/compose_trans.R")

#

## Data ----

#Import example data of output from MAGeCK analysis
# NOTE: This is a real output from an analysis, but with gene names randomized to preserve data sensitivity
example_mageck_data <- read.delim("Example_Data_Files/example_mageck_output.txt")

#Subset MAGeCK data to dataframe with correct plotting information

##Columns "neg.lfc" and "pos.lfc" contain identical information on limit fold change, 
## which is signed according to direction
plot_LFC_data <- data.frame(ID = example_mageck_data$id, LFC = example_mageck_data$neg.lfc)

##MAGeCK separates positive (resistance) and negative (sensitization) scores, 
## so select the correct score according to the sign of LFC
plot_LFC_data$Score <- ifelse(plot_LFC_data$LFC < 0,
                               example_mageck_data$neg.score,
                               example_mageck_data$pos.score)

#

## Plot Significant genes ----

# Plotting Parameters

sc <- 10^-3                    #maximum value of MAGeCK score for colour highlight
lfc <- 1.5                     #minimum absolute value of LFC for colour highlight
sens_col <- "blue"             #color used for sensitizing genes
res_col <- "mediumseagreen"    #color used for resistance genes
arrow_y <- round(log10(min(plot_LFC_data$Score)))-1    #exp for y value for Sens/Res arrow annotation
#If getting warning
## "ggrepel: 1 unlabeled data points (too many overlaps). Consider increasing max.overlaps"
# Consider increasing this:
m.o <- 15

# Plot

ggplot(data = plot_LFC_data,                                                      #Start with basic color mapping
       mapping = aes(x=LFC, y=Score, 
                     color = ifelse((Score>sc) | (abs(LFC)<lfc), 'darkgray',      #color for not-sig genes
                                    ifelse(LFC>lfc, res_col, sens_col)))) +   #colors for Res, Sens genes
  
  #Plot type is geom_point for scatterplot
  geom_point() +
  
  #Label assigned colors
  scale_colour_manual(labels = c("Sensitizing", "", "Resistance"),
                      values=c(sens_col, 'darkgray', res_col)) +
  
  #Provide axis and plot title labels
  labs(y = "MAGeCK score",x = "Limit fold change", title = "Volcano Plot of Gene Significance") +
  
  #Remove legend (labels as annotation instead), set plot background to white, retain axis labels
  theme(legend.position = "none", 
        panel.background = element_rect(fill = "white"), 
        panel.border = element_blank(), axis.line = element_line()) +
  
  #Set y axis as log10 and reverse order, with breaks every 10^x
  ## Expand breaks as needed
  scale_y_continuous(trans = compose_trans("log10", "reverse"), 
                     breaks = c(1,10^-2,10^-4,10^-6,10^-8,10^-10)) +
  
  #Set x axis with buffer for labels
  xlim(min(plot_LFC_data$LFC)-1, max(plot_LFC_data$LFC)+1) +
  
  #Add labels for significant genes
  ## Sensitizing labels
  geom_text_repel(aes(label=ifelse((LFC<(-lfc)) & (Score<sc), ID,'')),
                  min.segment.length = 0, size = 3, max.overlaps = m.o) +
  ## Resistance labels
  geom_text_repel(aes(label=ifelse((LFC>lfc) & (Score<sc), ID,'')), 
                  min.segment.length = 0, size = 3, max.overlaps = m.o) +
  
  #Add annotation arrows for direction and position of sensitizing and resistance genes
  ## Sensitizing
  annotate("segment", x=0, xend=min(plot_LFC_data$LFC)-0.8, 
           y=10^arrow_y, yend=10^arrow_y, 
           color = sens_col, 
           arrow = arrow(type = "closed", angle = 15)) +
  annotate("text", x = min(plot_LFC_data$LFC)/2, y = 11^arrow_y, 
           label = "Sensitizing", color = sens_col) +
  ## Resistance
  annotate("segment", x=0, xend=max(plot_LFC_data$LFC)+0.8, y=10^arrow_y, yend=10^arrow_y, 
           color = res_col, 
           arrow = arrow(type = "closed", angle = 15)) +
  annotate("text", x = max(plot_LFC_data$LFC)/2, y = 11^arrow_y, 
           label = "Resistance", color = res_col)

#

## Plot Specific Genes ----

# Plotting Parameters

hc <- "mediumvioletred"     #highlight color for genes of interest
gene_desc <- "Genes of Interest"    #label to use for genes of interest, e.g. "Glycosylases"

## For example, will take random selection of k genes
k <- 15
#Choosing from resistance genes for visual
set.seed(volcano)
genes_label <- example_mageck_data[sample(17000:nrow(example_mageck_data), k),"id"]
## If supplied with list of genes, verify all genes in dataset
length(which(genes_label %in% plot_LFC_data$ID))

## ggplot layers the points by row order, so ensure genes to highlight get plotted last so they will be visible
index_in <- which(plot_LFC_data$ID %in% genes_label)
index_out <- which(!(plot_LFC_data$ID %in% genes_label))

plot_LFC_ordered <- plot_LFC_data[c(index_out, index_in),]

# Plot

ggplot(data = plot_LFC_ordered, mapping = aes(x = LFC, y = Score,
                                              color = ifelse(ID %in% genes_label, hc,"darkgray"))) +
  
  #Plot type is geom_point for scatterplot
  geom_point() +
  
  #Label assigned colors
  scale_colour_manual(labels = c(paste0("Genes not in\n", gene_desc, " set"), gene_desc),
                      values=c("gray", hc)) +
  
  #Provide axis, plot title and colour legend labels
  labs(y = "MAGeCK score",x = "Limit fold change", 
       title = paste0("Volcano Plot Highlighting ", gene_desc),
       color = "Gene Colour") +
  
  #Set legend, plot background to white, retain axis labels
  ## To use only axis lines, copy from above
  theme(legend.position = "right",
        legend.box.background = element_rect(color="black"),
        panel.background = element_rect(fill = "white", color = "black")) +
  
  ## UPDATE THIS PART ##
  # Currently there is something about the text labels for the color in the legend that adds a small letter 'a' on top of the coloured point. Currently covering it up with a larger point for the color legend. If I can figure out the parameter that is adding the shape to the legend in the first place, I can just remove it.
  guides(color = guide_legend(override.aes = list(size = 3) ) ) +
  ## --- ##
  
  #Set y axis as log10 and reverse order, with breaks every 10^x
  ## Expand breaks as needed
  scale_y_continuous(trans = compose_trans("log10", "reverse"), 
                     breaks = c(1,10^-2,10^-4,10^-6,10^-8,10^-10)) +
  
  #Add labels for genes of interest
  geom_text_repel(aes(label=ifelse(ID %in% genes_label, ID, NA)), 
                  size = 3, min.segment.length = 0, max.overlaps = k) #+
  
  #To add annotation arrows for direction and position of sensitizing and resistance genes,
  # Cut and paste from above, ensuring environment has variables arrow_y, sens_col, and res_col

  

#

## END ----




#