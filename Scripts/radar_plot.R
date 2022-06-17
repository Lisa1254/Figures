# Radar Plot

#'
#'This script will make a series of radar plots of Zscore rank from DrugZ output.
#'Customization for 
#'

## Libraries & Functions ----

library(fmsb) #Has radarchart function
library(scales) #Allows transparency added to colors

#Modification on radarchart source code to better suit the biological application
source("Functions/radarchart_ca.R")

#

## Data ----

#Import example data of output from DrugZ analysis
# NOTE: This is a real output from an analysis, but randomized to preserve data sensitivity, then subset for first 1500 genes to keep file size manageable. The rankings were adjusted 1:1500 for each screen in order of previous ranking

gene_ranks <- read.delim("Example_Data_Files/example_gene_ranking_1500.txt")

#Name columns by screen name
# NOTE: fmsb plots in column order, counter-clockwise starting at the top.
order_screens <- c("NT_RKO", "NT_COL", "AraC", "HU", "Aph")
colnames(gene_ranks) <- order_screens


#Choose genes of interest
#Chose 5 randomly from data for interest of plotting
genes <- c("AACS", "AAMP", "AHNAK2", "ARMCX3", "ANGPT2")

#Set up dataframe of rankings, with first row = "Max" or the circumference of the plot, and second row = "Min" or the centre of the plot
# NOTE: here since the highest rank value is the minimum importance, we are putting the highest rank as the "Max" value.
plot_data <- t(data.frame(Circumference = rep(1,5), Centre = rep(max(gene_ranks),5)))
colnames(plot_data) <- order_screens

#For this to run, be sure that the order of screens is the same in both the gene_ranks and plot_data dataframes
for  (g in genes) {
  plot_data <- rbind(plot_data,
                     gene_ranks[which(rownames(gene_ranks) == g), ])
}
rownames(plot_data)[-c(1,2)] <- genes

#If using fmsb radarchart function, best to use a log scale
plot_data_log <- log10(plot_data)


#

## Plot Single Gene (Log Values, using fmsb radarchart function) ----

#Identify which gene will be plotted (index by row number of plot_data_log, so should be 3 or greater)
gene_ind <- 3

#Since the plotting data is log-transformed, I'd still like to label the figure with the actual rank. So here a vector of the labels is computed.
max_rank_log <- max(plot_data_log)
axis_labels <- c(max(plot_data), 
                 round(10^(0.75*max_rank_log)), 
                 round(10^(0.5*max_rank_log)),
                 round(10^(0.25*max_rank_log)), 
                 "1")

radarchart(plot_data_log[c(1,2,gene_ind),],
           title = rownames(plot_data_log)[gene_ind],
           
           # Variable options (Vertices of radar)
           #vlabels = c("New", "names", "for", "each", "vertex"),
           vlcex=0.8,     #controls the font size of variable labels
           
           # Polygon options
           #pty = 16,    #A vector to specify point symbol: Default 16 (closed circle),
           pcol="darkgreen",    #line color
           pfcol=scales::alpha("darkgreen",0.5), #scales::alpha allows transparency value
           #plwd=4,    #line types. Can be a numeric vector 1:6 or a character vector 
                        #c(“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”).
                        #To remove the line, use plty = 0 or plty = “blank”.
           plty=1,    #line width
           
           # Grid Options
           cglcol="navy",      #line color
           #cglty=1,           #line type (see above, default = dashed)
           cglwd=1,            # line width
           
           # Axis options
           #centerzero = T,     #Logical. If true, this function draws charts with scaling
                                #originated from (0,0). If false, charts originated
                                #from (1/segments).
           axislabcol="blue",    #color of axis label and numbers. Default is “blue”
           caxislabels = axis_labels, #defined above
           calcex=0.8,       #Font size magnification for caxislabels.
           axistype = 1 #The type of axes, specified by any of 0:5. 0 means no axis label. 
                        #1 means center axis label only. 2 means around-the-chart label only.
                        #3 means both center and around-the-chart (peripheral) labels.
                        #4 is *.** format of 1, 5 is *.** format of 3. Default is 0.
           
            
           )

#

## Plot multiple genes (Log Values, using fmsb radarchart function) ----

#Using above data and arguments, add more rows to the df input of the plot function. It is not recommended to plot more than 2-3 data types on a single plot because it gets confusing

gene_ind1 <- 4
gene_ind2 <- 5

radarchart(plot_data_log[c(1,2,gene_ind1, gene_ind2),],
           title = paste0("Comparison of ", rownames(plot_data_log)[gene_ind1], " and ", rownames(plot_data_log)[gene_ind2]),
           
           # Variable options (Vertices of radar)
           vlcex=0.8,
           
           # Polygon options
           pty = c(15,16),    
           pcol=c("darkgreen", "darkgoldenrod"),
           pfcol=c(scales::alpha("darkgreen",0.5), NA), 
           plwd=4,
           plty=1,
           
           # Grid Options
           cglcol="navy",      #line color
           #cglty=1,           #line type (see above, default = dashed)
           cglwd=1,            # line width
           
           # Axis options
           axislabcol="blue",    #color of axis label and numbers. Default is “blue”
           caxislabels = axis_labels, #defined above
           calcex=0.8,       #Font size magnification for caxislabels.
           axistype = 1 
           
           
)


#
## Plot with high rank highlight (Log Values, using fmsb radarchart function) ----

#Add row to plot for highlight data
#Choosing for HL the inner 4 sections with gray, outermost section left white. Taking formula from above for calculating axis labels
plot_data_h <- rbind(plot_data[1:2,], 
                     HL=rep(round(10^(0.25*max_rank_log)),5), 
                     plot_data[3:nrow(plot_data),])
plot_data_h_log <- log10(plot_data_h)

#Define gene index (now 4 or greater from the additional row)
#Could use for loop with 1:length(genes) + 3
gene_ind <- 4

radarchart(plot_data_h_log[c(1:3,gene_ind),],
           title = rownames(plot_data_h_log)[gene_ind],
           
           # Variable options (Vertices of radar)
           vlcex=0.8,
           
           # Polygon options
           pty = c(NA,16),    
           pcol=c(NA, "darkgreen"),
           pfcol=c(scales::alpha("gray78", 0.7),scales::alpha("darkgreen",0.5)), 
           plwd=4,
           plty=1,
           
           # Grid Options
           cglcol="navy",      #line color
           #cglty=1,           #line type (see above, default = dashed)
           cglwd=1,            # line width
           
           # Axis options
           axislabcol="blue",    #color of axis label and numbers. Default is “blue”
           caxislabels = axis_labels, #defined above
           calcex=0.8,       #Font size magnification for caxislabels.
           axistype = 1 
           
           
)

#

## Single Gene using Custom Axis fuction ----

radarchart_ca(plot_data[3,],
              axistype=1,
              num_ax=c(1500,400,100,10,1),
              title = rownames(plot_data)[3],
              pcol="darkgreen" , 
              pfcol=scales::alpha("lightgreen", 0.3), 
              plty=1,
              cglcol="navy",
              axislabcol="blue",
              cglwd=1, vlcex=0.8, calcex=0.8
)

#

## Plot multiple genes using Custom Axis function ----

radarchart_ca(plot_data[3:4,],
              axistype=1, 
              num_ax=c(1500,400,100,10,1),
              title = paste0(rownames(plot_data)[3], " & ", rownames(plot_data)[4]),
              pcol=c("darkorange4", "darkblue"), 
              pfcol=c(scales::alpha("darkorange", 0.3), scales::alpha("blue", 0.3)), 
              plty=1,
              cglcol="navy",
              axislabcol="blue",
              cglwd=1, vlcex=0.8, calcex=0.8
)

legend(x = "bottomleft",
       legend = c(rownames(plot_data)[3], rownames(plot_data)[4]), # Vector with the name of each group
       #fill = c(scales::alpha("darkorange", 0.3), scales::alpha("blue", 0.3)),   # Creates boxes in the legend with the specified colors
       col = c("darkorange4", "darkblue"),
       lty='solid',
       pch=19,
       title = "Gene")


#
