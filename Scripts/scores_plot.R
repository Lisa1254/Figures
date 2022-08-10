# Scores Plot

#'
#'This script will make a plot to compare two kinds of scores for a gene/residue.
#'Customization ?in progress
#'This script was drafted for comparing functional scores and structural deltadeltaG scores for amino acid residues in a protein under study. Can be applied to any relevant comparison of two numeric value scales for specific points, for example, comparing the p-values of two methods of significance testing for a set of genes.
#'

## Libraries, Functions, & Parameters ----

#Note: plyr is also used as plyr::mapvalues, and should be installed, but does not need to be loaded to the environment. If you do load plyr, it should be loaded before dplyr
library(dplyr) #for manipulating dataframes
library(stringr) #for extracting string information (in functions)
library(zoo) #for rolling mean used in plot
library(ggplot2) #for plot
library(ggrepel) #for scatterplot labels

#Functions
source("Functions/parse_mutStrings.R")

#This file imports three functions:
#parse_p_hgvs() and split_mutation_string() both take a dataframe that contains a mutation string, and returns the same dataframe, but with 3 new columns:
# ir = initial residue
# pos = position of residue
# fr = final residue, i.e. the amino acid after mutation
#parse_p_hgvs() expects the mutation string to be in the nomenclature of a protein by HGVS, while split_mutation_string() expects a specific number of characters for each initial and mutated amino acids
#mut2self.pro() takes as input a dataframe with the initial amino acid and final after mutation, and outputs a vector with the row index of all scores that represent mutations to the same amino acid as initially expressed. The function recognizes "=" as final residue indicating mutation to self. AA can be expressed as 1 or 3 letters.

#Express start and end position of residues for plotting
start_pos <- 1
end_pos <- 100


#
## Process Data ----

#Read in data
#Data here is taken as a subsection of real protein data used in the lab, but with the amino acids changed
f_scores <- read.delim("Example_Data_Files/example_protein_functionalScores.txt")
ddg_scores <- read.delim("Example_Data_Files/example_protein_ddGScores.txt")
#Note, add following arguments if necessary: header = FALSE, col.names = c("hgvs_pro", "score")

#Separate initial residue, position, final residue from mutation code. 
#The functional score data is in hgvs nomenclature
f_scores <- parse_p_hgvs(f_scores, 1)
#The ddg scores as saved from a mutatex style output of foldX are formatted in hgvs nomenclature with 1-letter codes, without the "p." at the beginning, but this will work still
ddg_scores <- parse_p_hgvs(ddg_scores, 1)

#When using the ddg scores from the Position Scan command, there can be digits in the 3 letter amino acid scores (e.g. H1S and H2S for protonated states of histidine), and the string will include the chain. If your mutation strings are in this style (e.g. METB1K), then you can use this function
#ddg_scores <- split_mutation_string(ddg_scores, 1, 3, 1, chain=TRUE)

#Get score mean for each position
#Functional scores will have mutate to "=" to indicate mutation to self, and "Ter" to indicate termination. For this analysis we'll remove both.
f_scores_plot <- f_scores %>%
  filter(fr != "=") %>%
  filter(fr != "Ter") %>%
  group_by(pos) %>%
  summarise(avg=mean(score)) %>%
  arrange(pos)

#When using the ddg scores from the Position Scan command, the final residue for histidine can be expressed as "o" or "e" for its protonated version, and these might not be have been recognized and removed previously for mutation to self. Can use the following to replace "o" and "e" with
#ddg_scores[which(ddg_scores$fr %in% c("o", "e")), "fr"] <- "H"
#Then use mut2self.pro function to identify which rows have mutation to self and remove
#ddg_selfMut <- mut2self.pro(ddg_scores, 3, 5)
#ddg_scores <- ddg_scores[-ddg_selfMut, ]

#ddg scores as calculated through mutatex wrapper of foldX's BuildModel command do not include any mutations to self, and do not have any Terminations marked
ddg_scores_plot <- ddg_scores %>%
  group_by(pos) %>%
  summarise(avg=mean(score)) %>%
  arrange(pos)

#Check that all positions are accounted for
#Functional Scores
seq(start_pos,end_pos)[-which(seq(start_pos,end_pos) %in% unique(f_scores_plot$pos))]
#ddG Scores
seq(start_pos,end_pos)[-which(seq(start_pos,end_pos) %in% unique(ddg_scores_plot$pos))]

#If something is missing, can can add NA values, or impute as desired
# e.g. if the last residue (in this case 100) is missing, can do:
# f_scores_plot <- rbind(f_scores_plot, c(100,NA))


#

## Plot Using Rolling Mean ----

#Define window
window = 5
space_padding = (window-1)/2
#Get rolling mean values
rm_f_score <- rollmean(f_scores_plot$avg, k = window)
rm_ddg_score <- rollmean(ddg_scores_plot$avg, k = window)
#Set up as dataframe
#Use space padding to centre values on same amino acids as RM was calculated from 
ggdf.RM <- data.frame(num=seq(start_pos,end_pos), 
                      F_Score=c(rep(NA, space_padding), rm_f_score, rep(NA, space_padding)), 
                      ddG_Score=c(rep(NA, space_padding), rm_ddg_score, rep(NA, space_padding)))

#Define y-axis scales
# in this example, f_score is primary axis
ylim.prim <- c(-1*max(ggdf.RM$F_Score, na.rm=T), -1*min(ggdf.RM$F_Score, na.rm=T))   # in this example, f_score
ylim.sec <- c(min(ggdf.RM$ddG_Score, na.rm=T), max(ggdf.RM$ddG_Score, na.rm=T))

#Set up scaling parameters for secondary y-axis
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

#Plot
ggplot(ggdf.RM, mapping = aes(x=num)) +
  geom_point(mapping = aes(y=-1*F_Score), color = "black") +
  geom_line(mapping = aes(y=-1*F_Score), color = "black") +
  geom_point(mapping = aes(y=a+b*ddG_Score), color = "red", shape = 15) +
  geom_line(mapping = aes(y=a+b*ddG_Score), color = "red") +
  scale_y_continuous("Functional Score", 
                     breaks = c(-1, -0.8, -0.6, -0.4),
                     labels = c(1, 0.8, 0.6, 0.4),
                     sec.axis=sec_axis(~(.-a)/b, name="ddG Score")) +
  #Here replace with ~./scaleFactor
  theme(axis.line.y.right = element_line(color = "red"), 
        axis.ticks.y.right = element_line(color = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.title.y.right = element_text(color = "red")) +
  scale_x_continuous("AA Position")

#
## Plot using Scores ----

#Very similar to above with rolling mean. Using the scores instead of rolling mean for this type of data is less informative as a visualization. Although there is more data being presented, the smoothing of the rolling mean provides an easier to interpret image

#Collect into single dataframe here
ggdf.Scores <- data.frame(num=seq(start_pos, end_pos),
                          F_Score=f_scores_plot$avg,
                          ddG_Score=ddg_scores_plot$avg)

#Define y-axis scales
# in this example, f_score is primary axis
ylim.prim <- c(-1*max(ggdf.Scores$F_Score, na.rm=T), -1*min(ggdf.Scores$F_Score, na.rm=T))   
ylim.sec <- c(min(ggdf.Scores$ddG_Score, na.rm=T), max(ggdf.Scores$ddG_Score, na.rm=T))

#Set up scaling parameters for secondary y-axis
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

#Plot
ggplot(ggdf.Scores, mapping = aes(x=num)) +
  geom_point(mapping = aes(y=-1*F_Score), color = "black") +
  geom_line(mapping = aes(y=-1*F_Score), color = "black") +
  geom_point(mapping = aes(y=a+b*ddG_Score), color = "red", shape = 15) +
  geom_line(mapping = aes(y=a+b*ddG_Score), color = "red") +
  scale_y_continuous("Functional Score", 
                     breaks = c(-1, -0.8, -0.6, -0.4),
                     labels = c(1, 0.8, 0.6, 0.4),
                     sec.axis=sec_axis(~(.-a)/b, name="ddG Score")) +
  theme(axis.line.y.right = element_line(color = "red"), 
        axis.ticks.y.right = element_line(color = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.title.y.right = element_text(color = "red")) +
  scale_x_continuous("AA Position")



#
## Scatterplot Correlation ----

#Combine the ddg_scores and f_scores into a single dataframe using preferred mutation code format for labels.
#Here I'll be using 3-letter amino acid code, position, one letter amino acid code.
aa1 <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
aa3<- c("ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR")

#Map 3 letter code in "fr" column of functional scores to 1 letter code
map_fr_fun <- plyr::mapvalues(toupper(f_scores$fr), aa3, aa1)
#Map 1 letter code in "ir" column of ddg scores to 3 letter code
map_ir_ddg <- plyr::mapvalues(ddg_scores$ir, aa1, aa3)

f_scores$mutation_code <- paste0(toupper(f_scores$ir), f_scores$pos, map_fr_fun)
ddg_scores$mutation_code <- paste0(map_ir_ddg, ddg_scores$pos, ddg_scores$fr)

#Merge scores by mutation code, only retains mutations with data in both inputs
merge_scores <- merge(f_scores[,c("mutation_code", "score")], ddg_scores[,c("mutation_code", "score")], by = "mutation_code")
colnames(merge_scores)[c(2,3)] <- c("f_score", "ddg_score")

#From here I saved the dataframe as txt file to explore in my plot labelling shiny so that it was easy to identify which mutations met different thresholds, or to locate the scores of a specific mutation.
#As an example of what that exploration might produce, a scatterplot is demonstrated below

#Preview shape of plot:
ggplot(merge_scores, mapping = aes(x=ddg_score, y=f_score)) +
  geom_point() +
  scale_y_continuous(trans = "reverse")

#Set minimum ddG score for labelling: 
ddg_min <- 30
#Set maximum Functional score for labelling:
fun_max <- 0

#Choosing to label if a point meets either threshold, rather than both
merge_scores$label <- ifelse((merge_scores$ddg_score >= ddg_min) | (merge_scores$f_score <= fun_max),"Yes", "No")

#Also in preview noticed an outlier at f_score > 3, adding that to the labelling set
merge_scores[which(merge_scores$f_score > 3), "label"] <- "Yes"

#In progress, updating this figure
ggplot(merge_scores, mapping = aes(x=ddg_score, y=f_score)) +
  geom_point(aes(color=label)) +
  scale_color_manual(values = c("gray","steelblue")) +
  scale_y_continuous(trans = "reverse") +
  geom_text_repel(aes(label = ifelse(label == "Yes", mutation_code, "")), size = 2.5, max.overlaps = 25, min.segment.length = 0) +
  theme(legend.position = "none") +
  labs(x="ddG Score", y="Functional Score")


