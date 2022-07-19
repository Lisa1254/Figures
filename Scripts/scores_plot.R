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

#If using rolling window for the mean, define window
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

#
#This hasn't been edited yet, and is just pasted from my prev work to generalize for use here
#

#Scatterplot of all variants ddG & functional score

library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)

# Import functional scores
f_scores <- read.csv("GNB1L_functional scores.csv")

#Separate initial residue, position, final residue from mutation code
#This uses strsplit and strextract commands and the positional digits
#In this data chart, all mutation codes begin with "p." so that needs to be removed
ir_f <- unlist(lapply(f_scores[,1], function(x){gsub("p.", "", strsplit(x, "[[:digit:]]+")[[1]][1])}))
fr_f <- unlist(lapply(f_scores[,1], function(x){strsplit(x, "[[:digit:]]+")[[1]][2]}))
pos_f <- unlist(lapply(f_scores[,1], function(x){str_extract(x, "[[:digit:]]+")}))

#Add separated information to dataframe
f_scores <- cbind(f_scores, ir=ir_f, pos=as.numeric(pos_f), fr=fr_f)

# Import ddG scores
ddg_scores <- read.delim("GNB1L_foldX_output_add227.txt", col.names=c("Mutation", "score"))

split_mutation_string <- function(df, col, nchar1, nchar2, chain=FALSE) {
  ir <- substr(df[,col], 1, nchar1)
  fr <- substr(df[,col], nchar(df[,col])- nchar2 + 1, nchar(df[,col]))
  if (chain == TRUE) {
    pos <- as.numeric(substr(df[,col], nchar1+2, nchar(df[,col])-nchar2))
  } else {
    pos <- as.numeric(substr(df[,col], nchar1+1, nchar(df[,col])-nchar2))
  }
  df <- cbind(df, ir, pos, fr)
  return(df)
}

#Use on ddg table
ddg_scores <- split_mutation_string(ddg_scores, 1, 3, 1, chain=TRUE)

#Remove self comparisons in odd H residues (copied from moving window script, but may not be necessary here)
rows2del <- ddg_scores %>%
  filter(ir %in% c("HIS", "H1S", "H2S")) %>%
  filter(fr %in% c("H", "e", "o")) %>%
  pull(Mutation)

ddg_scores <- ddg_scores[-which(ddg_scores$Mutation %in% rows2del),]

#Get plotting data

#Make common mutation code for 2 types of input

aa1 <- c('G', 'A', 'L', 'V', 'I', 'P', 'R', 'T', 'S', 'C', 'M', 'K', 'E', 'Q', 'D', 'N', 'W', 'Y', 'F', 'H')
aa3<- c('GLY', 'ALA', 'LEU', 'VAL', 'ILE', 'PRO', 'ARG', 'THR', 'SER', 'CYS', 'MET', 'LYS', 'GLU', 'GLN', 'ASP', 'ASN', 'TRP', 'TYR', 'PHE', 'HIS')

map_fr_fun <- mapvalues(toupper(f_scores$fr), aa3, aa1)

f_scores$mutation_code <- paste0(toupper(f_scores$ir), f_scores$pos, map_fr_fun)
ddg_scores$mutation_code <- paste0(ddg_scores$ir, ddg_scores$pos, ddg_scores$fr)

merge_scores <- merge(f_scores[,c("mutation_code", "score")], ddg_scores[,c("mutation_code", "score")], by = "mutation_code")
colnames(merge_scores)[c(2,3)] <- c("f_score", "ddg_score")

ggplot(merge_scores, mapping = aes(x=f_score, y=ddg_score)) +
  geom_point()

write.table(merge_scores, file = "Rolling_window_analysis/commonScores_func_ddg.txt", sep = "\t")

#Using instead
f_scores2 <- read.csv("Rolling_window_analysis/GNB1L_GFP_scores.csv", comment.char = "#")
#Modifications done as above, but in different window

map_fr2_fun <- mapvalues(toupper(f_scores2$fr), aa3, aa1)
f_scores2$mutation_code <- paste0(toupper(f_scores2$ir), f_scores2$pos, map_fr2_fun)

merge_scores2 <- merge(f_scores2[,c("mutation_code", "score")], ddg_scores[,c("mutation_code", "score")], by = "mutation_code")
colnames(merge_scores2)[c(2,3)] <- c("f_score2", "ddg_score")

ggplot(merge_scores2, mapping = aes(x=ddg_score, y=f_score2)) +
  geom_point() +
  scale_y_continuous(trans = "reverse")

write.table(merge_scores2, file = "Rolling_window_analysis/commonScores_funcGFP_ddg.txt", sep = "\t")


