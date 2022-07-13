# Figures
Repository for code for common figures made in R for CRISPR analysis
  
## Contents

### Scripts
1. Volcano Plot (using LFC from Mageck Scores)  
    a. Highlighting significant genes  
    b. Highlighting specific genes  
2. Radar Plot (using Zscore rank from DrugZ output)  
    a. Ranking of single gene  
    b. Ranking of multiple genes  
    c. Ranking of gene with background colour highlight  
    d. Use of custom scaling function with 1 or more genes  
3. Heatmap (using binary table)  
    a. Basic Overlap Plot (gene/GO)  
    b. Basic Overlap Plot (gene/Screen) aka FINGERPRINT PLOT  
  
### Functions
1. compose_trans: copy of source code for compose_trans from scales package, which I couldn't get to update properly  
2. radarchart_ca: modification of radarchart from fmsb to allow for custom axis and grid background highlighting to better demonstrate biological significance  
  
## Currently Working on
* Adding script for a figure that compares two types of scores at a given series of points, where each score type has a different scale.  
  
## To Do
* Add annotation sections for heatmap  
* In radar script, update customization at top  
* In radarchart_ca, see if I can switch the order of numerical axis input to circumference (small number) first  
  