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
3. Heatmap (using binary table)  
    a. Basic Overlap Plot (gene/GO)  
    b. Basic Overlap Plot (gene/Screen) aka FINGERPRINT PLOT  
  
### Functions
1. compose_trans: copy of source code for compose_trans from scales package, which I couldn't get to update properly  
2. radarchart_ca: modification of radarchart from fmsb to allow for custom axis and grid background highlighting to better demonstrate biological significance  
  
## Currently Working on
radarchart_ca custom axis function  
* Want the input of numerical axis to have circumference number first, not centre  
* Want to adjust drawing of grid segments to include fill for outer segments  
* Can remove all the extra annotation styles, since this will be custom only (if no input value labels, use as.character for numerical input)
  
## To Do
* Add annotation sections for heatmap  
* In radar script, update customization at top  
  