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
4. Scores (2 datasets with same X-axis, different Y-axis)  
    a. Plot Using Rolling Mean  
    b. Plot using Scores  
    c. Scatterplot Correlation  
  
### Functions
1. compose_trans: copy of source code for compose_trans from scales package, which I couldn't get to update properly  
2. radarchart_ca: modification of radarchart from fmsb to allow for custom axis and grid background highlighting to better demonstrate biological significance  
3. parse_mutStrings.R: script contains three functions for working with common mutation expressions:  
    a. parse_p_hgvs: uses HGVS nomenclature as expected input format to parse initial residue, position, final residue into dataframe columns  
    b. split_mutation_string: if input mutation expression does not follow HGVS nomenclature, but is consistent for using 1/3 letter abbreviations at initial and final residue positions  
    c. mut2self.pro: detects rows in dataframe with mutation to self from initial and final residue, including using "=" and if 1/3 letter codes are used  
  
## Currently Working on
  
  
## To Do
* Add annotation sections for heatmap  
* In radar script, update customization at top  
* In radarchart_ca, see if I can switch the order of numerical axis input to circumference (small number) first  
  