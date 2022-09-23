#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 09:59:18 2022

@author: lhoeg
"""

# %% Libraries and data

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#import matplotlib.gridspec as gridspec
plt.rcParams['font.size']=14

f_scores = pd.read_csv(r"~/Durocher Lab Dropbox/lisahoeg@gmail.com/GitRepositories/Figures/Example_Data_Files/example_protein_functionalScores.txt", sep="\t")
ddg_scores = pd.read_csv(r"~/Durocher Lab Dropbox/lisahoeg@gmail.com/GitRepositories/Figures/Example_Data_Files/example_protein_ddGScores.txt", sep="\t")
c_scores = pd.read_csv(r"~/Durocher Lab Dropbox/lisahoeg@gmail.com/GitRepositories/Figures/Example_Data_Files/example_protein_conservationScores.txt", sep="\t")

#Delete these after script is done. Just for testing with original data
ddg_fun_scores = pd.read_csv(r"~/Durocher Lab Dropbox/lisahoeg@gmail.com/Projects/5_YZ_Corr_220629/Rolling_window_analysis/avg_scores_ddgM_functional.txt", sep='\t')
con_scores = pd.read_csv(r"~/Durocher Lab Dropbox/lisahoeg@gmail.com/Projects/5_YZ_Corr_220629/Rolling_window_analysis/scores_conservation.txt", sep='\t')


# %% Process Data

#This works for f_scores since those have the "p." at the front, but not for the ddg scores which don't. Maybe ad an if statement for re.search? Can import re to help. Or maybe I can just do something with import re afterall.
def parse_p_hgvs(df,col):
    ir = df.iloc[:,col].str.split("[0-9]").str.get(0).str.split(".").str.get(1)
    fr = df.iloc[:,col].str.split("[0-9]").str.get(1)
    pos = df.iloc[:,col].str.extract(r'(\d+)')
    out_df = pd.concat([df, ir, pos, fr], axis=1, ignore_index=True)
    return out_df

f_scores = parse_p_hgvs(f_scores,0)
ddg_scores2 = parse_p_hgvs(ddg_scores, 0)

ir_ddg = ddg_scores.iloc[:,0].str.split("[0-9]").str.get(0)
ir_f = f_scores.iloc[:,0].str.split("[0-9]").str.get(0).str.extract(r'([A-Z]+)')

#df_scores = ddg_fun_scores.merge(con_scores, on="num")

#rm5_scores = df_scores.rolling(5, center = True).mean()
