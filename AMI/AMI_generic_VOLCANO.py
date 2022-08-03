# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:16:41 2022

@author: Lu Ri
"""

import numpy as np
from scipy.stats import ttest_ind, ttest_rel
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# function for AMI vs HC (non-paired, the calculation of t and fold is different)
def volcano_data(cond1,cond2, phe_list, paired = False):
    # cond1 is "control" condition, use as denominator   
    s = []
    features = ['cnt','NE','GzB','CD31']
    data = f_data.loc[f_data['Condition'].isin([cond1,cond2])].copy()
    for phe_info in phe_list:
        phe = phe_info[0]
        phe_num = phe_info[1]
        df_all = data.loc[data['Phenotype']==phe].copy()
        for feat in features:
            if feat == 'cnt':
                df = df_all.groupby(by = ['Condition','Sample_NO'])[['Phenotype']].count()
            else:
                df = df_all.groupby(by = ['Condition','Sample_NO'])[[feat]].median()
                if df[df.columns[0]].max() <1:
                    continue
            df.columns = ['v']
            df = df.reset_index()
            #df.columns = ['Condition','Sample_NO','v']
            df = df.sort_values(by = ['Condition','Sample_NO'])

            if paired:
                # sort first because want toensure the paired t test is done with correct pair
                t = ttest_rel(df.loc[df['Condition']==cond1]['v'], df.loc[df['Condition']==cond2]['v'],nan_policy='omit')[1]
                # still using the fold change of all pop (more robust), instead of mean of each pair (very heterogenuous)
                f =np.nan if df.loc[df['Condition']==cond1]['v'].mean() ==0 else df.loc[df['Condition']==cond2]['v'].mean()/df.loc[df['Condition']==cond1]['v'].mean()     
            else:
                # unpaired
                t = ttest_ind(df.loc[df['Condition']==cond1]['v'], df.loc[df['Condition']==cond2]['v'],equal_var=True,nan_policy='omit')[1]
                f = np.nan if df.loc[df['Condition']==cond1]['v'].mean() ==0 else df.loc[df['Condition']==cond2]['v'].mean()/df.loc[df['Condition']==cond1]['v'].mean()                                                                                                         
            
            df_out = pd.DataFrame([[phe_num,phe,feat, t,f]],columns = ['Phe_NO','Phenotype','ch','p-val','fold'])
            s.append(df_out)
    return pd.concat(s,ignore_index = True)



def volcano_plot(s, mypalette):
    df = s.dropna()
    x_shift = 0.01
    y_shift = 0.01    
    df['p-val'] = -np.log10(df['p-val'])
    df['fold'] = np.log2(df['fold'])
    df_to_annotate = df.loc[df['p-val']>-np.log10(0.05)].copy()
    df_to_annotate = df_to_annotate.loc[(df_to_annotate['fold'].abs()>=1)|(df_to_annotate['p-val']>2)].copy()
    
    # df_to_highlight = df_to_annotate.copy()
    c1 = df_to_annotate.loc[df_to_annotate['Phe_NO'] == 0].copy()
    c2 =df_to_annotate.loc[df_to_annotate['Phe_NO'] == 16].copy()
    c3 = df_to_annotate.loc[df_to_annotate['Phe_NO'] == 24].copy()
    # df_to_highlight['sum_of_xy'] = abs(df_to_highlight['p-val'])+abs(df_to_highlight['fold'])
    # df_to_highlight = df_to_highlight.sort_values(by = ['sum_of_xy'],ascending = False).iloc[:3,:].copy()
    
    
    
    fig, ax = plt.subplots(1,1,sharex = True, sharey = True, figsize = (5,5),dpi = 200)
    dot_size = 20
    alpha = 0.8
    sns.scatterplot(ax = ax, data =df, x = 'fold',y='p-val',hue = 'ch',  palette = mypalette, s = dot_size, alpha = alpha)


    for i in range(len(df_to_annotate)):
        ax.text((df_to_annotate.iloc[i,:]['fold']+x_shift),(df_to_annotate.iloc[i,:]['p-val']+y_shift),df_to_annotate.iloc[i,:]['Phe_NO'], fontsize=12)
    
    for color, c in zip(['gold','darkgrey','peru'],[c1,c3,c2]):
        for i in range(len(c)):
            ax.arrow((c.iloc[i,:]['fold']+0.3),(c.iloc[i,:]['p-val']+0.3),
                     -0.05,-0.05,head_width = 0.1,
                     color = color)
    
    # for i in range(len(df_to_highlight)):
    #     ax.arrow((df_to_highlight.iloc[i,:]['fold']+0.3),(df_to_highlight.iloc[i,:]['p-val']+0.3),
    #              -0.05,-0.05,head_width = 0.1,
    #              color = 'firebrick')
    
    ymax =4 
    ax.axhline(y = -np.log10(0.05),ls = '--', lw = 0.5,c = 'k')
    ax.axhline(y = -np.log10(0.01),xmin = 1/3,xmax = 2/3,ls = '--', lw = 0.5)
    ax.axvline(x = 1,ymax = 2/ymax, ls = '--', lw = 0.5)
    ax.axvline(x = -1,ymax = 2/ymax,ls = '--', lw = 0.5)
    ax.set_xlabel('log2(fold-change)')
    ax.set_ylabel('-log10(p-value)')
    ax.legend(loc = 'upper left')
    ax.set_xlim([-3,3])
    ax.set_ylim([0,ymax])

    fig.tight_layout()
 
    return fig, [c1,c2,c3]