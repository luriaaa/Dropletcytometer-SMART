# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 18:12:39 2022

@author: Lu Ri
"""
import pandas as pd
import numpy as np
from numpy import sqrt, log10, log, log2, ceil
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, ttest_rel
from matplotlib.ticker import MaxNLocator,PercentFormatter



def bar_all_features(df, mypalette, x = 'Condition', y = 'v', conds = ['HC','AMI T1','AHF T1','ARREST'], plotname = 'unnamed'):
    
    #get the xtab for each, if x,y test significant display, run the script in a sequence   
    def unit_conversion(data,y):
        df = data.copy()
        ch = df['ch'][0]
        if 'CD' in ch:
            k = unit_maps[ch]['k']
            b = unit_maps[ch]['b']
            df[y] = df[y].apply(lambda x: 10**((log10(max(x,0.01))-b)/k))
            df[y] *= 11.5/1000000 # convert to pg
        elif ch == 'cnt':
            pass
        else:
            k = unit_maps[ch]['k']
            b = unit_maps[ch]['b']
            df[y] = (df[y]-b)/k
            df[y] *= 11.5/1000000 # convert to pg
        return df
    
    
    def do_barplot(x, y,df, ax):
    
        def set_limit(v):
            if v <=100:
                return max(1.2,v)
            else:
                factor = max(2, len(str(int(ceil(v))))-2)
                return ceil(v/(10**factor))*10**factor

        f,p = do_test_unpaired(df,x,y)
        print(p)

        plt_df = unit_conversion(df,y)    # the exported data is already converted....
        # plt_df = df.copy()
        
        if plt_df[y].min()<0:
            plt_df[y] -= plt_df[y].min() # in case there are negative values
    
        # variable useful for annotating
        #xticks = plt_df[x].unique().tolist() # should be sorted
        bar_plot2 = sns.barplot(data=plt_df, x=x, y=y, ax = ax, alpha = 0.8,errwidth = 2,palette=mypalette)
        ylims = list(ax.get_ylim())
        
        def mutate_bar_width(bar, newwidth):
            x = bar.get_x()
            width = bar.get_width()
            centre = x + width/2.
            bar.set_x(centre - newwidth/2.)
            bar.set_width(newwidth)
            return 0
            
        mutate_bar_width(bar_plot2.patches[0], 0.5)
        #mutate_bar_width(bar_plot2.patches[len(cond)], 0.5)
       
        # sns.stripplot(data=plt_df, x=x, y=y, color = 'dimgray', size = 5,alpha = 0.8, ax = ax)
    
    
        # ax.set_ylim([0,100])
        # ax.yaxis.set_major_formatter(PercentFormatter())
        ax.yaxis.set_major_locator(MaxNLocator(5))   
         # limit y major tick to 5
        
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.xaxis.set_ticklabels([])

        return p,plt_df
    
    ############################################# t-test 
    def do_test_unpaired(df_to_test,group_by,val):
        df_to_test = df.copy()
        conds = df[group_by].unique().tolist()
        s = len(conds)
        f =  np.zeros((s,s))
        p = np.zeros((s,s))
        
        for i, c1 in enumerate(conds):
            for j, c2 in enumerate(conds):
                if i<j:
                    f[i,j] = log2(df.loc[df[group_by]==c2][val].mean()/df.loc[df[group_by]==c1][val].mean())

                    p[i,j] = ttest_ind(df.loc[df[group_by]==c1][val],df.loc[df[group_by]==c2][val], equal_var=True, nan_policy='omit')[1]
        def to_df(m):
            df = pd.DataFrame(m,columns = conds.copy(),index = conds.copy())
            return df.iloc[:-1,:].copy()
        
        return to_df(f),to_df(p)
    
    def do_test_paired(df_to_test,group_by,val):
        df = df.copy()
        conds = df[group_by].unique().tolist()
        s = len(conds)
        f =  np.zeros((s,s))
        p = np.zeros((s,s))
        
        for i, c1 in enumerate(conds):
            for j, c2 in enumerate(conds):
                if i<j:
                    f[i,j] = log2(df.loc[df[group_by]==c2][val].mean()/df.loc[df[group_by]==c1][val].mean())
                    p[i,j] = ttest_rel(df.loc[df[group_by]==c1][val],df.loc[df[group_by]==c2][val], nan_policy='omit')[1]

        def to_df(m):
            df = pd.DataFrame(m,columns = conds.copy(),index = conds.copy())
            return df.iloc[:-1,:].copy()
        
        return to_df(f),to_df(p)
    
    
    ############################################# annotate method
    def annotate_plot(ax,p, xs, ylims, ast_offset = 0.08, linespacing = 0.3):
        def p_marking_rule(p):
            def pmap(ind_p):
                if ind_p == 0:
                    return ''
                elif ind_p >0.05:
                    return ''
                elif ind_p > 0.01:
                    return '*'
                elif ind_p > 0.001:
                    return '**'
                else:
                    return '***'
            return p.applymap(pmap)
    
        p_mapped = p_marking_rule(p)
        y_int = ylims[1]-ylims[0]
        
        def plot_line_ast(ax,start,end,text, yline,yast,color):
            ax.plot([start,end],[yline,yline],linewidth = 1.5,color =color)
            ax.text(start+0.5*(end-start),yast,text, ha="center", va="center",fontsize=18,color = color)
            return 0
        
        y_marker = ylims[0]     
        # upper lines
        for x1 in ax.get_xticks()[1:]:
            if p[xs[x1]][xs[0]] <0.05:
                y_line = y_marker +linespacing*y_int*2
                y_star = y_line+ast_offset*y_int*2
                plot_line_ast(ax,0,x1,p_mapped[xs[x1]][xs[0]],y_line,y_star,'black')
                y_marker = y_line+linespacing*y_int*2
        
        # even_upper_lines
        for x1 in ax.get_xticks()[1:]:
            for x2 in ax.get_xticks()[2:]:
                if x2>x1:
                    if p[xs[x2]][xs[x1]] <0.05:
                        y_line = y_marker+linespacing*y_int*2
                        y_star = y_line+ast_offset*y_int*2
                        plot_line_ast(ax,x1,x2,p_mapped[xs[x2]][xs[x1]],y_line,y_star,'tab:green') #
                        y_marker = y_line+linespacing*y_int*2
        ax.set_ylim(ylims[0],y_marker)
        # your ylims will get extended again
        return 0

    order = dict(zip(conds,[i for i in range(len(conds))]))
    df = df.loc[df[x].isin(conds)].sort_values(by=[x], key=lambda x: x.map(order)).reset_index(drop = True)    
    #fig,ax = plt.subplots(1,1, figsize = (5,8), dpi = 200)
        
    fig, (ax_ann, ax) = plt.subplots(2, 1, figsize = (4.8,8), dpi = 200, gridspec_kw={'height_ratios': [1, 5]},sharex = True)
    fig.tight_layout()
    
    p,plt_df = do_barplot(x,y,df.copy(),ax)
    annotate_plot(ax_ann, p, conds, list(ax.get_ylim()))       
    #finalize adjusting the proplot
    ax.xaxis.set_tick_params(which='minor',bottom=False)
    ax.grid(False)
    
    ax_ann.grid(False)
    ax_ann.axis('off')

    return fig, ax