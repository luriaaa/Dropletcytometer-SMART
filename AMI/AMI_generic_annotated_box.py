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
from matplotlib.patches import FancyBboxPatch
from matplotlib.path import get_path_collection_extents


def box_all_features(df, mypalette, x = 'Condition', y = 'v', conds = ['HC','AMI T1','AHF T1','ARREST'], plotname = 'unnamed'):
    
    
    #get the xtab for each, if x,y test significant display, run the script in a sequence   
    def unit_conversion(data,y):
        df = data.copy()
        ch = df['ch'][0]
        if 'CD' in ch:
            k = unit_maps[ch]['k']
            b = unit_maps[ch]['b']
            df[y] = df[y].apply(lambda x: 10**((log10(max(x,0.01))-b)/k))
            df[y] *= 11.5/1000 # convert to fg
        elif ch == 'cnt':
            pass
        else:
            k = unit_maps[ch]['k']
            b = unit_maps[ch]['b']
            df[y] = (df[y]-b)/k
            df[y] *= 11.5/1000 # convert to fg
        return df
    
    def do_boxplot(x, y,df, ax):
        
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
        
        # here, you can define the min max of the df
        medianprops = dict(linestyle='', linewidth=4, color='black', alpha = 1)
        meanpointprops = dict(marker='', markeredgecolor='black',
                              markerfacecolor='black')
        means = plt_df.groupby(by= [x])[y].mean().sort_index(key=lambda x: x.map(order))
        
        
        box_plot = sns.boxplot(data=plt_df, x=x, y=y,  showcaps = False, showmeans = True, #showfliers = False,
                          medianprops = medianprops,meanprops = meanpointprops,palette=mypalette, ax = ax)
        
        ylims = list(ax.get_ylim())
        new_patches = []
        for patch in reversed(ax.artists):
            bb = patch.get_path().get_extents()
            color = patch.get_facecolor()
            p_bbox = FancyBboxPatch((bb.xmin+0.15*bb.width, bb.ymin),
                                    abs(bb.width*0.7), abs(bb.height),
                                    boxstyle="round,pad=0,rounding_size=0.2",
                                    ec="black", fc=color,linewidth = 2,alpha=0.8,
                                    mutation_aspect=0.03*set_limit(ylims[1]))
            #print(abs(bb.width*0.7))
            patch.remove()
            new_patches.append(p_bbox)
        for patch in new_patches:
            ax.add_patch(patch)    
        # sns.stripplot(data=plt_df, x=x, y=y, color = 'dimgray', size = 5,alpha = 0.8, ax = ax)
        
        vertical_offset = 0
        
        for xtick in box_plot.get_xticks():
            box_plot.text(xtick,means[xtick] + vertical_offset,f'{means[xtick]:.1f}', 
                    horizontalalignment='center',size='x-small',color='w',fontsize = 10,weight='semibold',bbox=dict(facecolor='black',ec = None))
        
        ax.yaxis.set_major_locator(MaxNLocator(5))   
        
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
    
    p,plt_df = do_boxplot(x,y,df.copy(),ax)
    annotate_plot(ax_ann, p, conds, list(ax.get_ylim()))       
    #finalize adjusting the proplot
    ax.xaxis.set_tick_params(which='minor',bottom=False)
    ax.grid(False)
    
    ax_ann.grid(False)
    ax_ann.axis('off')

    return fig, ax