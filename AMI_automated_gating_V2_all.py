# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 14:08:10 2022

@author: Lu Ri
"""

import glob
import os
import pandas as pd
import numpy as np
from numpy import sqrt, log10
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.collections import EllipseCollection
import re # to parse the string by specific char
import flowkit as fk

main_dir = 'D:/CAMP Biochemical/20220309 Baseline paired strict and without FSC/Orig/'
##### step 1 to read all the BG droplet and plot them as bar, to see what to do next ######
# bg_all = []
# for d in glob.glob(main_dir+'*/'):
#     spl_name = d[:-1].rpartition('\\')[-1]
#     bg_path = glob.glob(d+'*_BaselineP2P.pkl')[0]
#     bg = pd.read_pickle(bg_path)
#     channels = bg.columns.tolist()
#     top = bg.quantile(0.9)
#     bottom = bg.quantile(0.1)
#     bool_map = bg.apply(lambda row: (row>bottom)&(row<top),axis = 1)
#     bg_means = bg[bool_map.eq(1).all(axis = 1)].mean().to_frame().reset_index()
#     # wide_forms
#     # bg_means = bg[bool_map.eq(1).all(axis = 1)].mean().to_frame().T
#     bg_means['Sample_name'] = spl_name 
#     bg_all.append(bg_means)
# bg_all = pd.concat(bg_all, ignore_index = True)
# bg_all.columns = ['channel','v','Sample_name']

# fig, axes = plt.subplots(5,1,figsize = (15,8),sharex = True)
# for i,ax in enumerate(axes.flatten()):
#     sns.barplot(data = bg_all.loc[bg_all['channel']==channels[i]], x= 'Sample_name', y = 'v',ax = axes[i])
#     ax.tick_params(axis='x',rotation = 90)

#### autogating, but still need to do with the comped data and bg, because only when you are using comp subtracted, you know for sure what channel is positive and what is negative

def do_comp(df_orig,interp = False,  matrix_path = 'D:/CAMP Biochemical/flowjo matrix/faye_for_python_channel_name.csv'):
    df = df_orig.copy()
    if interp:
        df = df.ffill().add(df.bfill()).div(2)
    sample = fk.Sample(df)
    detectors = [sample.pnn_labels[i] for i in sample.fluoro_indices]
    comp_mat = fk.Matrix(
        'my_spill',
        matrix_path,
        detectors
    )
    sample.apply_compensation(comp_mat)
    comp = sample.as_dataframe(source='comp',col_names = detectors) 
    comp.columns = channels
    return comp    

def get_gating_thres(bg_comp, factor = 2.5): # factor is to divide the distance between mean and 90th percentile, compare 4 vs 5, 5 i think ok so 4 will be ok , now we pick the one that gives best t-test result
    bg = bg_comp.copy()
    top = bg.quantile(0.9)
    bottom = bg.quantile(0.1)
    bool_map = bg.apply(lambda row: (row>bottom)&(row<top),axis = 1)
    #bg_means = bg[bool_map.eq(1).all(axis = 1)].mean()
    bg_thres =  factor*bg[bool_map.eq(1).all(axis = 1)].std()
    #bg_thres = top + (top-bg_means)/factor
    #bg_thres = (top-bg_means)/factor
    return bg_thres

def do_gating(df_comp, bg_comp_thres_shifted, channels):
    df = df_comp.copy()
    new_col_list = [ch+'_gate' for ch in channels]
    gated_result = df > bg_comp_thres_shifted
    gated_result.columns = new_col_list
    gated_result = gated_result.astype(int)
    def categorize(row):
        cat = '|'
        for ch in new_col_list:
            if row[ch] ==1:
                cat+=ch.rpartition('_')[0]+'|'
        return cat
    gated_result['Phenotype'] = gated_result.apply(lambda row: categorize(row), axis = 1)
    new_col_list.append('Phenotype')
    return gated_result # only return gating part

def force_zero(comp_subtracted_withgating,channels):
    def set_zero(row):
        for ch in channels:
            if row[ch+'_gate']==0:
                row[ch] = 0
        return row
    df = comp_subtracted_withgating.copy()
    df = df.apply(lambda row: set_zero(row),axis = 1)
    return df

def shift_negatives(comp_subtracted_withgating,channels, shift):
    def set_zero(row):
        for ch in channels:
            if row[ch+'_gate']==0:
                row[ch] -= shift
        return row
    df = comp_subtracted_withgating.copy()
    df = df.apply(lambda row: set_zero(row),axis = 1)
    return df


def parse_name(spl_name, spl_no):
    if 'C0' in spl_name:
        cond = 'HC'
        no = spl_no[cond]
        spl_no[cond]+=1
    elif 'AMI' in spl_name:
        cond = spl_name[:spl_name.find('0')] + ' ' + spl_name.rpartition('_')[-1]
        no = spl_no[cond]
        spl_no[cond]+=1
    elif 'AHF' in spl_name:
        cond = spl_name[:spl_name.find('0')] + ' ' + spl_name.rpartition('_')[-1]
        no = spl_no[cond]
        spl_no[cond]+=1
    elif 'ARREST' in spl_name:
        cond = 'ARREST'
        no = spl_no[cond]
        spl_no[cond]+=1
    return dict({'Condition':cond,'Sample_NO':no}),spl_no
#####################################################################


channels = ['NE','GzB','CD66b','CD3','CD31']
cnt_100k = pd.read_csv('D:/CAMP Biochemical/20220709 AMI_automatic_gating/count_100k.csv')

spl_no = dict({'HC':1,
               'AMI T1': 1,'AMI T2': 1,'AMI T3': 1,
               'AHF T1':1, 'AHF T2':1, 'AHF T3':1, 
               'ARREST':1 })
def get_data(spl_no):
    data_all = []
    bg_thres_all = []
    for d in glob.glob(main_dir+'*/'):
        spl_name = d[:-1].rpartition('\\')[-1]
        #print(spl_name)
        if 'AMI007' not in spl_name:
            if 'C003' not in spl_name:
                if 'C004' not in spl_name:
                    meta,spl_no = parse_name(spl_name, spl_no)
                    #limit = cnt_100k[(cnt_100k['Condition'] == meta['Condition']) & (cnt_100k['Sample_NO'] == meta['Sample_NO']) ]['Cnt_100k'].values[0]
                    limit = 3000
                    bg_path = glob.glob(d+'*_BaselineP2P.pkl')[0]
                    bg = pd.read_pickle(bg_path).iloc[:limit,:].copy().reset_index(drop = True)
                    bg.columns = channels
                    
                    data_path = glob.glob(d+'*_P2P.pkl')[0]
                    data = pd.read_pickle(data_path).iloc[:limit,:].copy()
                    data.columns = channels
                    data_comp = do_comp(data)
                    
                    bg_comp = do_comp(bg,interp = True)
                    bg_thres = get_gating_thres(bg_comp) # small distance
                    b = bg_thres.to_frame().T
                    b['Sample_name'] = spl_name
                    bg_thres_all.append(b)
                    
                    
                    gating_result = do_gating(data_comp,(bg_comp+bg_thres),channels)
                    comp_subtracted = data_comp.reset_index(drop = True)-bg_comp.reset_index(drop = True)
                    
                    #comp_subtracted+=500
                    comp_subtracted = pd.concat([comp_subtracted,gating_result],axis = 1)
                    comp_subtracted = shift_negatives(comp_subtracted, channels,1500)
    
                    
                    for key in meta:
                        comp_subtracted[key] = meta[key]
                    #comp_subtracted['Sample_name'] = spl_name
                    #comp_subtracted[channels].to_csv(export_dir1 + spl_name+'comp_subtracted.csv',index = False)
                    data_all.append(comp_subtracted)
    
    data_all = pd.concat(data_all,ignore_index = True)
    return data_all

f_data = get_data(spl_no)

#data_all.to_csv(export_dir2 + 'all_data_2.5std.csv',index = False)

#bg_thres_all = pd.concat(bg_thres_all,ignore_index = True)