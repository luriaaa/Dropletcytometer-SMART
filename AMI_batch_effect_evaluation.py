# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 22:19:17 2022

@author: Lu Ri
"""


import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import ptitprince as pt


f_data = []
for k in data_all.keys():
    if (k == 'AMI') or (k =='AHF'):
        for key in data_all[k].keys():
            frame = data_all[k][key].copy()
            frame['Condition'] = k + ' ' + key
            f_data.append(frame)
    else:
        frame = data_all[k].copy()
        frame['Condition'] = k
        f_data.append(frame)
    
f_data = pd.concat(f_data, ignore_index = True)

channels = f_data.columns.tolist()[:5]



bat = []
bat_u = []
for s in samples:
    bat_u.append(f_data_unsub.loc[(f_data_unsub['Condition'] == s[0])& (f_data_unsub['Sample_NO'] == s[1])].copy())
    bat.append(f_data.loc[(f_data['Condition'] == s[0])& (f_data['Sample_NO'] == s[1])].copy())

bat = pd.concat(bat,ignore_index=True)
bat_u = pd.concat(bat_u,ignore_index=True)

bat['experiment'] = bat['Condition']+'_'+bat['Sample_NO'].astype(str)
bat_u['experiment'] = bat_u['Condition']+'_'+bat_u['Sample_NO'].astype(str)

for ch in channels:
    fig, ax = plt.subplots(figsize=(5, 6))
    fig.tight_layout()
    ax=pt.RainCloud(x = 'experiment', y = ch, data = bat_u, bw =  .2, box_showfliers = False,
                     width_viol = .6, ax = ax, orient = 'h')
    ax.set_xlim([-1000,15000])
    fig.savefig(export_dir+'batch_evalu'+ch+'.png',bbox_inches = '')