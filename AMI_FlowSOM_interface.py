# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 23:44:11 2022

@author: Lu Ri
"""

import pandas as pd
import numpy as np

som_data = pd.read_csv('D:/CAMP Biochemical/20220709 AMI_automatic_gating/gated_single_csv/all_data_2.5std.csv')
som_data['Condition'].unique().tolist() # to check the sequence
som_data = som_data.loc[som_data['Condition'].isin(conditions)].copy().iloc[:,:5].reset_index(drop = True)
# som_data.loc[som_data['Condition'].isin(conditions)]['Condition'].unique().tolist()
# f_data.loc[f_data['Condition'].isin(conditions)]['Condition'].unique().tolist() # to check the sequence
manual_gating_info = f_data.loc[f_data['Condition'].isin(conditions)].sort_values(by = ['Condition','Sample_NO']).iloc[:,5:].reset_index(drop = True)

som_result = pd.read_csv('D:/CAMP Biochemical/20220523 AMI_TOTAL_ANALYSIS/All_automated_figures/Cluster_of_8/AMI_all_AMI_hc_donor_no_force_zero_thresed..ExtNode.csv')
som_data['SOM'] = som_result['FlowSOM']//100-1
som_data = pd.concat([som_data,manual_gating_info],axis = 1)


som_grid = pd.read_csv('D:/CAMP Biochemical/20220523 AMI_TOTAL_ANALYSIS/All_automated_figures/Cluster_of_8/ClusterExplorer_Node_AMI_all_AMI_hc_donor_no_force_zero_thresed.fcs_runInfo.csv')
#sns.barplot(data = som_grid, x= 'Cluster',y = ' Size')
fig,ax = plt.subplots(figsize = (10,16),dpi = 200)
sns.heatmap(som_grid[['Cluster']+channels[1:]+['CD66b','CD3']].set_index('Cluster').T, square = True, linewidths=.5,ax = ax)
ax.set_xticklabels([i for i in range(8)],rotation = 0)
ax.yaxis.tick_right()
som_list = sorted(som_data['SOM'].unique().tolist())

#import matplotlib as mpl
cmap = plt.cm.gist_rainbow # define the colormap
palette = [cmap((i+1)*32-1)for i in range(len(som_list))]

ss_cnt = []
for cond, g in som_data.groupby(by = 'Condition'):
    s_cnt =[0 for _ in range(len(som_list))]
    for i, som in enumerate(som_list):
        s_cnt[i] = len(g.loc[g['SOM'] == som])
    ss_cnt.append(s_cnt+[cond])
    
    
ss_cnt = [ss_cnt[-1]]+ss_cnt[:-1] # move HC to front
fig,axes = plt.subplots(1,4,figsize = (10,3),dpi = 200)
for i, ax in enumerate(axes.flatten()):
    ax.pie(ss_cnt[i][:-1],radius = 1.2,startangle = 90, explode=[0.01 for _ in range(len(som_list))],colors = palette)

manual_gating_in_som = []
for n, g in som_data.groupby(by = 'SOM'):
    tcc = len(g)
    for i, phe in enumerate(new_col_list):
        per = len(g.loc[g['Phenotype'] == phe])/tcc
        manual_gating_in_som.append([str(n),per,i,phe])
        
manual_gating_in_som = pd.DataFrame(manual_gating_in_som,columns = ['SOM','per','Phe_NO','Phenotype'])
manual_gating_in_som[['SOM','per','Phe_NO']].pivot(index = ['SOM'],columns = 'Phe_NO', values = 'per').plot.bar(stacked = True)
    

patient_dir = 'D:/CAMP Biochemical/20220523 AMI_TOTAL_ANALYSIS/ALL_subtracted/'
patient = []
for p in glob.glob(patient_dir+'AMI001_T1*'):
    cond = os.path.basename(p).partition('_')[-1]
    cond = cond[:cond.find('c')]
    p = pd.read_csv(p,index_col=0)
    p['Condition'] = cond
    patient.append(p)
patient = pd.concat(patient,ignore_index = True)
patient.columns = f_data.columns.tolist()[:5]+['Condition']
#patient = patient.loc[patient['Sample_NO'] == patient_no].copy()
chans = channels[1:]
for ch in chans:
    patient = patient.loc[patient[ch]>patient[ch].quantile(0.01)].copy()
    patient = patient.loc[patient[ch]<patient[ch].quantile(0.99)].copy()
    # if patient[ch].min()<0:
    #     patient[ch] += abs(patient[ch].min())
    #patient = patient.loc[patient[ch]<patient[ch].quantile(0.9)].copy()
# for ch in chans:
#     patient[ch] = (0.00001+(patient[ch]-patient[ch].min())/(patient[ch].max()-patient[ch].min()))*100000
pal = dict(zip(['T1','T2','T3'],['tab:orange','tab:blue','tab:green']))
pp = sns.pairplot(data = patient[chans+['Condition']], hue = 'Condition', palette = pal, vars = chans,plot_kws = {'s': 10,'alpha':0.6})#,diag_kind="hist") corner=True, 


pal = dict(zip(['CD66b+CD3-','CD66b+CD3+','CD3+','CD66b-CD3-'],['tab:purple','tab:red','tab:green','gold']))
phe_data = patient[['CD66b','CD3']+['Condition']].copy()

def celltype(row):
    if row['CD3']>1000:
        if row['CD66b']>300:
            return 'CD66b+CD3+'
        else:
            return 'CD3+'
    elif row['CD66b']>300:
        return 'CD66b+CD3-'
    else:
        return 'CD66b-CD3-'
phe_data['Phenotype'] = phe_data.apply(lambda row: celltype(row),axis = 1)
pp = sns.pairplot(data = phe_data, hue = 'Phenotype', palette = pal, vars = ['CD66b','CD3'],plot_kws = {'s': 10,'alpha':0.6},diag_kind="hist",diag_kws = {'stat':'count'})#,diag_kind="hist") corner=True, diag_kws = {'common_norm':False}
# for ax in pp.axes.flat:
#     if ax is not None:
#         ax.set(xscale="log")
#         ax.set(yscale="log")



############## UMAP display ######################

umap = pd.read_csv('D:/CAMP Biochemical/20220523 AMI_TOTAL_ANALYSIS/All_automated_figures/FlowSOM_analysis/10-Jul-2022/UMAP_ID_GYYF/Results_UMAP_GYYF.csv')
x,y = umap.columns.tolist()
umap = pd.concat([umap,manual_gating_info],axis = 1)
umap['FlowSOM'] = som_result['FlowSOM']//100-1
umap['FlowSOM'] = umap['FlowSOM'].astype(str)
fig,axes = plt.subplots(2,5, figsize=(20,20))
fig.tight_layout()
axes = axes.flatten()
for i, hue in enumerate(['FlowSOM']+manual_gating_info.columns.tolist()[:5]+['Condition']):
    sns.scatterplot(data =umap, x = x,y = y,hue = hue,s = 5, alpha = 0.7,ax =axes[i])

sns.scatterplot(data =umap, x = x,y = y,hue = 'FlowSOM',s = 5, alpha = 0.7)