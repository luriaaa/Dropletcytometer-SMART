# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 13:29:58 2022

@author: Lu Ri
"""
from scipy.stats import pearsonr 

### 1. to scan for the featuers correlative to time for PCI (ami t1)
tpci = pd.read_excel('D:/CAMP Biochemical/20220523 AMI_TOTAL_ANALYSIS/All_automated_figures/clinical_params/Time_to_PCI.xlsx')
days = pd.read_excel('D:/CAMP Biochemical/20220523 AMI_TOTAL_ANALYSIS/All_automated_figures/clinical_params/Total_hospitalization_days.xlsx')
t1 = all_features.loc[all_features['Condition']=='AMI T1'].copy().sort_values(by = 'Sample_NO').reset_index(drop = True)
t2 = all_features.loc[all_features['Condition']=='AMI T2'].copy().sort_values(by = 'Sample_NO').reset_index(drop = True)
t3 = all_features.loc[all_features['Condition']=='AMI T3'].copy().sort_values(by = 'Sample_NO').reset_index(drop = True)
all_t = all_features.loc[all_features['Condition'].str.contains('AMI')].copy().sort_values(by = ['Condition','Sample_NO']).reset_index(drop = True)

import math
def fold_change(s_norm, s_dnorm):
    if s_dnorm == 0:
        fc = math.copysign(1 , s_norm - s_dnorm)
    elif s_norm == 0:
        fc = math.copysign(1, s_norm - s_dnorm)
    else:
        fc = np.log2(s_norm/s_dnorm)
    return fc

feat_list = []
for col in t1.columns.tolist()[2:]:
    r,p = pearsonr(tpci['Time to PCI'].values, t1[col].values)
    if (abs(r)>0.6) and (p<0.05):
        feat_list.append([col,r,p])
        
#feat_list = pd.DataFrame(feat_list,columns = ['Feature_name','r','p'])

for num, feat in enumerate(feat_list):
    fig,(ax,ax2) = plt.subplots(1,2, figsize = (10,5),dpi = 200)
    df = pd.DataFrame([tpci['Time to PCI'],t1[feat[0]], t1['Sample_NO'].astype(str)]).T
    sns.scatterplot(data = df, x = df.columns[0] , y= df.columns[1],hue = 'Sample_NO', ax = ax)
    ax.set_title(feat[0] + ' r=' + "{:.2f}".format(feat[1]) + ' p= '+ "{:.3f}".format(feat[2]))
    ax.set_xlabel('Time to PCI [min]')
    ax.get_legend().remove()
    
    
        
    # t1_val = t1[feat[0]].values
    # t2_val = t2[feat[0]].values
    # t3_val = t3[feat[0]].values
    
    post_df = pd.concat([all_t[feat[0]].to_frame(),days[days.columns.tolist()[:3]]],axis = 1)
    post_df = post_df.loc[~post_df['Condition'].str.contains('T1')].copy()
    post_df['Sample_NO'] = post_df['Sample_NO'].astype(str)
    sns.lineplot(data = post_df, x = 'Time[Days]' , y= feat[0],hue = 'Sample_NO', ax = ax2)
    
    # fc12 = []
    # fc13 = []
    
    # for i in range(len(t1_val)):
    #     fc12.append(fold_change(t2_val[i],t1_val[i]))
    #     fc13.append(fold_change(t3_val[i],t1_val[i]))
    # fc_df = pd.concat([t1['Sample_NO'].astype(str).to_frame(), pd.DataFrame(fc12,columns = ['T2vsT1']),pd.DataFrame(fc13,columns = ['T3vsT1'])],axis = 1)
    # fc_df = pd.melt(fc_df,id_vars = 'Sample_NO',var_name='Timepoints',value_name = 'Fold-change')
    # sns.stripplot(data = fc_df, x = 'Timepoints' , y= 'Fold-change',hue = 'Sample_NO', ax = ax2)
    
    
    # Shrink current axis by 20%
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    
    fig.savefig(export_dir+str(num)+'tpci.png',bbox_inches = 'tight')
    plt.close()
    
    

### 2. to scan for the features between AMI T1 and T2 that sample_no = 6 is different from rest of sample, use absolute fold change