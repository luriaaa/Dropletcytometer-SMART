# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 20:50:43 2022

make summary bigger, compiling everything, 

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


def do_gating(df, gatings, channels):
    new_col_list = []
    for ch in channels:
        new_col = ch+"_gate"
        new_col_list.append(new_col)
        thres = gatings[ch]
        fn = lambda x: 0 if x[ch]<thres else 1
        df[new_col]=df.apply(fn, axis=1) 
    def categorize(row):
        cat = '|'
        for ch in new_col_list:
            if row[ch] ==1:
                cat+=ch.rpartition('_')[0]+'|'
        return cat
    df['Phenotype'] = df.apply(lambda row: categorize(row), axis = 1)
    return 0

def convert_value(df_gated,unit_maps):
    cols = unit_maps.columns.tolist()
    for i,col in enumerate(cols):
        if i<2:
            df_gated[col] = (df_gated[col] - unit_maps[col]['b'])/unit_maps[col]['k']
        else:
            df_gated[col] = 10**((log10(df_gated[col]) - unit_maps[col]['b'])/unit_maps[col]['k'])
    return 0

def prepare_data(convert = False, limit = False):
    sample_no = {'AMI':{'T1':0,'T2':0,'T3':0},'AHF':{'T1':0,'T2':0,'T3':0},'HC':0,'ARREST':0}
    for p in paths:
        s_name = os.path.basename(p)[:os.path.basename(p).find('comp')] 
        if s_name in gatings.columns:
            data = pd.read_csv(p, index_col=0)
            #data = filter_data(data)
            #data[data<0] = 0 # enable if you want to set negative values to 0
            
            channels = [col.partition('_')[0] for col in data.columns.to_list()]
            data.columns = channels
            gate = gatings[s_name]
            #data['Phenotype'] = data.apply(lambda row: Neutrophil_CD31_gating(row),axis = 1)
            #pad_zeros(data)
            do_gating(data,gate,channels)
            c_limit = len(data)
            if convert:
                convert_value(data,unit_maps)
            if 'AMI' in s_name:
                if '007' not in s_name:
                    tp = s_name.partition('_')[-1]
                    sample_no['AMI'][tp]+=1
                    data['Sample_NO'] = sample_no['AMI'][tp]
                    if limit:
                        c_limit = cnt_100k.query("Condition == 'AMI' and Sample_NO == %s and Timepoint == '%s'" % (sample_no['AMI'][tp],tp))['Cnt_100k'].values[0]
                    data_all['AMI'][tp].append(data.iloc[:c_limit,:].copy())
            elif 'AHF' in s_name:
                tp = s_name.partition('_')[-1]
                sample_no['AHF'][tp]+=1
                data['Sample_NO'] = sample_no['AHF'][tp]
                if limit:
                    c_limit = cnt_100k.query("Condition == 'AHF' and Sample_NO == %s and Timepoint == '%s'" % (sample_no['AHF'][tp],tp))['Cnt_100k'].values[0]
                data_all['AHF'][tp].append(data.iloc[:c_limit,:].copy())
            elif'C' in s_name:
                if '003' not in s_name and '004' not in s_name:
                    sample_no['HC']+=1
                    data['Sample_NO'] = sample_no['HC']
                    if limit:
                        c_limit = cnt_100k.query("Condition == 'HC' and Sample_NO == %s and Timepoint == 'T1'" % (sample_no['HC']))['Cnt_100k'].values[0]
                    data_all['HC'].append(data.iloc[:c_limit,:].copy())
            elif 'ARREST' in s_name:
                sample_no['ARREST']+=1
                data['Sample_NO'] = sample_no['ARREST']
                if limit:
                    c_limit = cnt_100k.query("Condition == 'ARREST' and Sample_NO == %s and Timepoint == 'T1'" % (sample_no['ARREST']))['Cnt_100k'].values[0]
                data_all['ARREST'].append(data.iloc[:c_limit,:].copy())

    condition = 'AMI'
    for key in data_all[condition].keys():
        data_all[condition][key] = pd.concat(data_all[condition][key], axis = 0, ignore_index = True)
        
    condition = 'AHF'
    for key in data_all[condition].keys():
        data_all[condition][key] = pd.concat(data_all[condition][key], axis = 0, ignore_index = True)  
        
    condition = 'ARREST'
    data_all[condition] = pd.concat(data_all[condition], axis = 0, ignore_index = True)  
    condition = 'HC'
    data_all[condition] = pd.concat(data_all[condition], axis = 0, ignore_index = True)  
    return 0

# need to generate custom RGB color bar value, and custom color bar, as a dictionary
def stack_data(force_zero = True):
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
    def set_zeros(row):
        for ch in channels:
            if row[ch+'_gate'] == 0:
                row[ch] = 0
        return row
    if force_zero:
        f_data = f_data.apply(lambda row: set_zeros(row),axis = 1)
    return f_data

# designed specific for 100k data, because we sum the counts directly
def summarize_data_100k(method = 'mean'):
    channels = f_data.columns.tolist()[:5]
    phenotypes=sorted(f_data['Phenotype'].unique())
    conditions = sorted(f_data['Condition'].unique())
    
    ############################################### ver 1 exclusive phenotypee ########################################

    counts = []
    features = []
    def fill_data(df,prefix,basename):
        df.columns = ['v']
        df['Feature_name'] = '['+prefix+']'+basename
        return df.reset_index()
    
    def long_to_wide(s):
        wide = s.pivot(index = ['Condition','Sample_NO'],columns = 'Feature_name', values = 'v').fillna(0).reset_index()   
        return wide.loc[:, ~(wide == 0).all()].copy()
    
    def wide_to_long(s):
        return pd.melt(s, id_vars=['Condition','Sample_NO'], value_vars=[cat for cat in s.columns.tolist() if (cat!='Condition')&(cat!='Condition')],
                var_name='Feature_name',value_name = 'v')
    
    for phe in phenotypes:
        f_data_p = f_data.loc[f_data['Phenotype']==phe].copy()
        count_df = pd.DataFrame(f_data_p.groupby(by = ['Condition','Sample_NO']).count().iloc[:,0])

        counts.append(fill_data(count_df,'cnt',phe))
        
        for ch in channels:
            features.append(fill_data(f_data_p[[ch,'Condition','Sample_NO']].groupby(by = ['Condition','Sample_NO']).agg(method),ch,phe))
    
        
    type1_data = []
    for d in [counts,features]:
        type1_data.append(pd.concat(d, ignore_index = True))
    type1_data = pd.concat(type1_data,ignore_index = True)
    type1_data = wide_to_long(long_to_wide(type1_data))
    
    
    ################################################ ver 3 2-level features ############################################
    def do_counting(df):
        cnt = pd.DataFrame([df.iloc[:,0].count()], columns = ['v'])
        return cnt
    
    def do_percentage(df, df_main):
        return 0
    
    def get_value(df,keyword):
        if len(keyword) == 0:
            return pd.DataFrame([0], columns = ['v'])
        else:            
            return pd.DataFrame([df[keyword].agg(method)], columns = ['v'])
            
    def cd31_count(df):
        df_len = len(df)
        df1 = df.loc[df['CD31_gate']==1].copy()
        pos_cnt = pd.DataFrame([df1.iloc[:,0].count()], columns = ['v'])
        neg_cnt = pd.DataFrame([df_len-df1.iloc[:,0].count()], columns = ['v'])
        return pos_cnt, neg_cnt
    
    def cd31_cd31(df):
        df1 = df.loc[df['CD31_gate']==1].copy()
        return pd.DataFrame([df1['CD31'].agg(method)], columns = ['v'])
    
    def cd31_var(df,keyword):
        df1 = df.loc[df['CD31_gate']==1].copy()
        df2 = df.loc[df['CD31_gate']==0].copy()
        if len(keyword) == 0:
            return pd.DataFrame([0], columns = ['v']),pd.DataFrame([0], columns = ['v'])
        else:            
            return pd.DataFrame([df1[keyword].agg(method)], columns = ['v']),pd.DataFrame([df2[keyword].agg(method)], columns = ['v'])
   
    def associate_info(df,name,index,level):
        df = df.copy()
        df['Feature_name'] = name
       # df['Feature_lvl'] = level
        df['Condition'] = index[0]
        df['Sample_NO'] = index[1]
        return df
    
    # lastly compile them as a dataframe
    cnt = []
    val = []
    cd31_pos_cnt = [] # cd31 neg count will be generated later, as [pos_cnt,neg_cnt], so we first show enzyme vs phenotype, then show enzyme vs heart marker, 
    cd31_neg_cnt = [] # neg count is convenient to check the CD31- ne/gzb level
    cd31_pos_cd31_val =[] # cd31 var to light up the positive cd31, when plotting the pie, you can force the color of the negative to black, no need extra dataframe for neg_cd31
    cd31_pos_val = [] # ne or gzb
    cd31_neg_val = [] # ne or gzb, for cd31 negative cell
    
    def collect_data(f,n,kw,f_name,fl):
        # f = df at that level, n = index, kw = kw_level, f_name , fl = feature level
        cnt.append(associate_info(do_counting(f),f_name,n,fl))
        val.append(associate_info(get_value(f,kw),f_name,n,fl))
        
        pos, neg = cd31_count(f)
        cd31_pos_cnt.append(associate_info(pos,f_name,n,fl))
        cd31_neg_cnt.append(associate_info(neg,f_name,n,fl))
        
        cd31_pos_cd31_val.append(associate_info(cd31_cd31(f),f_name,n,fl))
        
        pvar, nvar = cd31_var(f,kw)
        cd31_pos_val.append(associate_info(pvar,f_name,n,fl))
        cd31_neg_val.append(associate_info(nvar,f_name,n,fl))
        return 0

    # we can examine what to plot for the  cd31 map later
    for n, g in f_data.groupby(by = ['Condition', 'Sample_NO']):
        # level 1 counts and features
        ne = g.loc[g['NE_gate']==1].copy()
        gzb = g.loc[g['GzB_gate']==1].copy()
        double_neg = g.loc[(g['GzB_gate']==0) & (g['NE_gate']==0)].copy()
        
        f_lv1 = [ne,gzb,double_neg]
        f_name_lv1 = ['NE+','GzB+','NE-GzB-']
        kw_lv1 = ['NE','GzB','']
        
        for i, f in enumerate(f_lv1):
            collect_data(f,n,kw_lv1[i],f_name_lv1[i],1)
            # level 2 counts and features
            cd66b = f.loc[f['CD66b_gate']==1].copy()
            # cd66bcd3 =  f.loc[(f['CD66b_gate']==1)&(f['CD3_gate']==1)].copy()
            cd3 = f.loc[f['CD3_gate']==1]
            cd3 = cd3.loc[cd3['CD66b_gate']==0].copy()
            rest_of_cell = f.loc[(f['CD66b_gate']==0) & (f['CD3_gate']==0)].copy()
            f_lv2 = [cd66b,cd3,rest_of_cell]
            f_name_lv2 = [f_name_lv1[i]+ name for name in ['CD66b+','CD3+','CD66b-CD3-']]
            for j, f2 in enumerate(f_lv2):
                collect_data(f2,n,kw_lv1[i],f_name_lv2[j],2)


    def mutate_feature(l_df, prefix):
        def set_prefix(row):
            if 'auto' not in prefix:
                return '['+prefix+']'+row['Feature_name']
            else:
                index = row['Feature_name'][:6].find('+')
                if index == -1:
                    return 'none' # later we can remove those "rubbish features" by eliminating feature_name =='none'
                else:
                    return '['+prefix[:prefix.find('auto')]+row['Feature_name'][:index]+']'+row['Feature_name']
        l_df = pd.concat(l_df,ignore_index = True)
        l_df['Feature_name'] = l_df.apply(lambda row: set_prefix(row),axis = 1)            
        return l_df.loc[~l_df['Feature_name'].str.contains('none')].copy()
    
    #data_list = [cnt,val,cd31_pos_cnt,cd31_neg_cnt,cd31_pos_cd31_val,cd31_pos_val,cd31_neg_val]    
    # do some post processing to integrate data and change the feature_name
    # if i remember correctly, CD31 neg value(enzymes) dont make too much significance, so we can omit them
    type2_data = []
    for p, df in zip(['cnt','auto','CD31+cnt','CD31-cnt', 'CD31+auto','CD31-auto','CD31','CD31+ cells'],[cnt,val,cd31_pos_cnt,cd31_neg_cnt,cd31_pos_val,cd31_neg_val, cd31_pos_cd31_val,cd31_pos_val]):
        if p == 'CD31+ cells':
            f = mutate_feature(df,'')
            f['Feature_name'] = f.apply(lambda row: '('+p+')'+row['Feature_name'],axis = 1)
            type2_data.append(f)
            
        else:
            type2_data.append(mutate_feature(df,p))
    type2_data = pd.concat(type2_data,ignore_index = True)
    type2_data = wide_to_long(long_to_wide(type2_data))
    
    # goal is to compile type2_data into a single df, by mutating the feature name, you can see the cd31_pos_cd31_valk and cd31_pos_val are literally 2 channel of the same category
    
        
    # type 3 data is the CD66b+ CD3+ CD66b-CD3- data, followed by the CD31+ count, 
    type3_data = []
    cnt = []
    feats = []
    feat_cnt = []
        

    cd3 = f_data.loc[(f_data['CD3_gate']==1)&(f_data['CD66b_gate']==0)].copy()
    cd66b = f_data.loc[(f_data['CD66b_gate']==1)].copy()
    cd66bonly = f_data.loc[(f_data['CD3_gate']==0)&(f_data['CD66b_gate']==1)].copy()
    cd66bcd3 = f_data.loc[(f_data['CD3_gate']==1)&(f_data['CD66b_gate']==1)].copy()
    double_neg = f_data.loc[(f_data['CD3_gate']==0) & (f_data['CD66b_gate']==0)].copy()
    for ph, df in zip(['CD3+','CD66b+','CD66b+CD3-','CD66b+CD3+','CD66b-CD3-'],[cd3,cd66b,cd66bonly,cd66bcd3,double_neg]):
        cnt.append(fill_data(pd.DataFrame(df.groupby(by = ['Condition','Sample_NO']).count().iloc[:,0]),'cnt',ph))
        for ch in ['NE','GzB','CD31']:
            feats.append(fill_data(df[[ch,'Condition','Sample_NO']].groupby(by = ['Condition','Sample_NO']).agg(method),ch,ph))
            feat_cnt.append(fill_data(df[[ch+'_gate','Condition','Sample_NO']].groupby(by = ['Condition','Sample_NO']).agg('sum'),ch+'+cnt',ph))
    
    
    for d in [cnt,feats,feat_cnt]:
        type3_data.append(pd.concat(d, ignore_index = True))
    type3_data = pd.concat(type3_data,ignore_index = True)
    type3_data = wide_to_long(long_to_wide(type3_data))
    
    cd31_data = []
    cnt = []
    cnt_neg = []
    cd31_v = []
    for n, g in f_data.groupby(by = ['Condition', 'Sample_NO']):
        cd31 = g.loc[g['CD31_gate'] == 1].copy()
        cd31_neg = g.loc[g['CD31_gate'] == 0].copy()
        cnt.append(fill_data(pd.DataFrame(cd31.groupby(by = ['Condition','Sample_NO']).count().iloc[:,0]),'cnt','CD31+'))
        cnt_neg.append(fill_data(pd.DataFrame(cd31_neg.groupby(by = ['Condition','Sample_NO']).count().iloc[:,0]),'cnt','CD31-'))
        cd31_v.append(fill_data(pd.DataFrame(cd31[['CD31','Condition','Sample_NO']].groupby(by = ['Condition','Sample_NO']).agg(method)),'CD31','CD31+'))
    for d in [cnt,cd31_v,cnt_neg]:
        cd31_data.append(pd.concat(d, ignore_index = True))
    cd31_data = pd.concat(cd31_data,ignore_index = True)
    cd31_data = wide_to_long(long_to_wide(cd31_data))
    
    cd66b_cd3_data = []
    per_cd3_in_cd66b = []
    per_negzb = []  
    ###### these data can be retrieved by the type1_data dataframe
    # double_pos_cd31_cnt = []
    # double_pos_cd31neg_cnt = []
    # double_pos_cd31_ne = []
    # double_pos_cd31neg_ne = []
    
    # we can examine what to plot for the  cd31 map later
    def negzb_per(df):
        norm = len(df.loc[(df['NE_gate']==1)&(df['GzB_gate']==1)])
        denorm = len(df)
        return norm/denorm

    cd3 = f_data.loc[(f_data['CD3_gate']==1)&(f_data['CD66b_gate']==0)].copy()
    cd66b = f_data.loc[(f_data['CD66b_gate']==1)].copy()
    cd66bonly = f_data.loc[(f_data['CD66b_gate']==1)&(f_data['CD3_gate']==0)].copy()
    cd66bcd3 = f_data.loc[(f_data['CD3_gate']==1)&(f_data['CD66b_gate']==1)].copy()
    double_neg = f_data.loc[(f_data['CD3_gate']==0) & (f_data['CD66b_gate']==0)].copy()
    
    per = cd66b.groupby(by = ['Condition','Sample_NO'])['CD3_gate'].sum()/cd66b.groupby(by = ['Condition','Sample_NO'])['CD66b_gate'].sum()
    per_cd3_in_cd66b.append(fill_data(pd.DataFrame(per),'per','CD3+(CD66b+)'))
    for phe, df in zip(['CD3+','CD66b+CD3-','CD66b+CD3+','CD66b-CD3-'],[cd3,cd66bonly,cd66bcd3,double_neg]):
        per = df.groupby(by = ['Condition','Sample_NO']).apply(lambda df: negzb_per(df))
        per_negzb.append(fill_data(pd.DataFrame(per),'per','NE+GzB+'+'('+phe+')'))
        
    # make sense if we indeed find out that if cell has both cd66b+cd3+ you expect them to secrete both NE GzB right?
    per_necd66bcd31= []

    
    
    necd66b = f_data.loc[(f_data['NE_gate']==1)&(f_data['CD66b_gate']==1)].copy()
    ne_cd66bcd31 = necd66b.loc[(necd66b['CD31_gate']==1)].copy()
    ne_cd66bcd31neg = necd66b.loc[(necd66b['CD31_gate']==0)].copy()
    for feat, df in zip(['(CD31+)NE+CD66b+','(CD31-)NE+CD66b+'],[ne_cd66bcd31,ne_cd66bcd31neg]):
        per_necd66bcd31.append(fill_data(pd.DataFrame(df.groupby(by = ['Condition','Sample_NO']).apply(lambda df: df['CD3_gate'].sum()/len(df))),'per',feat+'CD3+'))
        per_necd66bcd31.append(fill_data(pd.DataFrame(df.groupby(by = ['Condition','Sample_NO']).apply(lambda df: 1-(df['CD3_gate'].sum()/len(df)))),'per',feat+'CD3-'))

        
    for d in [per_cd3_in_cd66b,per_negzb,per_necd66bcd31]:
        cd66b_cd3_data.append(pd.concat(d, ignore_index = True))
    cd66b_cd3_data = pd.concat(cd66b_cd3_data,ignore_index = True)
    cd66b_cd3_data = wide_to_long(long_to_wide(cd66b_cd3_data))
        
    return type1_data, type2_data,type3_data,cd31_data,cd66b_cd3_data

def compute_total(cnt_df,val_df):
    tot_df = val_df.copy().reset_index(drop = True)
    tot_df['v'] *= cnt_df.copy().reset_index(drop = True)['v']
    tot_df['Feature_name'] = tot_df.apply(lambda row: '(total)'+row['Feature_name'],axis = 1)
    return tot_df 

def compute_percentage(norm_df,denorm_df):
    perc_df = norm_df.copy().reset_index(drop = True)
    denorm_df = denorm_df.copy().reset_index(drop = True)
    perc_df['v'] /= denorm_df['v'] # rare boundary cases, debug later
    tag = denorm_df['Feature_name'].unique()[0].rpartition(']')[-1]
    perc_df['Feature_name'] = perc_df.apply(lambda row: '(perc in '+ tag+')'+row['Feature_name'],axis = 1)
    return perc_df
    

def get_perc():
    perc_data = []
    enzs = ['NE+','GzB+']
    phes = ['CD66b+','CD3+','CD66b-CD3-']
    
    for ch in enzs :
        for phe in phes:
            n = '[cnt]'+ch+phe
            d = '[cnt]'+ch
            # print(n,d)
            perc_data+=[compute_percentage(type2_data.loc[type2_data['Feature_name']==n].copy(),
                                          type2_data.loc[type2_data['Feature_name']==d].copy())]
            
    for phe in phes:
        n = '[CD31+cnt]'+phe
        d = '[cnt]CD31+'
        perc_data+=[compute_percentage(type3_data.loc[type3_data['Feature_name']==n].copy(),
                                      cd31_data.loc[cd31_data['Feature_name']==d].copy())]
    
    for phe in phes:
        for ch in enzs+['CD31+']:
            n = '['+ch+'cnt]'+phe
            d = '[cnt]'+phe
            perc_data+=[compute_percentage(type3_data.loc[type3_data['Feature_name']==n].copy(),
                                          type3_data.loc[type3_data['Feature_name']==d].copy())]
 
    perc_data = pd.concat(perc_data,
                                     join = 'inner',ignore_index = True).fillna(0).reset_index(drop = True)
    return perc_data


def long_to_wide(s):
    wide = s.pivot(index = ['Condition','Sample_NO'],columns = 'Feature_name', values = 'v').fillna(0)
    return wide.loc[:, ~(wide == 0).all()].copy() # clean up columns with all values == 0


def get_sum(gating_,chan,conds= ['HC','AMI T1','AMI T2','AMI T3']):
    data = f_data.loc[f_data['Condition'].isin(conds)].copy()
    for key in gating_:
        data = data.loc[data[key]==gating_[key]].copy()
    summary = data.groupby(by = ['Condition','Sample_NO'])[[chan]].sum()
    summary.columns = ['v']
        
    return summary.reset_index()


def filter_rare_populations(counts):
    ct = counts[counts.index.str.contains('[',regex = False)].copy()#.reset_index()
    cols = ct.columns.tolist()
    ct['Condition'] = ct.index
    #ct['Sample_NO'] = ct['Condition'].apply(lambda x: x[(x.find('[')+1):(x.find(']'))])
    ct['Condition'] = ct['Condition'].apply(lambda x: x[:x.find('[')])
    
    new_cols = []
    for i, col in enumerate(cols):
        colvar = ct.groupby(by = ['Condition'])[col].quantile(0.5)
        if colvar.max()>=10:
            new_cols.append([col,i])
    return new_cols

