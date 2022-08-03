# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 17:06:30 2022

@author: Lu Ri
"""
def order_columns(counts):
    cols = counts.columns.tolist()
    def indices(l, filtr=lambda x: bool(x)):
        return [i for i,x in enumerate(l) if filtr(x)]
    
    ### main phenotype clusters
    cd66bcd3 =  [x for i, x in enumerate(cols) if i in set(indices(cols,filtr=lambda x: '|CD66b|CD3|' in x))]    
    rem = [x for i, x in enumerate(cols) if i not in set(indices(cols,filtr=lambda x: '|CD66b|CD3|' in x))]    
    cd66b = [x for i, x in enumerate(rem) if i in set(indices(rem,filtr=lambda x: '|CD66b|' in x))]
    cd3 = [x for i, x in enumerate(rem) if i in set(indices(rem,filtr=lambda x: '|CD3|' in x))]
    double_neg = [x for i, x in enumerate(rem) if i in set(indices(rem,filtr=lambda x: ('|CD3|' not in x) and ('|CD66b|' not in x)))]
    
    prim_clusters = [cd66b,cd66bcd3,cd3,double_neg] ### in order
    ### secondary phenotype clusters cd31+ and cd31- (much simpler than ne/gzb)
    ordered_cols = []
    
    for clus in prim_clusters:
        cd31_pos = [x for i, x in enumerate(clus) if i in set(indices(clus,filtr=lambda x: '|CD31|' in x))]
        cd31_neg = [x for i, x in enumerate(clus) if i in set(indices(clus,filtr=lambda x: '|CD31|' not in x))]
        sec_clusters = [cd31_pos,cd31_neg]
        for sc in sec_clusters:
            ordered_cols += counts[sc].sum().sort_values(ascending = False).index.tolist()
        
    ## then order under secondary clusters, by the counts 
    return ordered_cols

def normalize(df_orig):
    df = df_orig.copy()
    cols = df.columns.tolist()
    for col in cols:
        if df[col].max()!=0:
            df[col] = 50+1000*(df[col]-df[col].min())/(df[col].max()-df[col].min()) ## to get some difference from 0
    return df
    
def prepare_dot_data(f_data, conds = ['HC','AMI T1','AMI T2','AMI T3'],channels = ['NE','GzB','CD31']):
    
    # you can set a column named 'order' to help sort decending at last
    dummyconds = [(i+1)*' ' for i in range(len(conds)-1)]
    counts = []
    features = []
    for num, cond in enumerate(conds):
        s_data = f_data.loc[(f_data['Condition']==cond)].copy()
        for spl_no in s_data['Sample_NO'].unique().tolist():

            data = s_data.loc[(s_data['Sample_NO']==spl_no)].copy()
            cnt = data.groupby(by = 'Phenotype')[['Phenotype']].count()
            chan_sums = cnt.copy()
            chan_sums['Phenotype'] = 0            
            if num <len(conds)-1:
                dummycond = chan_sums.copy().T.rename(index = {'Phenotype':dummyconds[num]})
                
            for ch in channels:
                chvar = data.groupby(by = 'Phenotype')[[ch]].median()
                chvar.columns = ['Phenotype']
                chan_sums += chvar
            cnt = cnt.T#.reset_index(drop = True)
            chan_sums = chan_sums.T#.reset_index(drop = True)
            new_idx = cond + '[' + str(spl_no) + ']'
            cnt = cnt.rename(index = {'Phenotype':new_idx})
            chan_sums = chan_sums.rename(index = {'Phenotype':new_idx})
            counts.append(cnt)
            features.append(chan_sums)
        if num <len(conds)-1:
            counts.append(dummycond)
            #features.append(dummycond)
    counts = pd.concat(counts).fillna(0)
    features = pd.concat(features).fillna(0)
    features = normalize(features)  # automatically column_wise normalization
    features.reindex(counts.index).fillna(0)
    new_col_list = order_columns(counts.copy())
    counts = counts[new_col_list]
    features = features[new_col_list]
    
    return counts, features, new_col_list


def dot_plot(counts,features):

    sizes = [1,5,13]
    size_gating = [10,100,1000]
    def count_mapping(x):
        if x <size_gating[0]:
            return sizes[0] * x/size_gating[0]
        
        elif x <size_gating[1]:
            return sizes[0]+ ((x-size_gating[0])/(size_gating[1]-size_gating[0]))*(sizes[1]-sizes[0])
        
        elif x < size_gating[2]:
            return sizes[1]+ ((x-size_gating[1])/(size_gating[2]-size_gating[1]))*(sizes[2]-sizes[1])
        
        else:
            return sizes[2]+2*x/size_gating[2]
    
    def feature_mapping(x):
        return log10(x)
    
    counts = counts.applymap(count_mapping)
    # features = features.applymap(feature_mapping)
    
    num_col = counts.shape[1]
    num_row = counts.shape[0]
    x = np.arange(num_col)
    y = np.arange(num_row)
    X, Y = np.meshgrid(x, y)
    XY = np.column_stack((X.ravel(), Y.ravel()))
    WW = counts.to_numpy()/20

    ############################################################
    #mypalette = [cmc.cork_r(0),cmc.nuuk(200),cmc.roma(30),cmc.vik_r(20),cmc.bam(0)]
    mypalette = ['turquoise','green','gold','crimson','orchid']
    # master_pal = {'green':[cmc.cork_r,(0,115)],
    #               'off_green':[cmc.bam_r,(0,115)],
    #               'blue':[cmc.cork,(0,115)],
    #               'light_blue':[cmc.roma_r,(0,115)],
    #               'off_blue':[cmc.oleron,(0,115)], # recommended for 2 cond sequential
    #               'brown':[cmc.oleron,(128,256)],              
    #               'orange':[cmc.roma,(0,115)],
    #               'red':[cmc.vik_r,(0,115)],
    #               'off_red':[cmc.bilbao_r,(0,200)],
    #               'purple':[cmc.bam,(0,115)],
    #               'karki':[cmc.broc_r,(0,115)]}

    fig = plt.figure(figsize = (18,20),constrained_layout=True)
    gs = fig.add_gridspec(22,20,left=0.05, right=0.95,
                            hspace = 0.01,wspace=0.5)
    ax = fig.add_subplot(gs[4:, :18])
    ax_cbar = fig.add_subplot(gs[4:16, 18])
    ax_size = fig.add_subplot(gs[16:, 18:])
    ax_phe_leg= fig.add_subplot(gs[2:4, :18])
    ax_tree = fig.add_subplot(gs[0:2, :18])
    
    ################ prepare the phenotype label ##############
    new_col_list = counts.columns.tolist()
    gate_list = [c for c in f_data.columns.tolist()if '_gate' in c]+['Phenotype']
    channels = f_data.columns.tolist()[:5]
    df = f_data[gate_list]
    size = 100
    marker = 's'
    for i, phe in enumerate(new_col_list):
        #ax_phe_leg.text(i,4.5,str(i), ha = 'center', va = 'center',color = 'k')
        chan_code = df.loc[df['Phenotype']== phe].iloc[-1,:-1].tolist()
        ys = [i for i in range(len(chan_code))][::-1]
        for j, code in enumerate(chan_code):
            if code == 1:
                ax_phe_leg.scatter(i,ys[j],s = size, marker = marker, c = mypalette[j])
            else:
                ax_phe_leg.scatter(i,ys[j],s = size, marker = marker, facecolors = 'none', edgecolors = mypalette[j])  
                
    ######################## tree ##############################
    color = 'k'
    for x in np.arange(len(new_col_list)):
        ax_tree.plot([x,x], [0,1],color = color)  
        
    for xlims in [(np.array([0,3])+i*4)for i in range(8)]:
        ax_tree.plot(xlims, [1,1], color = color)
        
    for x in [1.5 + i *4 for i in range(8)]:
        ax_tree.plot([x,x], [1,2],color = color)
    
    for xlims in [(np.array([1.5,5.5])+i*8)for i in range(4)]:
        ax_tree.plot(xlims,[2,2],color = color)
    
    l = [3.5,11.5]
    for x in l:
        ax_tree.plot([x,x],[2,3],color = color)    
    ax_tree.plot(l,[3,3],color = color)
    cent1 = l[0]+0.5*(l[1]-l[0])
    ax_tree.plot([cent1,cent1], [3,4],color = color) 
    l = [19.5,27.5]
    for x in l:
        ax_tree.plot([x,x], [2,4], color = color)  
    ax_tree.plot([cent1,27.5], [4,4],color = color)
    
    
    ######################### plot main ##############################
    ec = EllipseCollection(WW, WW, 0, units='x', offsets=XY,cmap = 'magma', #cmc.batlowK,
                           transOffset=ax.transData)
    ec.set_array(features.to_numpy().ravel())
    ax.add_collection(ec)
    ax.set_xlim(-1,num_col)
    ax.set_ylim(-1,num_row)

    x_ticks = ['']+counts.columns.tolist()
    y_ticks = ['']+counts.index.tolist()
    ax.xaxis.set_ticks(np.arange(num_col+1)-1)
    ax.yaxis.set_ticks(np.arange(num_row+1)-1)
    ax.set_xticklabels([])
    ax.set_yticklabels(y_ticks)
    ax.xaxis.set_ticks_position('top') 
    cbar = fig.colorbar(ec, cax = ax_cbar)
    cbar.set_label('Normalized median expression [a.u.]')
    
    for i,s in enumerate(sizes):
        ax_size.scatter(0,i, s = 4*s**2, color = 'k')
        ax_size.text(0.3,i-0.05, str(size_gating[i]))
    
    ax_size.set_ylim(-0.25,len(sizes)+0.5)
    ax_size.set_xlim(-0.25,1.4)
        
    
    hlines = [9,19,29]
    for hl in hlines:
        ax.axhline(y = hl, lw = 0.5,color = 'silver')

    vlines = [7.5,15.5,23.5]
    for vl in vlines:
        ax.axvline(x = vl,lw = 0.5,color = 'silver')
        
    ax_phe_leg.set_xlim(-1,num_col)
    ax_phe_leg.set_ylim(-0.3,len(channels)-0.7)
    ax_phe_leg.set_yticklabels(channels[::-1])
    ax_phe_leg.yaxis.set_ticks(np.arange(len(channels)))
    ax_tree.set_xlim(-1,num_col)
    ax_tree.set_ylim(0,4.5)
    
    
    #ax_phe_leg.get_xaxis().set_visible(False)
    ax_phe_leg.xaxis.set_ticks_position('top') 
    ax_phe_leg.xaxis.set_ticks(np.arange(num_col+1)-1)
    ax_phe_leg.set_xticklabels(['']+[i for i in range(len(new_col_list))])
    ax_phe_leg.tick_params(axis=u'both', which=u'both',length=0) # to hide the ticks
    for key, spine in ax_phe_leg.spines.items():
        spine.set_visible(False)
    for a in [ax_cbar,ax_size,ax_tree]:
        a.axis('off')
    return fig, ax


export_dir = 'D:/CAMP Biochemical/20220523 AMI_TOTAL_ANALYSIS/All_automated_figures/'

counts,features,new_col_list = prepare_dot_data(f_data)

# fig, ax = dot_plot(counts.iloc[::-1], features.iloc[::-1])
# fig.savefig(export_dir+'dotplot.png',bbox_inches = 'tight')


    


