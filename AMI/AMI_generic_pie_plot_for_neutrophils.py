# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 21:23:13 2022

@author: Lu Ri
"""
def do_pieplot(piedata):
    fig, axes = plt.subplots(3,4, figsize = (12,18),dpi = 200)
    plt.subplots_adjust(hspace = 2)
    rows = piedata.index.tolist()
    cols = conditions

    for i in range(3):
        cat = rows[i].partition('[')[-1]
        cat = cat[:cat.find('+')]#[(rows[i].find('[')+1)]
        pal = mypalette[cat]
        for j in range(4):
            posvar = piedata[cols[j]][rows[i]]
            axes[i,j].pie([posvar,1-posvar],radius = 1.2,colors = [pal[cols[j]], 'ghostwhite'],startangle = 90)
    return fig

fig = do_pieplot(piedata)
fig.savefig(export_dir+'neutrophil_pie.png',bbox_inches = 'tight')