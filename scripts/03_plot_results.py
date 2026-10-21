#!/usr/bin/env python3
"""Reconstructed plotting code with portable paths.
Replots tables using the documented analysis choices; visual styling may differ
from historical figures. ENC values and omega are imported archived estimates.
"""
from pathlib import Path
from io import StringIO
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
from scipy.stats import linregress
from Bio import Phylo
root=Path(__file__).resolve().parents[1];res=root/'results';out=res/'figures';out.mkdir(exist_ok=True)
def save(fig,name):
    fig.tight_layout()
    for ext in ['pdf','png']:fig.savefig(out/f'{name}.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)
# Euclidean distance and Ward linkage; unscaled species mean RSCU.
d=pd.read_csv(res/'rscu_medias_14_especies.tsv',sep='\t')
m=d.pivot(index='species',columns='codon',values='mean_RSCU')
ri=leaves_list(linkage(pdist(m.to_numpy()),method='ward'))
ci=leaves_list(linkage(pdist(m.to_numpy().T),method='ward'))
m=m.iloc[ri,ci];fig,ax=plt.subplots(figsize=(18,7))
im=ax.imshow(m.to_numpy(),aspect='auto',cmap='YlOrRd',vmin=0,vmax=6)
ax.set_yticks(range(len(m)),m.index,fontstyle='italic',fontsize=9)
ax.set_xticks(range(len(m.columns)),m.columns,rotation=90,fontsize=7)
ax.set_xlabel('Codon');fig.colorbar(im,ax=ax,label='Mean RSCU');save(fig,'Figure2_RSCU')
# CodonW outputs remain archived inputs until an external rerun is completed.
e=pd.read_csv(root/'data/reference/Supplementary_ENC_GC3s_30.tsv',sep='\t')
e['species']=e.species.str.replace(r'^Helicoverpa armigera.*','Helicoverpa armigera',regex=True)
e=e.groupby('species',as_index=False)[['ENC','GC3s']].mean()
fig,ax=plt.subplots(figsize=(8,6));x=np.linspace(.001,.999,500)
ax.plot(x,2+x+29/(x*x+(1-x)**2),color='grey',label='Expected ENC')
for i,r in e.iterrows():ax.scatter(r.GC3s,r.ENC,color=plt.cm.tab20(i),label=r.species)
ax.set(xlabel='GC3s (fraction)',ylabel='ENC');ax.legend(bbox_to_anchor=(1.02,1),fontsize=7);save(fig,'Figure3_ENC')
s=pd.read_csv(res/'species_composition.tsv',sep='\t')
fig,ax=plt.subplots(figsize=(8,6))
for i,r in s.iterrows():ax.scatter(r.PR2_G3_G3C3,r.PR2_A3_A3T3,color=plt.cm.tab20(i),label=r.species)
ax.axhline(.5,color='grey',ls='--');ax.axvline(.5,color='grey',ls='--')
ax.set(xlim=(0,1),ylim=(0,1),xlabel='G3 / (G3 + C3)',ylabel='A3 / (A3 + T3)')
ax.legend(bbox_to_anchor=(1.02,1),fontsize=7);save(fig,'Figure4_PR2')
f=linregress(s.GC3,s.GC12);fig,ax=plt.subplots(figsize=(8,6))
for i,r in s.iterrows():ax.scatter(r.GC3,r.GC12,color=plt.cm.tab20(i),label=r.species)
x=np.linspace(s.GC3.min(),s.GC3.max(),100);ax.plot(x,f.intercept+f.slope*x)
ax.set(xlabel='GC3 (%)',ylabel='GC12 (%)');ax.legend(bbox_to_anchor=(1.02,1),fontsize=7)
ax.text(.03,.97,f'R squared = {f.rvalue**2:.4f}; p = {f.pvalue:.5f}',transform=ax.transAxes,va='top');save(fig,'Figure5_Neutrality')
g=pd.read_csv(res/'Supplementary_S1.tsv',sep='\t');fig,ax=plt.subplots(figsize=(9,5))
ax.bar(g.gene,g.omega,color='#365C70');ax.set(xlabel='Gene',ylabel='dN/dS (archived CODEML M0)');save(fig,'Figure6_M0')
# Trees are re-rendered as cladograms with original support labels retained.
for label,path,fmt in [('ML','PCG13_31_ML.treefile','newick'),('BI','PCG13_31_BI.con.tre','nexus')]:
    text=(root/f'data/reference/phylogeny_31/{path}').read_text()
    if fmt=='nexus':
        text=re.sub(r'(\))\[&prob=([^,\]]+)[^\]]*\]', lambda m:m[1]+m[2], text)
        text=re.sub(r'\[[^\]]*\]','',text)
    tree=Phylo.read(StringIO(text),fmt)
    outgroup=next(t for t in tree.get_terminals() if 'MN793322' in t.name)
    tree.root_with_outgroup(outgroup)
    for clade in tree.find_clades():clade.branch_length=1
    fig,ax=plt.subplots(figsize=(11,12));Phylo.draw(tree,axes=ax,do_show=False)
    ax.set_xlabel('Cladogram');save(fig,f'Figure1_{label}_cladogram')
print('Figures written to',out)
