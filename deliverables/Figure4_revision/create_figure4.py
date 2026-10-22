from pathlib import Path
import json
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
R=Path(__file__).resolve().parents[2];O=R/'analysis_v2/results';F=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
spc={'angustifolium':'#0072B2','myrtilloides':'#D55E00','highbush reference':'#009E73'}
spl={'angustifolium':'V. angustifolium','myrtilloides':'V. myrtilloides','highbush reference':'V. corymbosum'}
sitec={'Deblois':'#6A51A3','Jonesboro':'#4D4D4D','Hope':'#3A923A','Fort Kent':'#E69F00','Reference':'#D9D9D9'}
mgc={'Managed':'#36648B','Unmanaged':'#E5BA78','Reference':'#D9D9D9'}
qc=['#4477AA','#EE6677','#228833','#CCBB44']
m=pd.read_csv(O/'all_PCA.csv');q=pd.read_csv(O/'all_Q_final.csv').to_numpy();h=pd.read_csv(O/'all_ward_merge.csv');order=pd.read_csv(O/'all_ward_order.csv')['index'].to_numpy()-1
assert len(m)==129 and q.shape==(129,4) and np.allclose(q.sum(1),1,atol=1e-5)
fig=plt.figure(figsize=(13,11.5))
def title(x,y,letter,text):
 fig.text(x,y,letter,fontsize=17,weight='bold',va='top')
# A: cross-entropy
ax=fig.add_axes([.075,.68,.285,.25]);d=pd.read_csv(O/'all_cross_entropy.csv');s=d.groupby('K').cross_entropy.agg(['mean','std'])
for k,z in d.groupby('K'):
 ax.scatter(np.repeat(k,len(z))+np.linspace(-.07,.07,len(z)),z.cross_entropy,c='#B7BCC4',s=18,zorder=2)
ax.errorbar(s.index,s['mean'],s['std'],fmt='o-',capsize=4,color='#333333',lw=1.2,zorder=3)
ax.set_xticks(range(1,7));ax.set_xlabel('Number of ancestry components (K)');ax.set_ylabel('Cross-entropy');ax.set_ylim(.562,.587)
title(.035,.978,'A','Ancestry model selection')
# D: PCA
ax=fig.add_axes([.075,.345,.285,.24]);per=json.load(open(O/'analysis_results.json'))['all']['PC_percent']
for sp,col in spc.items():
 z=m[m.species==sp];ax.scatter(z.PC1,z.PC2,c=col,s=35,edgecolor='white',linewidth=.45,label=spl[sp],alpha=.9)
ax.axhline(0,color='#dddddd',lw=.7,zorder=0);ax.axvline(0,color='#dddddd',lw=.7,zorder=0)
ax.set_xlabel(f'PC1 ({per[0]:.2f}%)');ax.set_ylabel(f'PC2 ({per[1]:.2f}%)');title(.035,.632,'D','Principal component analysis')
ax.legend(frameon=False,prop={'size':7.5,'style':'italic'},loc='upper right',ncol=1)
# B: hierarchical clustering
ax=fig.add_axes([.365,.335,.63,.615],projection='polar');n=len(m);mx=h.height.max();ang=np.linspace(.058,2*np.pi-.058,n);ax.set_theta_offset(np.pi/2);ax.set_theta_direction(-1)
coords={-int(ind+1):(ang[pos],1.) for pos,ind in enumerate(order)}
for i,r in h.iterrows():
 t1,r1=coords[int(r.left)];t2,r2=coords[int(r.right)];rh=1-r.height/mx
 ax.plot([t1,t1],[r1,rh],c='#626972',lw=.55);ax.plot([t2,t2],[r2,rh],c='#626972',lw=.55)
 tt=np.linspace(t1,t2,max(8,int(abs(t2-t1)*70)));ax.plot(tt,np.full_like(tt,rh),c='#626972',lw=.55);coords[i+1]=((t1+t2)/2,rh)
step=(ang[-1]-ang[0])/(n-1)
for pos,ind in enumerate(order):
 r=m.iloc[ind];t=ang[pos];site='Reference' if r.species=='highbush reference' else r.site
 ax.scatter([t],[1],s=6,c=spc[r.species],zorder=4)
 ax.bar(t,.023,width=step*.9,bottom=1.018,color=sitec[site],edgecolor='none');ax.bar(t,.023,width=step*.9,bottom=1.05,color=mgc[r.management],edgecolor='none')
 rot=90-np.degrees(t);ha='left'
 if np.pi<t<2*np.pi:rot+=180;ha='right'
 ax.text(t,1.088,r['sample'],rotation=rot,rotation_mode='anchor',ha=ha,va='center',fontsize=7.0,color=spc[r.species])
ax.set_ylim(0,1.34);ax.set_axis_off();title(.415,.978,'B','Hierarchical clustering · 129 samples')
fig.legend(handles=[Patch(facecolor=c,label=(r'$\it{V.\ corymbosum}$' if k=='Reference' else k)) for k,c in sitec.items()],loc='upper left',bbox_to_anchor=(.435,.348),title='Inner ring: location',frameon=False,ncol=2,fontsize=8,title_fontsize=9)
fig.legend(handles=[Patch(facecolor=c,label=(r'$\it{V.\ corymbosum}$' if k=='Reference' else k)) for k,c in mgc.items()],loc='upper left',bbox_to_anchor=(.745,.348),title='Outer ring: management',frameon=False,ncol=1,fontsize=8,title_fontsize=9)
# C: ancestry proportions in clustering order
ax=fig.add_axes([.075,.080,.88,.135]);bottom=np.zeros(n)
for j in range(4):
 v=q[order,j];ax.bar(np.arange(n),v,bottom=bottom,width=1,color=qc[j],linewidth=0);bottom+=v
ax.set_xlim(-.5,n-.5);ax.set_ylim(0,1);ax.set_ylabel('Ancestry proportion');ax.set_yticks([0,.5,1]);ax.set_xticks(np.arange(n));ax.set_xticklabels(m.iloc[order]['sample'],rotation=90,fontsize=5.5);ax.tick_params(axis='x',length=0,pad=2)
for tick,sp in zip(ax.get_xticklabels(),m.iloc[order].species):tick.set_color(spc[sp])
title(.035,.258,'C','Ancestry proportions · K = 4')
for ext in ['png','pdf','svg']:fig.savefig(F/f'Figure4_ABCD.{ext}',dpi=450,facecolor='white')
fig.savefig(F/'Figure4_ABCD_Docs.png',dpi=220,facecolor='white')
plt.close(fig)
m.iloc[order].assign(**{f'component_{j+1}':q[order,j] for j in range(4)}).to_csv(F/'Figure4_sample_order_and_ancestry.csv',index=False)
print('Saved Figure4_ABCD PNG/PDF/SVG')
