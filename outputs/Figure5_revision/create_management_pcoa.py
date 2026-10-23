from pathlib import Path
import pandas as pd, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[2]; O=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
t=pd.read_csv(R/'outputs/management_Hope34_excluded/corrected_pcoa_coordinates.csv').rename(columns={'Coordinate1':'PCo1','Coordinate2':'PCo2'})
t['site']=t.group.str.rsplit('_',n=1).str[0];t['management']=t.group.str.rsplit('_',n=1).str[1]
summary=json.loads((R/'outputs/management_Hope34_excluded/recalculation_summary.json').read_text());p=summary['pcoa_axis_percentages']
assert len(t)==136 and 'Hope34' not in t['sample'].tolist()
colors={'Deblois':'#6A51A3','Jonesboro':'#4D4D4D','Hope':'#3A923A'}
f,ax=plt.subplots(figsize=(7.5,5.4));f.subplots_adjust(left=.12,right=.97,top=.97,bottom=.27)
for (site,m),v in t.groupby(['site','management'],sort=False):
 ax.scatter(v.PCo1,v.PCo2,s=40,c=colors[site],marker='o' if m=='Managed' else '^',edgecolors='white',linewidths=.45)
ax.axhline(0,color='#dddddd',lw=.6,zorder=0);ax.axvline(0,color='#dddddd',lw=.6,zorder=0)
ax.set_xlabel(f'Coordinate 1 ({p[0]:.2f}%)');ax.set_ylabel(f'Coordinate 2 ({p[1]:.2f}%)');ax.set_aspect('equal',adjustable='datalim');ax.margins(.12)
for sp in ax.spines.values():sp.set_visible(True);sp.set_linewidth(.8);sp.set_color('#333333')
h=[Line2D([],[],marker='o',color='none',markerfacecolor=c,markeredgecolor='none',label=site) for site,c in colors.items()]
l=ax.legend(handles=h,loc='upper center',bbox_to_anchor=(.5,-.16),ncol=3,frameon=False,fontsize=10);ax.add_artist(l)
h=[Line2D([],[],marker=mk,color='none',markerfacecolor='#666666',markeredgecolor='none',label=m) for m,mk in [('Managed','o'),('Unmanaged','^')]]
ax.legend(handles=h,loc='upper center',bbox_to_anchor=(.5,-.25),ncol=2,frameon=False,fontsize=10)
for ext in ['png','pdf','svg']:f.savefig(O/f'Figure5_management_PCoA.{ext}',dpi=400,facecolor='white')
f.savefig(O/'Figure5_preview.png',dpi=150,facecolor='white')
print(len(t),p)
