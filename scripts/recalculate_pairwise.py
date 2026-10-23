from pathlib import Path
import json, hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'outputs/Management_differentiation_preview';O.mkdir(exist_ok=True)
f=R/'outputs/management_Hope34_excluded/corrected_management_input.csv';t=pd.read_csv(f);X=t.iloc[:,2:].to_numpy();g=t.group.to_numpy();ch=np.array(['_'.join(c.split('_')[:-2]) for c in t.columns[2:]])
assert len(t)==136 and 'Hope34' not in t['sample'].tolist() and len(np.unique(ch))==12
D=np.array([(X!=r).sum(axis=1) for r in X],float)
def components(d,labels):
 n=len(d);labs=np.unique(labels);k=len(labs);ns=np.array([(labels==a).sum() for a in labs]);total=d.sum()/(2*n);within=sum(d[np.ix_(labels==a,labels==a)].sum()/(2*ns[i]) for i,a in enumerate(labs));msw=within/(n-k);msa=(total-within)/(k-1);n0=(n-(ns**2).sum()/n)/(k-1);va=(msa-msw)/n0
 return va,msw
def phi(d,labels):
 va,vw=components(d,labels);return va/(va+vw)
labels=list(dict.fromkeys(g));matrix=pd.DataFrame(np.zeros((6,6)),index=labels,columns=labels)
for i,a in enumerate(labels):
 for j,b in enumerate(labels[:i]):
  keep=np.isin(g,[a,b]);matrix.iloc[i,j]=matrix.iloc[j,i]=phi(D[np.ix_(keep,keep)],g[keep])
matrix.to_csv(O/'pairwise_PhiPT_corrected.csv')
