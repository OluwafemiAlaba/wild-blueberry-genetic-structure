from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[1];O=R/'outputs/management_Hope34_excluded';O.mkdir(exist_ok=True)
f=R/'inputs/amova/managed_unmanaged.xlsx';w=openpyxl.load_workbook(f,read_only=True,data_only=True);rr=list(w['Sheet1'].values);ids=np.array([r[0] for r in rr[3:]]);groups=np.array([r[1] for r in rr[3:]]);X=np.array([r[2:] for r in rr[3:]],int);headers=list(rr[2][2:]);D=np.array([(X!=r).sum(axis=1) for r in X],float)
def ordination(d):
 n=len(d);J=np.eye(n)-np.ones((n,n))/n;G=-.5*J@d@J;v,U=np.linalg.eigh(G);ix=np.argsort(v)[::-1];v=v[ix];U=U[:,ix];return U[:,:2]*np.sqrt(v[:2]),v,100*v[:2]/v[v>1e-8].sum()
oldcoord,oldv,oldpct=ordination(D);saved=np.array(list(w['PCoA'].iter_rows(min_row=33,max_row=169,min_col=2,max_col=3,values_only=True)),float)
cor=[float(np.corrcoef(oldcoord[:,i],saved[:,i])[0,1]) for i in range(2)];assert all(abs(c)>0.999999999 for c in cor)
keep=ids!='Hope34';assert sum(~keep)==1
excluded={'sample':'Hope34','historical_group':str(groups[~keep][0]),'reason':'Author confirmed V. myrtilloides; excluded from V. angustifolium management analysis.'}
ids=ids[keep];groups=groups[keep];X=X[keep];D=D[np.ix_(keep,keep)];n=len(ids);labs=list(dict.fromkeys(groups));g=np.array([labs.index(v) for v in groups]);sizes=np.bincount(g);k=len(labs);sst=D.sum()/(2*n);n0=(n-(sizes*sizes).sum()/n)/(k-1)
def stat(gg):
 ssw=sum(D[np.ix_(gg==a,gg==a)].sum()/(2*sizes[a]) for a in range(k));ssa=sst-ssw;msw=ssw/(n-k);msa=ssa/(k-1);va=(msa-msw)/n0;return ssa,ssw,msa,msw,va,va/(va+msw)
a=stat(g);rng=np.random.default_rng(20261009);perms=9999;exceed=0
for i in range(perms):exceed+=stat(rng.permutation(g))[-1]>=a[-1]
p=(exceed+1)/(perms+1)
coords,eigs,pct=ordination(D)
for i,c in enumerate(cor):coords[:,i]*=np.sign(c)
df=pd.DataFrame({'sample':ids,'group':groups,'Coordinate1':coords[:,0],'Coordinate2':coords[:,1]});df.to_csv(O/'corrected_pcoa_coordinates.csv',index=False)
inputdf=pd.DataFrame(X,columns=headers);inputdf.insert(0,'group',groups);inputdf.insert(0,'sample',ids);inputdf.to_csv(O/'corrected_management_input.csv',index=False)
table=pd.DataFrame([['Among groups',k-1,a[0],a[2],a[4],100*a[5]],['Within groups',n-k,a[1],a[3],a[3],100*(1-a[5])],['Total',n-1,sst,None,a[4]+a[3],100]],columns=['Source','df','SS','MS','Variance_component','Percent']);table.to_csv(O/'corrected_amova.csv',index=False)
out={'source_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'excluded':excluded,'n':n,'markers':X.shape[1],'groups':dict(zip(labs,map(int,sizes))),'PhiPT':a[-1],'P':p,'permutations':perms,'exceedances':int(exceed),'seed':20261009,'permutation_scheme':'Unrestricted permutation of sample group labels, preserving group sizes; reproduces historical six-group design, not a site-controlled management test.','pcoa_axis_percentages':pct.tolist(),'historical_pcoa_axis_correlations':cor,'historical_pcoa_axis_percentages_reproduced':oldpct.tolist(),'distance':'Count of unequal genotype codes; PCoA G=-0.5 J D J (equivalent to classical PCoA of square-root mismatch distance).','limits':['Unknown missing-call filling method inherited from source matrix.','TP canopy representative and independence questions remain.','No other samples or genotype values changed.','PCoA axis units use unscaled principal coordinates; historical workbook used a common rescaling.']};(O/'recalculation_summary.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
