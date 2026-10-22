from pathlib import Path
import numpy as np,pandas as pd,json,itertools
R=Path(__file__).resolve().parent;O=R/'results';out={};rng=np.random.default_rng(20261013)
for name,K in [('all',4),('angust',2)]:
 m=pd.read_csv(O/f'{name}_PCA.csv');G=np.load(O/f'{name}_kernels.npz')['G'];q=pd.read_csv(O/f'{name}_Q_final.csv').values
 stab=[]
 for run in range(1,6):
  z=np.loadtxt(R/'data'/f'{name}.snmf'/f'K{K}'/f'run{run}'/f'{name}_r{run}.{K}.Q');best=min(itertools.permutations(range(K)),key=lambda p:np.mean((q-z[:,p])**2));stab.append(float(np.sqrt(np.mean((q-z[:,best])**2))))
 out[name]={'Q_missingness_r':[float(np.corrcoef(q[:,k],m.original_missing)[0,1]) for k in range(K)],'Q_aligned_RMSE_vs_five_masked_runs':stab}
 # Assess sensitivity to sample missingness using residual ordination.
 S=np.column_stack([np.ones(len(m)),m.original_missing]);Q,_=np.linalg.qr(S);E=np.eye(len(m))-Q@Q.T;Gc=E@G@E;v,U=np.linalg.eigh(Gc);ix=np.argsort(v)[::-1];v=np.maximum(v[ix],0);C=U[:,ix]*np.sqrt(v);mm=m.copy();mm['PC1']=C[:,0];mm['PC2']=C[:,1];mm.to_csv(O/f'{name}_PCA_missing_adjusted.csv',index=False);out[name]['adjusted_PC_percent']=(v[:2]/v.sum()*100).tolist()
 # Freedman-Lane residual permutation; site blocks maintained for management.
 def covtest(labels,site=None):
  n=len(m);cats=list(dict.fromkeys(labels));Z=np.column_stack([(np.array(labels)==x).astype(float) for x in cats[1:]])
  S=np.column_stack([np.ones(n),m.original_missing])
  if site is None:blocks=[np.arange(n)]
  else:
   cats2=list(dict.fromkeys(site));S=np.column_stack([S]+[(np.array(site)==x).astype(float) for x in cats2[1:]]);blocks=[np.where(np.array(site)==x)[0] for x in cats2]
  Q0,_=np.linalg.qr(S);E=np.eye(n)-Q0@Q0.T;W=E@Z;u,s,_=np.linalg.svd(W,full_matrices=False);rank=int(sum(s>1e-8));U=u[:,:rank];Gr=E@G@E;tr=np.trace(Gr);ss=np.trace(U.T@Gr@U);df=n-S.shape[1]-rank;F=(ss/rank)/((tr-ss)/df);hits=0
  for b in range(9999):
   perm=np.arange(n)
   for ind in blocks:perm[ind]=rng.permutation(ind)
   # P^T times design bases, equivalent to permuting residual rows in Y.
   up=U[perm];qp=Q0[perm];ssb=np.trace(up.T@Gr@up);trb=np.trace(Gr)-np.trace(qp.T@Gr@qp);fb=(ssb/rank)/((trb-ssb)/df);hits+=fb>=F
  return {'F':float(F),'R2_total':float(ss/np.trace(G)),'partial_R2':float(ss/tr),'P':float((hits+1)/10000),'df1':rank,'df2':df}
 if name=='all':out[name]['species_adjusted_missing']=covtest(m.species)
 else:out[name]['site_adjusted_missing']=covtest(m.site);out[name]['management_adjusted_site_missing']=covtest(m.management,m.site)
# Holm correction for the three exploratory within-site tests.
a=json.loads((O/'analysis_results.json').read_text());tests=a['angust']['within_site'];ordered=sorted(tests,key=lambda k:tests[k]['P']);running=0
for i,k in enumerate(ordered):running=max(running,min(1,(3-i)*tests[k]['P']));tests[k]['P_Holm']=running
(O/'analysis_results.json').write_text(json.dumps(a,indent=2));(O/'diagnostics.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
