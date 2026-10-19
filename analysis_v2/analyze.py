from pathlib import Path
import pandas as pd,numpy as np,json
R=Path(__file__).resolve().parent;O=R/'results';M=pd.read_csv(R/'data/analysis_samples.csv');z=np.load(R/'data/genotypes.npz');X=z['X'][:,z['ld']];V=pd.read_csv(R/'data/variant_manifest.csv').query('LD_retained').reset_index(drop=True);res={}
def pca(A):
 p=np.mean(A,axis=0)/2;ok=(p>0)&(p<1)&(np.std(A,axis=0)>0);Y=(A[:,ok]-2*p[ok])/np.sqrt(2*p[ok]*(1-p[ok]));Y-=Y.mean(axis=0);G=Y@Y.T;v,U=np.linalg.eigh(G);ix=np.argsort(v)[::-1];v=np.maximum(v[ix],0);U=U[:,ix];C=U*np.sqrt(v);return C,v/v.sum()*100,G
for name,sel in [('all',np.ones(len(M),bool)),('angust',M.species.eq('angustifolium').values),('myrt',M.species.eq('myrtilloides').values)]:
 meta=M[sel].reset_index(drop=True);raw=X[sel];mean=np.nanmean(raw,axis=0);filled=np.where(np.isnan(raw),mean,raw)
 if name!='myrt':A=np.loadtxt(R/'data'/f'{name}_final.lfmm_imputed.lfmm')
 else:A=filled
 assert A.shape==raw.shape;assert np.all(A[np.isfinite(raw)]==raw[np.isfinite(raw)])
 C,pct,G=pca(A);Cm,pm,Gm=pca(filled)
 for j in range(C.shape[1]):
  if C[np.argmax(abs(C[:,j])),j]<0:C[:,j]*=-1
 out=meta.copy()
 for j in range(min(10,len(meta))):out[f'PC{j+1}']=C[:,j]
 out.to_csv(O/f'{name}_PCA.csv',index=False);np.savez(O/f'{name}_kernels.npz',G=G,Gmean=Gm,C=C)
 # RV coefficient compares centered genotype kernels without arbitrary eigenvector signs.
 rv=float(np.sum(G*Gm)/np.sqrt(np.sum(G*G)*np.sum(Gm*Gm)))
 missingcorr=[float(np.corrcoef(C[:,i],meta.original_missing)[0,1]) for i in range(2)]
 res[name]={'n':len(meta),'PC_percent':pct[:10].tolist(),'imputed_vs_mean_kernel_RV':rv,'PC_missingness_correlations':missingcorr}
 np.save(O/f'{name}_PC_percent.npy',pct)
# Kernel-based multivariate tests; genotype standardization per species.
rng=np.random.default_rng(20261012)
def modeltest(G,labels,strata=None,perms=9999):
 n=len(G);group=np.asarray(labels);cats=list(dict.fromkeys(group));Z=np.column_stack([np.ones(n)]+[(group==v).astype(float) for v in cats[1:]])
 if strata is None:S=np.ones((n,1));blocks=[np.arange(n)]
 else:
  st=np.asarray(strata);sc=list(dict.fromkeys(st));S=np.column_stack([(st==v).astype(float) for v in sc]);blocks=[np.where(st==v)[0] for v in sc]
 Q0=np.linalg.qr(S)[0];rank0=np.linalg.matrix_rank(S);H0=Q0[:,:rank0]@Q0[:,:rank0].T;E=np.eye(n)-H0;W=E@Z[:,1:];u,s,_=np.linalg.svd(W,full_matrices=False);rank=int(sum(s>1e-8));U=u[:,:rank];Gr=E@G@E;total=float(np.trace(G));tr=float(np.trace(Gr));effect=float(np.trace(U.T@Gr@U));df2=n-rank0-rank;F=(effect/rank)/((tr-effect)/df2);gt=0
 for b in range(perms):
  perm=np.arange(n)
  for ix in blocks:perm[ix]=rng.permutation(ix)
  up=U[perm];ss=float(np.trace(up.T@Gr@up));fp=(ss/rank)/((tr-ss)/df2);gt+=fp>=F
 return {'F':F,'df1':rank,'df2':df2,'R2_total':effect/total,'partial_R2':effect/tr,'P':(gt+1)/(perms+1),'permutations':perms,'exceedances':int(gt)}
for name in ['all','angust','myrt']:
 meta=pd.read_csv(O/f'{name}_PCA.csv');D=np.load(O/f'{name}_kernels.npz');G=D['G']
 if name=='all':res[name]['species_test']=modeltest(G,meta.species)
 elif name=='angust':
  res[name]['site_test']=modeltest(G,meta.site);res[name]['management_adjusted_site']=modeltest(G,meta.management,meta.site);res[name]['management_meanfill_sensitivity']=modeltest(D['Gmean'],meta.management,meta.site)
  # Pairwise observed-call Euclidean distances, centered; sensitivity without imputed genotypes.
  raw=X[M.species.eq('angustifolium').values];dist=np.empty((len(raw),len(raw)))
  for i in range(len(raw)):
   diff=raw-raw[i];dist[i]=np.nanmean(diff*diff,axis=1)
  J=np.eye(len(raw))-1/len(raw);Gi=-.5*J@dist@J;res[name]['management_pairwise_observed_sensitivity']=modeltest(Gi,meta.management,meta.site)
  # Dispersion diagnostic: distances to each site-by-management centroid, scalar site-adjusted test.
  C=D['C'];dists=np.empty(len(meta));labs=meta.site+'_'+meta.management
  for label in labs.unique():
   ix=np.where(labs==label)[0];dists[ix]=np.linalg.norm(C[ix]-C[ix].mean(axis=0),axis=1)
  y=dists-dists.mean();res[name]['dispersion_site_management']=modeltest(np.outer(y,y),meta.management,meta.site)
  res[name]['within_site']={}
  for site in meta.site.unique():
   ix=np.where(meta.site==site)[0];g=G[np.ix_(ix,ix)];j=np.eye(len(ix))-1/len(ix);res[name]['within_site'][site]=modeltest(j@g@j,meta.management.iloc[ix],perms=4999)
 elif name=='myrt':
  keep=meta.site.ne('Hope').values;g=G[np.ix_(keep,keep)];n=sum(keep);j=np.eye(n)-1/n;res[name]['two_site_test']=modeltest(j@g@j,meta.site[keep])
# Within-species descriptive diversity from observed, not imputed, calls.
div=[]
for (sp,site,mg),sub in M.groupby(['species','site','management']):
 ix=sub.index.to_numpy();a=X[ix];ho=np.sum(a==1,axis=1)/np.sum(np.isfinite(a),axis=1);p=np.nanmean(a,axis=0)/2;he=2*p*(1-p);boots=np.mean(rng.choice(ho,(2000,len(ho)),replace=True),axis=1)
 div.append({'species':sp,'site':site,'management':mg,'n':len(ix),'Ho':float(ho.mean()),'Ho_low':float(np.quantile(boots,.025)),'Ho_high':float(np.quantile(boots,.975)),'He_recoded':float(np.nanmean(he))})
pd.DataFrame(div).to_csv(O/'diversity_observed_calls.csv',index=False)
# Method QC: technical repeat concordance from original calls, not population cohort.
raw=pd.read_csv(R.parent/'inputs/dose.txt',sep='\t');rep=[]
for a,b in [('S7','S7_Seq_Control_rep1'),('S7','S7_Seq_Control_rep3'),('S7_Seq_Control_rep1','S7_Seq_Control_rep3')]:
 ok=raw[a].notna()&raw[b].notna();rep.append({'a':a,'b':b,'n':int(ok.sum()),'concordance':float((raw.loc[ok,a]==raw.loc[ok,b]).mean())})
res['technical_replicates']=rep
(O/'analysis_results.json').write_text(json.dumps(res,indent=2,default=lambda x:x.item()));print(json.dumps(res,indent=2,default=lambda x:x.item()))
