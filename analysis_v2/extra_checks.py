from pathlib import Path
import numpy as np,pandas as pd,json
R=Path(__file__).resolve().parent;O=R/'results';rng=np.random.default_rng(20261014);out={}
for name,group in [('all','species'),('angust','site')]:
 m=pd.read_csv(O/f'{name}_PCA.csv');C=np.load(O/f'{name}_kernels.npz')['C'];labels=m[group].values;d=np.zeros(len(m));cats=list(dict.fromkeys(labels))
 for g in cats:
  ix=labels==g;d[ix]=np.linalg.norm(C[ix]-C[ix].mean(axis=0),axis=1)
 def f(y):
  total=((y-y.mean())**2).sum();ss=sum(sum(labels==g)*(y[labels==g].mean()-y.mean())**2 for g in cats);return (ss/(len(cats)-1))/((total-ss)/(len(y)-len(cats)))
 obs=f(d);hits=sum(f(rng.permutation(d))>=obs for _ in range(9999));out[name+'_dispersion']={'F':float(obs),'P':(hits+1)/10000}
a=json.loads((O/'analysis_results.json').read_text());t=[a['all']['species_test'],a['angust']['site_test'],a['myrt']['two_site_test'],a['angust']['management_adjusted_site']];order=sorted(range(4),key=lambda i:t[i]['P']);running=0
for j,i in enumerate(order):running=max(running,min(1,(4-j)*t[i]['P']));t[i]['P_Holm_primary']=running
(O/'analysis_results.json').write_text(json.dumps(a,indent=2));(O/'extra_checks.json').write_text(json.dumps(out,indent=2));print(out)
# Export imputed matrices with sample and marker identifiers.
v=pd.read_csv(R/'data/variant_manifest.csv').query('LD_retained');col=v.CHROM.astype(str)+':'+v.POS.astype(str)+':'+v.REF+':'+v.ALT
for name in ['all','angust']:
 m=pd.read_csv(O/f'{name}_PCA.csv');x=np.loadtxt(R/'data'/f'{name}_final.lfmm_imputed.lfmm',dtype=int);pd.DataFrame(x,index=m['sample'],columns=col).rename_axis('sample').to_csv(R/'data'/f'{name}_imputed_named.csv.gz',compression='gzip')
