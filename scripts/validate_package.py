"""Check study dimensions, sample alignment and observed-call preservation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,openpyxl
R=Path(__file__).resolve().parents[1];checks={}
d=pd.read_csv(R/'inputs/dose.txt',sep='\t');assert d.shape==(7407,277);checks['source_shape']=[7407,272]
m=pd.read_csv(R/'analysis/data/analysis_samples.csv');v=pd.read_csv(R/'analysis/data/variant_manifest.csv');z=np.load(R/'analysis/data/genotypes.npz');assert len(m)==129 and len(v)==6230 and sum(v.LD_retained)==4640
assert np.array_equal(z['samples'],m['sample']);assert np.array_equal(np.flatnonzero(v.LD_retained),z['ld']);x=z['X'][:,z['ld']]
a=np.loadtxt(R/'analysis/data/all_final.lfmm_imputed.lfmm');assert a.shape==(129,4640) and set(np.unique(a))<={0,1,2};obs=np.isfinite(x);assert np.array_equal(a[obs],x[obs]);checks['imputed_observed_calls_preserved']=int(obs.sum());checks['imputed_shape']=list(a.shape)
q=pd.read_csv(R/'analysis/results/all_Q_final.csv').to_numpy();assert q.shape==(129,4) and np.allclose(q.sum(1),1,atol=1e-5)
ce=pd.read_csv(R/'analysis/results/all_cross_entropy.csv');assert ce.groupby('K').size().to_dict()=={k:5 for k in range(1,7)};checks['mean_cross_entropy']=ce.groupby('K').cross_entropy.mean().to_dict()
for stem,n in [('angust_geographic',135),('myrt_without_hope',28),('managed_unmanaged',137),('myrt_regions',30)]:
 w=openpyxl.load_workbook(R/'inputs/amova'/f'{stem}.xlsx',read_only=True,data_only=True);rows=list(w['Sheet1'].values);assert len(rows)-3==n and len(rows[2])-2==6491;assert all(len(row)==6493 for row in rows[2:]);checks[stem]={'n':n,'markers':6491}
g=pd.read_csv(R/'outputs/management_Hope34_excluded/corrected_management_input.csv');assert g.shape==(136,6493) and 'Hope34' not in g['sample'].tolist()
w=openpyxl.load_workbook(R/'inputs/amova/managed_unmanaged.xlsx',read_only=True,data_only=True);rows=list(w['Sheet1'].values)[3:];rows=[row for row in rows if row[0]!='Hope34'];assert g['sample'].tolist()==[r[0] for r in rows];assert np.array_equal(g.iloc[:,2:].to_numpy(),np.asarray([r[2:] for r in rows]));checks['management_retained_calls_unchanged']=True
h=pd.read_csv(R/'analysis/results/all_ward_order.csv');assert sorted(h['index'])==list(range(1,130));assert h['sample'].tolist()==m['sample'].iloc[h['index'].to_numpy()-1].tolist()
checks['status']='passed';checks['scope']='Stored data dimensions, alignment, observed-call preservation, historical cohort sizes and Hope34 exclusion. Not a raw-read or full stochastic workflow rerun.'
(R/'metadata/validation.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
