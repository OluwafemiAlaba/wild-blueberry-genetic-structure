from pathlib import Path
import pandas as pd,numpy as np,openpyxl,re,json,hashlib
R=Path(__file__).resolve().parent;D=pd.read_csv(R.parent/'inputs/dose.txt',sep='\t'); records=[]
def parent(s):return re.sub(r'_(N|S|E|W|C)$','',s) if s.startswith('TP_') else s
w=openpyxl.load_workbook(R.parent/'inputs/amova/managed_unmanaged.xlsx',read_only=True,data_only=True)
mg={parent(r[0]):r[1].split('_')[-1] for r in w['Sheet1'].iter_rows(min_row=4,max_col=2,values_only=True)}
for f,species in [('angust_geographic','angustifolium'),('myrt_regions','myrtilloides')]:
 w=openpyxl.load_workbook(R.parent/'inputs/amova'/f'{f}.xlsx',read_only=True,data_only=True)
 for sample,site in w['Sheet1'].iter_rows(min_row=4,max_col=2,values_only=True):
  site={'Wyman':'Deblois','BBHF':'Jonesboro','FortKent':'Fort Kent'}.get(site,site)
  records.append([sample,parent(sample),species,site,mg.get(parent(sample),'Unspecified') if species=='angustifolium' else 'Unmanaged'])
for s in ['Bluejay','Nelson','MU6937','Bluegold','Patrioit','Bluecorp','Sunrise']:records.append([s,s,'highbush reference','Jonesboro','Reference'])
M=pd.DataFrame(records,columns=['sample','parent','species','site','management']);assert M.parent.is_unique;assert len(M)==172
miss=D[M['sample']].isna().mean();M['original_missing']=M['sample'].map(miss);M['retained']=M.original_missing<=.30;M.to_csv(R/'data/cohort_manifest.csv',index=False);M=M[M.retained].reset_index(drop=True)
# Exclude duplicated positions and retain biallelic SNPs.
unique=~D.duplicated(['CHROM','POS'],keep=False);snv=D.REF.isin(list('ACGT'))&D.ALT.isin(list('ACGT'));base=D[unique&snv].copy();X=base[M['sample']].to_numpy(float).T;p=np.nanmean(X,axis=0)/2;missing=np.mean(np.isnan(X),axis=0);keep=(missing<=.20)&(np.minimum(p,1-p)>=.05)&(np.nanstd(X,axis=0)>0);V=base.loc[keep,['SNP','CHROM','POS','REF','ALT']].reset_index(drop=True);X=X[:,keep];V['missing']=missing[keep]
# Prune LD within 500 kb; mean-fill missing calls only for LD screening.
means=np.nanmean(X,axis=0);Z=np.where(np.isnan(X),means,X);Z-=Z.mean(axis=0);sd=np.sqrt((Z*Z).sum(axis=0));Z/=sd;ret=[]
for j in range(X.shape[1]):
 prior=[i for i in ret if V.at[i,'CHROM']==V.at[j,'CHROM'] and abs(V.at[j,'POS']-V.at[i,'POS'])<=500000]
 if not prior or np.max((Z[:,prior].T@Z[:,j])**2)<=.2:ret.append(j)
V['LD_retained']=False;V.loc[ret,'LD_retained']=True;V.to_csv(R/'data/variant_manifest.csv',index=False);M.to_csv(R/'data/analysis_samples.csv',index=False);np.savez_compressed(R/'data/genotypes.npz',X=X,ld=np.array(ret),samples=M['sample'].values.astype(str))
np.savetxt(R/'data/all.lfmm',np.nan_to_num(X[:,ret],nan=9).astype(int),fmt='%d')
a=M.species.eq('angustifolium').to_numpy();np.savetxt(R/'data/angust.lfmm',np.nan_to_num(X[a][:,ret],nan=9).astype(int),fmt='%d');M[a].to_csv(R/'data/angust_samples.csv',index=False)
out={'input_sha256':hashlib.sha256((R.parent/'inputs/dose.txt').read_bytes()).hexdigest(),'source_records':len(D),'source_samples':len(D.columns)-5,'planned_n':172,'sample_QC_n':len(M),'excluded_samples':miss[miss>.30].to_dict(),'unique_biallelic_SNV_records':len(base),'filtered_variants':len(V),'ld_pruned_variants':len(ret),'missing_rate_filtered':float(np.isnan(X).mean()),'groups':M.groupby(['species','site','management']).size().to_dict()};out['groups']={str(k):int(v) for k,v in out['groups'].items()};(R/'results/preparation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
