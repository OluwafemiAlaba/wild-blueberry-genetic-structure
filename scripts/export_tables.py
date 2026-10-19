"""Export saved geographic and recalculated management results, without new inference."""
from pathlib import Path
import json,csv
import pandas as pd
import openpyxl
R=Path(__file__).resolve().parents[1];O=R/'outputs/tables';O.mkdir(parents=True,exist_ok=True)
main=[];full=[];counts=[]
for stem,taxon,panel in [('angust_geographic','V. angustifolium','A'),('myrt_without_hope','V. myrtilloides','B')]:
 w=openpyxl.load_workbook(R/'inputs/amova'/f'{stem}.xlsx',read_only=True,data_only=True)
 vals={}
 for row in w['PhiPT'].values:
  for j,v in enumerate(row):
   if isinstance(v,str) and v in ['Among Pops','Within Pops','Total','PhiPT']:vals[v]=row[j+1:j+6]
 n=int(w['PCoA']['B5'].value);phi,p=vals['PhiPT'][:2]
 for source,key in [('Among locations','Among Pops'),('Within locations','Within Pops'),('Total','Total')]:
  a=vals[key];row=dict(panel=panel,species=taxon,n=n,source=source,df=a[0],SS=a[1],MS=a[2],variance_component=a[3],variation_percent=100*a[4],PhiPT=phi,P=p)
  full.append(row)
  if source!='Total':main.append({k:row[k] for k in ['panel','species','n','source','variance_component','variation_percent','PhiPT','P']})
 groups=pd.Series([row[1] for row in w['Sheet1'].iter_rows(min_row=4,max_col=2,values_only=True)]).value_counts(sort=False)
 for group,number in groups.items():counts.append(dict(analysis=stem,species=taxon,source_group=group,n=int(number)))
 # Keep cell layout: geographical pairwise output upper/lower triangles have distinct meanings.
 with (O/f'Table_S{1 if panel=="A" else 2}_saved_pairwise_cells.csv').open('w') as f:csv.writer(f).writerows(w['PhiPTP'].values)
a=pd.read_csv(R/'outputs/management_Hope34_excluded/corrected_amova.csv');s=json.loads((R/'outputs/management_Hope34_excluded/recalculation_summary.json').read_text())
for t in a.itertuples():
 row=dict(panel='C',species='V. angustifolium',n=136,source=t.Source,df=t.df,SS=t.SS,MS=t.MS,variance_component=t.Variance_component,variation_percent=t.Percent,PhiPT=s['PhiPT'],P=s['P']);full.append(row)
 if t.Source!='Total':main.append({k:row[k] for k in ['panel','species','n','source','variance_component','variation_percent','PhiPT','P']})
for group,n in s['groups'].items():counts.append(dict(analysis='management_Hope34_excluded',species='V. angustifolium',source_group=group,n=n))
pd.DataFrame(main).to_csv(O/'Table1_AMOVA.csv',index=False)
pd.DataFrame(full).to_csv(O/'Table_S4_full_AMOVA.csv',index=False)
pd.DataFrame(counts).to_csv(O/'Table_S3_analysis_sample_counts.csv',index=False)
m=pd.read_csv(R/'deliverables/Management_differentiation_preview/pairwise_PhiPT_corrected.csv',index_col=0)
m.to_csv(O/'Table2_pairwise_PhiPT.csv')
(O/'README.txt').write_text('These exports support manuscript tables; they are not formatted manuscript replacements. Table S3 export contains verified analysis-specific counts, not unverified reported collection totals or GPS. PhiPTP exports preserve GenAlEx cell layout and labels; do not treat all numeric cells as PhiPT or as P values. Main Table 1 panels A-B and S1-S2 use saved geographic results; C uses recalculated management results.\n')
print('Exported six table-support files.')
