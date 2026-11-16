"""Independent reconstruction from Phase2 retained original RP/QV rows.
Only authorized boundary thresholds, CHN exclusion and weight concentration.
Does not import/run previous pipelines, fit models or produce new bootstrap CIs.
"""
from pathlib import Path
import csv,json,hashlib,math,shutil
from collections import defaultdict,Counter
REPO=Path(__file__).resolve().parents[2];HERE=REPO/'.reproduction/phase3_1';P2=REPO/'.reproduction/phase2';P3=REPO/'.reproduction/phase3'
for sub in ['SI_TABLES','metadata','checks','references','FIGURES_MAIN','FIGURES_SI']:(HERE/sub).mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
roots=[P2,P3]
before={str(p.relative_to(HERE.parent)):sha(p) for root in roots for p in root.rglob('*') if p.is_file()}
(HERE/'checks/FROZEN_PHASE1_3_SHA256.json').write_text(json.dumps(before,indent=2),encoding='utf-8')
checks=[]
def ck(name,a,b,tol=1e-10):
 ok=math.isclose(a,b,rel_tol=tol,abs_tol=tol) if isinstance(a,(int,float)) and isinstance(b,(int,float)) else a==b
 checks.append(dict(check=name,observed=a,expected=b,passed=ok))
 if not ok:raise AssertionError((name,a,b))
def write(rows,path):
 with (HERE/path).open('w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

variables=['p_kg','crop_usd2015','herbicide_kg','insecticide_kg','fungicide_kg']
original=P2/'data/original_source_rows.csv'
if not original.exists():original=HERE/'metadata/PHASE2_original_source_rows.csv'
annual=defaultdict(dict);rowkeys={}
with original.open(encoding='utf-8-sig',newline='') as f:
 for r in csv.DictReader(f):
  if r['variable'] not in variables or r['Value']=='':continue
  key=(r['iso3'],int(r['Year']));v=r['variable'];value=float(r['Value'])*1000
  ck('unit '+r['row_key'],r['Unit'],'1000 USD' if v=='crop_usd2015' else 't')
  if v=='crop_usd2015':assert r['Item Code']=='2041' and r['Element Code']=='58'
  else:assert r['Element Code']=='5157'
  assert v not in annual[key],(key,v)
  annual[key][v]=value;rowkeys[(key,v)]=r['row_key']
# Enumeration independently reconstructs all eligible common-state comparisons.
entities=sorted({k[0] for k in annual});periods={}
for iso in entities:
 for st in range(1990,2025,5):
  joint=[y for y in range(st,st+5) if 'p_kg' in annual[(iso,y)] and annual[(iso,y)].get('crop_usd2015',0)>0]
  common=[y for y in joint if all(v in annual[(iso,y)] for v in variables)]
  means={v:sum(annual[(iso,y)][v] for y in common)/len(common) if common else None for v in variables}
  periods[(iso,st)]=dict(joint=joint,common=common,means=means)
reconstructed=[];audit=[]
for iso in entities:
 for st in range(1995,2025,5):
  p0=periods[(iso,st-5)];p1=periods[(iso,st)]
  if min(len(p0['joint']),len(p1['joint']),len(p0['common']),len(p1['common']))<3:continue
  m0=p0['means'];m1=p1['means']
  if m0['p_kg']<=0:continue
  gp=m1['p_kg']/m0['p_kg']-1;gy=m1['crop_usd2015']/m0['crop_usd2015']-1
  if not(gp<0 and gy>0):continue
  inc={v:(m1[v]-m0[v])/m0['p_kg'] for v in variables[2:]};maxinc=max(0,*inc.values());id=f'{iso}_{st-5}_{st}_w5'
  r=dict(id=id,iso3=iso,previous_start=st-5,period_start=st,common_support0=';'.join(map(str,p0['common'])),common_support1=';'.join(map(str,p1['common'])),cp0=m0['p_kg'],cp1=m1['p_kg'],cy0=m0['crop_usd2015'],cy1=m1['crop_usd2015'],cgp=gp,cgy=gy,max_increment_over_P0=maxinc,threshold_defined_event=int(maxinc>.01))
  for v in variables[2:]:r[v+'0']=m0[v];r[v+'1']=m1[v];r[v+'_delta_over_P0']=inc[v]
  reconstructed.append(r)
  for side,p in [('0',p0),('1',p1)]:
   for y in p['common']:
    for v in variables:audit.append(dict(comparison_id=id,iso3=iso,side=side,year=y,variable=v,value_kg_or_real_USD=annual[(iso,y)][v],original_row_key=rowkeys[((iso,y),v)]))
ck('reconstructed common absolute count',len(reconstructed),189);ck('reconstructed reporting entities',len({r['iso3'] for r in reconstructed}),106)
frozen_table=P2/'tables/absolute_common_support_189.csv'
if not frozen_table.exists():frozen_table=HERE/'SI_TABLES/absolute_common_support_189.csv'
with frozen_table.open(encoding='utf-8-sig',newline='') as f:old={r['id']:r for r in csv.DictReader(f)}
ck('independently reconstructed sample keys',set(r['id'] for r in reconstructed),set(old))
for r in reconstructed:
 for v in ['cp0','cp1','cy0','cy1','cgp','cgy','max_increment_over_P0']+[v+s for v in variables[2:] for s in ['0','1']]:ck(r['id']+' '+v,r[v],float(old[r['id']][v]))
 for v in ['common_support0','common_support1']:ck(r['id']+' '+v,r[v],old[r['id']][v])
 ck(r['id']+' class event',r['threshold_defined_event'],int(old[r['id']]['substantial_gt1pctP0'] in ['True','true','1','1.0']))
write(reconstructed,'SI_TABLES/absolute_common_support_189_independently_reconstructed.csv');write(audit,'SI_TABLES/targeted_original_row_audit.csv')

ref=REPO/'data/metadata/USER_REFERENCE_boundary_sensitivity.csv'
if not ref.exists():ref=HERE/'metadata/USER_REFERENCE_boundary_sensitivity.csv'
reference=list(csv.DictReader(ref.open(encoding='utf-8-sig',newline='')));thresholds=[]
for t in [0,.005,.01,.02,.05]:
 eligible=[r for r in reconstructed if r['cgp'] < -t and r['cgy'] > t];k=sum(r['threshold_defined_event'] for r in eligible)
 row=dict(minimum_strict_change_fraction=t,eligible_n=len(eligible),class_gt_1pct_P0_k=k,conditional_rate=k/len(eligible),eligible_countries=len({r['iso3'] for r in eligible}),new_CI='not estimated; descriptive sensitivity only')
 rr=next(r for r in reference if math.isclose(float(r['minimum_strict_change_fraction']),t))
 for col in ['eligible_n','class_gt_1pct_P0_k','conditional_rate']:ck(f'reference match t={t} '+col,float(row[col]),float(rr[col]))
 thresholds.append(row)
write(thresholds,'SI_TABLES/TableS9_boundary_threshold_sensitivity.csv')
if ref.resolve()!=(HERE/'metadata/USER_REFERENCE_boundary_sensitivity.csv').resolve():shutil.copy2(ref,HERE/'metadata/USER_REFERENCE_boundary_sensitivity.csv')

weighted=[]
for weight in ['cy0','cp0']:
 total=sum(r[weight] for r in reconstructed);eventtotal=sum(r[weight]*r['threshold_defined_event'] for r in reconstructed)
 for exclude in ['', 'CHN']:
  d=[r for r in reconstructed if r['iso3']!=exclude];den=sum(r[weight] for r in d);num=sum(r[weight]*r['threshold_defined_event'] for r in d)
  weighted.append(dict(weight=weight,sample='all_common_absolute' if not exclude else 'exclude_CHN',n=len(d),countries=len({r['iso3'] for r in d}),events=sum(r['threshold_defined_event'] for r in d),weighted_event_numerator=num,weight_denominator=den,weighted_incidence=num/den,excluded_country='none' if not exclude else exclude,removed_weight_share_of_full_sample=1-den/total,new_CI='not estimated; full-sample CI remains frozen Phase2 only'))
for w,ex,e in [('cy0','all_common_absolute',.15630935391918524),('cy0','exclude_CHN',.311582),('cp0','all_common_absolute',.31307237832474316),('cp0','exclude_CHN',.379341)]:
 r=next(r for r in weighted if r['weight']==w and r['sample']==ex);ck(w+' '+ex,r['weighted_incidence'],e,tol=5e-7 if ex=='exclude_CHN' else 1e-10)
write(weighted,'SI_TABLES/TableS10_CHN_weight_exclusion_sensitivity.csv')
bycountry=[]
for iso in sorted({r['iso3'] for r in reconstructed}):
 d=[r for r in reconstructed if r['iso3']==iso]
 row=dict(iso3=iso,comparisons=len(d),events=sum(r['threshold_defined_event'] for r in d))
 for w in ['cy0','cp0']:
  row[w+'_sum']=sum(r[w] for r in d);row[w+'_share']=row[w+'_sum']/sum(r[w] for r in reconstructed);row[w+'_event_sum']=sum(r[w]*r['threshold_defined_event'] for r in d)
 bycountry.append(row)
bycountry.sort(key=lambda r:-r['cy0_share']);ck('CHN cumulative crop-value share',next(r['cy0_share'] for r in bycountry if r['iso3']=='CHN'),.498337,tol=5e-7)
write(bycountry,'SI_TABLES/TableS11_country_baseline_weight_concentration.csv')

# Compact checks exclude row-by-row unit records but retain all data comparisons.
write([r for r in checks if not r['check'].startswith('unit ')],'checks/INDEPENDENT_RECOMPUTATION_CHECKS.csv')
annual_snapshot=P2/'data/annual_reverified.csv'
if not annual_snapshot.exists():annual_snapshot=HERE/'SI_TABLES/annual_reverified.csv'
summary=dict(source_file=str(original.relative_to(HERE.parent)),source_sha256=sha(original),input_phase2_annual_file_sha256=sha(annual_snapshot),comparison_rows=len(reconstructed),countries=106,original_rows_audited=len(audit),checks_passed=len(checks),boundary_reference_all_rows_match=True,thresholds=thresholds,weighted_sensitivity=weighted,CHN_weight_details=next(r for r in bycountry if r['iso3']=='CHN'),new_models=0,new_bootstrap_runs=0)
(HERE/'checks/RECOMPUTATION_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
