"""Independent raw-extract arithmetic and bootstrap integrity audit."""
from pathlib import Path
from collections import defaultdict
from decimal import Decimal,localcontext
import csv,json,hashlib,math
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[1];BASE=R/'00_BASELINE/repository' if (R/'00_BASELINE/repository').exists() else R.parents[1];checks=[]
def ck(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
annual=defaultdict(dict)
with (BASE/'data/processed/original_source_rows.csv').open(encoding='utf-8-sig',newline='') as f:
    for r in csv.DictReader(f):
        if r['variable'] in ['p_kg','crop_usd2015','herbicide_kg','insecticide_kg','fungicide_kg'] and r['Value']:
            annual[(r['iso3'],int(r['Year']))][r['variable']]=Decimal(r['Value'])*1000
cross=pd.read_csv(BASE/'data/metadata/country_crosswalk.csv').dropna(subset=['iso3']);mapiso=cross.set_index('fao_area_code').iso3.to_dict();land=defaultdict(dict)
with (BASE/'data/source_extracts/RL_selected_original_rows.csv').open(encoding='utf-8-sig',newline='') as f:
    for x in csv.DictReader(f):
        if int(x['Area Code']) not in mapiso or not x['Value'] or x['Item Code'] not in ['6621','6650']:continue
        ck('land unit '+x['Area Code']+' '+x['Year']+' '+x['Item Code'],x['Unit']=='1000 ha')
        land[(mapiso[int(x['Area Code'])],int(x['Year']))][x['Item Code']]=Decimal(x['Value'])*1000
for k,v in land.items():
    if len(v)==2 and all(x>=0 for x in v.values()):annual[k]['cropland_ha']=sum(v.values())
records=pd.read_csv(R/'01_ANALYSIS/comparison_records.csv');actual_events=[];available=0
for _,r in records.iterrows():
    cols=['p_kg','crop_usd2015','cropland_ha']+(['herbicide_kg','insecticide_kg','fungicide_kg'] if r.support=='PYAHIF' else [])
    means=[];good=True
    for side,start in [('0',int(r.previous_start)),('1',int(r.period_start))]:
        use=[y for y in range(start,start+5) if all(v in annual[(r.iso3,y)] for v in cols) and annual[(r.iso3,y)]['crop_usd2015']>0 and annual[(r.iso3,y)]['cropland_ha']>0 and annual[(r.iso3,y)]['p_kg']>=0]
        good=good and len(use)>=3
        ck('years '+r.id+' '+r.support+' '+side,';'.join(map(str,use))==str(r['years'+side]) if use else pd.isna(r['years'+side]))
        with localcontext() as ctx:
            ctx.prec=40;m={v:sum(annual[(r.iso3,y)][v] for y in use)/len(use) for v in cols} if use else {};means.append(m)
    good=good and means[0].get('p_kg',0)>0;ck('support '+r.id+' '+r.support,good==r.support_available)
    if not good:continue
    m0,m1=means;I0=m0['p_kg']/m0['cropland_ha'];I1=m1['p_kg']/m1['cropland_ha'];event=m1['p_kg']<m0['p_kg'] and m1['crop_usd2015']>m0['crop_usd2015'] and I1>I0
    ck('intensity identity '+r.id+' '+r.support,math.isclose(float(I1/I0-1),r.gI,rel_tol=1e-9,abs_tol=1e-12))
    ck('event '+r.id+' '+r.support,bool(event)==bool(r.joint_divergence))
    if r.support=='PYA':
        available+=1
        if event and r.frozen_state=='absolute':actual_events.append(r.id);ck('denominator arithmetic '+r.id,r.ga<r.gp<0 and r.gy>0)
    else:
        for c in cols[3:]:ck('class intensity '+r.id+' '+c,bool(m1[c]/m1['cropland_ha']>m0[c]/m0['cropland_ha'])==bool(r[c+'_intensity_up']))
ck('land-supported887',available==887);ev=pd.read_csv(R/'02_TABLES/joint_divergence_events.csv');ck('all18eventkeys',set(actual_events)==set(ev.id) and len(ev)==18)
est=pd.read_csv(R/'02_TABLES/intensity_cluster_estimates.csv')
for _,r in est.iterrows():
    ck('ratio '+r.cohort+' '+r.event,math.isclose(r.events/r.n,r.estimate,abs_tol=1e-12) if r.n else pd.isna(r.estimate))
    if r.countries>=20:ck('CI bounds '+r.cohort+' '+r.event,0<=r.ci_low<=r.ci_high<=1)
# Verify principal CI via direct country resampling,matching the archived multinomial draw exactly.
sel=records[records.support.eq('PYA')&records.support_available&records.frozen_state.eq('absolute')];g=sel.groupby('iso3').joint_divergence.agg(['sum','size']);row=est[est.cohort.eq('frozen_absolute_joint')].iloc[0];rng=np.random.default_rng(int(row.seed));cnt=rng.multinomial(len(g),np.ones(len(g))/len(g),size=4000);rep=(cnt@g['sum'].to_numpy())/(cnt@g['size'].to_numpy());lo,hi=np.quantile(rep,[.025,.975]);ck('principal CI exact',np.allclose([lo,hi],[row.ci_low,row.ci_high],rtol=0,atol=1e-14))
pd.DataFrame(checks).to_csv(R/'07_QA/INDEPENDENT_ARITHMETIC_CHECKS.csv',index=False,encoding='utf-8-sig')
result=dict(checks=len(checks),passed=all(x['passed'] for x in checks),raw_land_unit='1000ha x1000',independent_decimal_original_rows=True,events=18,denominator=199,clusters=108,new_ci=[float(lo),float(hi)],minimum_positive_gI=float(ev.gI.min()),median_positive_gI=float(ev.gI.median()),excluded=records[records.support.eq('PYA')&~records.support_available][['id','iso3','frozen_state']].to_dict('records'))
(R/'07_QA/INDEPENDENT_VALIDATION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
