"""Phase 2: independent read-only re-verification of Phase 1 inputs."""
import pathlib,json,hashlib,zipfile,argparse,sys,time,os
import pandas as pd
import numpy as np
REPOSITORY=pathlib.Path(__file__).resolve().parents[2]
HERE=pathlib.Path(os.environ.get('DECOUPLING_RUN_ROOT',str(REPOSITORY/'.reproduction/phase2')))
P1=pathlib.Path(os.environ.get('PHASE1_ROOT',str(HERE.parent/'phase1'))).resolve()
TABLE=HERE/'tables'
for d in ['tables','figures','metadata','logs','references','data']: (HERE/d).mkdir(parents=True,exist_ok=True)
MAJOR=['herbicide_kg','insecticide_kg','fungicide_kg']
ITEM={1357:'p_kg',1320:'herbicide_kg',1309:'insecticide_kg',1331:'fungicide_kg',1355:'other_kg',1341:'growthreg_kg',1345:'rodent_kg',1316:'mineral_legacy_kg',1354:'mineral_kg'}
def save(d,name,folder=TABLE):d.to_csv(folder/name,index=False,encoding='utf-8-sig')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def classify(gp,gy):
 if not np.isfinite(gp) or not np.isfinite(gy):return 'undetermined'
 if gy>0 and gp<0:return 'absolute'
 if gy>0 and 0<=gp<gy:return 'relative'
 return 'coupled_or_other'
def detailed(gp,gy):
 if not np.isfinite(gp) or not np.isfinite(gy):return 'undefined'
 if gy>0 and gp<0:return 'absolute'
 if gy>0 and gp<gy:return 'relative'
 if gy>0:return 'growth_coupling'
 if gy<0:return 'output_decline_p_decrease' if gp<0 else 'output_decline_p_nondecrease'
 return 'output_stagnant'
def load_reverified():
 archive=pathlib.Path(os.environ.get('PHASE1_ARCHIVE',str(HERE.parent/'pesticide-decoupling_phase1.zip'))).resolve();entries=[]
 with zipfile.ZipFile(archive) as z:
  for rel in ['data/raw/RP.zip','data/raw/QV.zip','data/processed/panel_country_year.csv','scripts/pipeline.py','metadata/run_metrics.json']:
   p=P1/rel;az='pesticide-decoupling/'+rel;ok=sha(p)==hashlib.sha256(z.read(az)).hexdigest();assert ok,rel
   entries.append(dict(input_file=rel,sha256=sha(p),matches_archive=ok,archive_sha256=sha(archive)))
 save(pd.DataFrame(entries),'phase1_input_snapshot.csv',HERE/'metadata')
 panel=pd.read_csv(P1/'data/processed/panel_country_year.csv');cw=pd.read_csv(P1/'metadata/country_crosswalk.csv');cw=cw.dropna(subset=['iso3'])
 raw=[]
 z=zipfile.ZipFile(P1/'data/raw/RP.zip');name=next(n for n in z.namelist() if 'Normalized' in n);rp=pd.read_csv(z.open(name),encoding='utf-8-sig')
 rp=rp[rp['Element Code'].eq(5157)&rp['Item Code'].isin(ITEM)&rp.Year.between(1990,2024)].copy();rp['variable']=rp['Item Code'].map(ITEM)
 rp=rp.merge(cw[['fao_area_code','iso3']],left_on='Area Code',right_on='fao_area_code',validate='many_to_one');assert rp.Unit.eq('t').all();rp['analysis_value']=rp.Value*1000
 rp['source_file']='Phase1/data/raw/RP.zip';rp['source_sha256']=sha(P1/'data/raw/RP.zip');raw.append(rp)
 z=zipfile.ZipFile(P1/'data/raw/QV.zip');name=next(n for n in z.namelist() if 'Normalized' in n);parts=[]
 for x in pd.read_csv(z.open(name),encoding='utf-8-sig',chunksize=200000,low_memory=False):parts.append(x[x['Item Code'].eq(2041)&x['Element Code'].eq(58)&x.Year.between(1990,2024)].copy())
 qv=pd.concat(parts,ignore_index=True).merge(cw[['fao_area_code','iso3']],left_on='Area Code',right_on='fao_area_code',validate='many_to_one');assert qv.Unit.eq('1000 USD').all();assert qv.Element.str.contains('constant 2014-2016').all();qv['variable']='crop_usd2015';qv['analysis_value']=qv.Value*1000;qv['source_file']='Phase1/data/raw/QV.zip';qv['source_sha256']=sha(P1/'data/raw/QV.zip');raw.append(qv)
 source=pd.concat(raw,ignore_index=True);source['row_key']=source['Area Code'].astype(str)+'|'+source.Year.astype(str)+'|'+source['Item Code'].astype(str)+'|'+source['Element Code'].astype(str)
 assert not source.duplicated(['iso3','Year','variable']).any()
 for var in list(ITEM.values())+['crop_usd2015']:
  d=source[source.variable.eq(var)][['iso3','Year','analysis_value']].rename(columns={'Year':'year'});chk=panel[['iso3','year',var]].merge(d,on=['iso3','year'],how='left',validate='one_to_one');assert np.allclose(chk[var],chk.analysis_value,equal_nan=True),var
 save(source,'original_source_rows.csv',HERE/'data')
 # Reconstruct key measured values from original rows; covariates reused as authorized.
 wide=source.pivot(index=['iso3','Year'],columns='variable',values='analysis_value').reset_index().rename(columns={'Year':'year'})
 for v in list(ITEM.values())+['crop_usd2015']:
  if v not in wide:wide[v]=np.nan
 panel=panel.drop(columns=list(ITEM.values())+['crop_usd2015']).merge(wide,on=['iso3','year'],how='left',validate='one_to_one')
 quality=pd.read_csv(P1/'results/tables/country_quality.csv');panel=panel.merge(quality,on='iso3',how='left',validate='many_to_one')
 return panel,source

def comparisons(panel,width=5,minobs=3,start=1990,yvar='crop_usd2015',allowed_flags=None,endpoint=False):
 notes=pd.read_csv(P1/'metadata/country_notes_coded_changes.csv');proxy=pd.read_csv(P1/'metadata/country_proxy_notes.csv').set_index('iso3').proxy_or_measurement_note.to_dict();out=[]
 starts=list(range(start,2025-width+1,width))
 for iso,g in panel.groupby('iso3',sort=True):
  gr=g.set_index('year');allperiods={}
  for st in starts:
   d=g[g.year.between(st,st+width-1)&g.p_kg.notna()&g[yvar].gt(0)].copy()
   if allowed_flags is not None:d=d[d.p_kg_flag.isin(allowed_flags)]
   allperiods[st]=d
  for st in starts[1:]:
   aa=allperiods[st-width];bb=allperiods[st];n0=len(aa);n1=len(bb);valid=n0>=minobs and n1>=minobs
   if endpoint:
    aa=aa[aa.year.eq(st-1)];bb=bb[bb.year.eq(st+width-1)];valid=len(aa)==1 and len(bb)==1
   r=dict(iso3=iso,id=f'{iso}_{st-width}_{st}_w{width}',width=width,previous_start=st-width,previous_end=st-1,period_start=st,period_end=st+width-1,n_previous=n0,n_current=n1,region=g.region.iloc[0].strip(),income=g.income.iloc[0].strip(),fao_area=g.fao_area.iloc[0],fao_area_code=int(g.fao_area_code.iloc[0]),national_A_share=g.official_share.iloc[0],national_I_share=g.imputed_share.iloc[0],proxy_note=proxy.get(iso,'not_coded_no_claim_of_direct_use'),p_source='RP1357/5157 t x1000',y_source=f'QV2041/58 constant2014-2016 USD' if yvar=='crop_usd2015' else yvar,source_valid=bool(valid))
   changes=notes[notes.iso3.eq(iso)&notes.change_year.gt(st-width)&notes.change_year.le(st+width-1)];r['known_break']=len(changes)>0;r['known_break_years']=';'.join(changes.change_year.astype(str));r['questionnaire2023_overlap']=st+width-1>=2023 and st-width<2023
   for suffix,d in [('0',aa),('1',bb)]:
    r['support_years'+suffix]=';'.join(d.year.astype(str))
    for v in ['p_kg',yvar,'cropland_ha','gdp_pc_const','cereal_share','vegfruit_share']:
     r[(('y' if v==yvar else v)+suffix)]=d[v].mean() if len(d) else np.nan
     r[v+'_n'+suffix]=int(d[v].notna().sum())
    for flag in ['A','E','I','X']:r['p_'+flag+'_share'+suffix]=d.p_kg_flag.eq(flag).mean() if len(d) else np.nan
   p0=r['p_kg0'];y0=r['y0'];r['gp']=r['p_kg1']/p0-1 if valid and p0>0 else np.nan;r['gy']=r['y1']/y0-1 if valid and y0>0 else np.nan
   r['state']=classify(r['gp'],r['gy']);r['detail_state']=detailed(r['gp'],r['gy']);r['near_zero_y_growth']=bool(abs(r['gy'])<.01) if np.isfinite(r['gy']) else None;r['elasticity']=r['gp']/r['gy'] if np.isfinite(r['gy']) and abs(r['gy'])>=.01 else np.nan
   a=aa[aa[MAJOR].notna().all(axis=1)];b=bb[bb[MAJOR].notna().all(axis=1)];need=1 if endpoint else minobs;complete=valid and len(a)>=need and len(b)>=need and a.p_kg.mean()>0
   r['common_complete']=bool(complete);r['common_n0']=len(a);r['common_n1']=len(b);r['common_support0']=';'.join(a.year.astype(str));r['common_support1']=';'.join(b.year.astype(str))
   r['cp0']=a.p_kg.mean() if complete else np.nan;r['cp1']=b.p_kg.mean() if complete else np.nan;r['cy0']=a[yvar].mean() if complete else np.nan;r['cy1']=b[yvar].mean() if complete else np.nan
   r['cgp']=r['cp1']/r['cp0']-1 if complete else np.nan;r['cgy']=r['cy1']/r['cy0']-1 if complete else np.nan;r['common_state']=classify(r['cgp'],r['cgy'])
   r['common_known_negative_residual']=bool(((a.p_kg-a[MAJOR].sum(axis=1)) < -2000).any() or ((b.p_kg-b[MAJOR].sum(axis=1)) < -2000).any()) if complete else None
   for c in MAJOR:
    va=aa[aa[c].notna()];vb=bb[bb[c].notna()];ok=valid and len(va)>=need and len(vb)>=need;r[c+'_complete']=bool(ok)
    for suffix,d in [('0',va),('1',vb)]:r[c+'_own'+suffix]=d[c].mean() if ok else np.nan
    for suffix,d in [('0',a),('1',b)]:r[c+suffix]=d[c].mean() if complete else np.nan
    v0=r[c+'0'];v1=r[c+'1'];r[c+'_delta_t']=(v1-v0)/1000;r[c+'_growth']=v1/v0-1 if v0>0 else np.nan;r[c+'_entry']=bool(v0==0 and v1>0) if complete else None
    r[c+'_share_pp']=(v1/r['cp1']-v0/r['cp0'])*100 if complete and r['cp1']>0 else np.nan;r[c+'_delta_over_P0']=(v1-v0)/r['cp0'] if complete else np.nan
   if complete:
    delta=np.array([r[c+'1']-r[c+'0'] for c in MAJOR]);growth=np.array([r[c+'_growth'] for c in MAJOR]);zeros=[bool(r[c+'_entry']) for c in MAJOR]
    r['any_increase']=bool((delta>0).any());r['positive_class_baselines']=all(r[c+'0']>0 for c in MAJOR);r['any_zero_entry']=any(zeros)
    r['any_gt5']=bool((growth>.05).any()) if r['positive_class_baselines'] else np.nan;r['any_gt10']=bool((growth>.10).any()) if r['positive_class_baselines'] else np.nan
    r['positive_increment_over_P0']=sum(np.maximum(delta,0))/r['cp0'];r['max_increment_over_P0']=max(np.maximum(delta,0))/r['cp0'];r['substantial_gt1pctP0']=r['max_increment_over_P0']>.01;r['substantial_gt5pctP0']=r['max_increment_over_P0']>.05
    r['delta_total_t']=(r['cp1']-r['cp0'])/1000;r['delta_HIF_t']=delta.sum()/1000;r['residual0_t']=(r['cp0']-sum(r[c+'0'] for c in MAJOR))/1000;r['residual1_t']=(r['cp1']-sum(r[c+'1'] for c in MAJOR))/1000;r['delta_residual_t']=r['residual1_t']-r['residual0_t'];r['class_coverage0']=sum(r[c+'0'] for c in MAJOR)/r['cp0'];r['class_coverage1']=sum(r[c+'1'] for c in MAJOR)/r['cp1'] if r['cp1']>0 else np.nan
   else:
    for v in ['any_increase','positive_class_baselines','any_zero_entry','any_gt5','any_gt10','positive_increment_over_P0','max_increment_over_P0','substantial_gt1pctP0','substantial_gt5pctP0','delta_total_t','delta_HIF_t','residual0_t','residual1_t','delta_residual_t','class_coverage0','class_coverage1']:r[v]=np.nan
   out.append(r)
 return pd.DataFrame(out)

def main():
 print('Phase2 input verification; no Phase1 writes/downloads',flush=True);panel,source=load_reverified();save(panel,'annual_reverified.csv',HERE/'data');grid=comparisons(panel);save(grid,'period_grid_including_undetermined.csv');b=grid[grid.state.ne('undetermined')].copy();save(b,'country_period_all.csv')
 old=pd.read_csv(P1/'results/tables/decoupling_5yr_country.csv');j=b.merge(old,left_on=['iso3','period_start'],right_on=['iso3','period_start'],validate='one_to_one');assert len(j)==len(old)==893
 for left,right in [('gp','g_p'),('gy','g_y'),('p_kg0','p_previous'),('p_kg1','p_current'),('y0','y_previous'),('y1','y_current')]:assert np.allclose(j[left],j[right]),left
 assert (j.detail_state==j.classification).all();absb=b[b.state.eq('absolute')];common=b[b.common_state.eq('absolute')];assert len(absb)==200 and len(b[b.state.eq('relative')])==133 and len(common)==189 and int(common.any_increase.sum())==102
 ids=set(common.id);lost=absb[~absb.id.isin(ids)].copy();entered=common[~common.id.isin(absb.id)].copy()
 lost['loss_reason']=np.where(~lost.common_complete,'insufficient_HIF_common_years','changed_total_state_on_common_years');save(lost,'absolute_200_to_189_exclusions.csv');save(entered,'absolute_common_support_entrants.csv');save(common,'absolute_common_support_189.csv');save(absb,'absolute_all_200.csv')
 # Per-original-observation provenance for all reported absolute-total comparisons.
 rows=[]
 for _,r in absb.iterrows():
  d=source[source.iso3.eq(r.iso3)&source.Year.between(r.previous_start,r.period_end)].copy();d['comparison_id']=r.id;d['side']=np.where(d.Year<r.period_start,'previous','current');d['in_joint_py_support']=d.Year.isin([int(y) for y in (r.support_years0+';'+r.support_years1).split(';') if y]);d['in_common_HIF_support']=d.Year.isin([int(y) for y in (r.common_support0+';'+r.common_support1).split(';') if y]);rows.append(d)
 audit=pd.concat(rows,ignore_index=True);save(audit,'absolute_comparisons_original_row_audit.csv')
 long=[]
 for _,r in absb.iterrows():
  for c in MAJOR:
   rr=dict(comparison_id=r.id,iso3=r.iso3,previous_start=r.previous_start,previous_end=r.previous_end,period_start=r.period_start,period_end=r.period_end,category=c,own_support_complete=r[c+'_complete'],common_HIF_complete=r.common_complete,common_total_state=r.common_state,own_start_t=r[c+'_own0']/1000,own_end_t=r[c+'_own1']/1000,start_t=r[c+'0']/1000,end_t=r[c+'1']/1000,delta_t=r[c+'_delta_t'],growth_pct=r[c+'_growth']*100,share_change_pp=r[c+'_share_pp'],increment_pct_P0=r[c+'_delta_over_P0']*100,zero_entry=r[c+'_entry'],known_break=r.known_break,record_source='RP/'+str({v:k for k,v in ITEM.items()}[c])+'/5157; original_row_audit.csv')
   rr['category_status']='common_years_supported_semantic_audit_partial' if r.common_complete else ('own_class_only_common_HIF_insufficient' if r[c+'_complete'] else 'not_classifiable');long.append(rr)
 save(pd.DataFrame(long),'absolute_category_changes_long.csv')
 counts=dict(country_years=int(panel.joint_py.sum()),countries=panel[panel.joint_py].iso3.nunique(),comparisons=len(b),absolute=len(absb),relative=int(b.state.eq('relative').sum()),common_absolute=len(common),any_increase=int(common.any_increase.sum()),exclusions=len(lost),entrants=len(entered));(HERE/'metadata/reproduction_counts.json').write_text(json.dumps(counts,indent=2),encoding='utf-8');print(json.dumps(counts,indent=2),flush=True)
if __name__=='__main__':main()
