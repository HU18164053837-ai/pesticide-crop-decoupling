"""Reproducible Phase 1 analysis. Run with --root PATH; no interpolation."""
import argparse,pathlib,json,zipfile,hashlib,datetime,sys,urllib.request,platform
import numpy as np
import pandas as pd
import pycountry
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
ap=argparse.ArgumentParser();ap.add_argument('--root',default=str(pathlib.Path(__file__).resolve().parents[1]));ap.add_argument('--download',action='store_true');args=ap.parse_args()
R=pathlib.Path(args.root);RAW=R/'data/raw';T=R/'results/tables';F=R/'results/figures';M=R/'metadata'
for p in [RAW,T,F,M,R/'data/processed',R/'reports']:p.mkdir(parents=True,exist_ok=True)
def save(d,name,folder=T):d.to_csv(folder/name,index=False,encoding='utf-8-sig')
def readjson(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def fetch(url,path):
 if not path.exists():
  with urllib.request.urlopen(url,timeout=180) as resp:path.write_bytes(resp.read())
if args.download:
 fetch('https://bulks-faostat.fao.org/production/datasets_E.json',RAW/'datasets_E.json')
 cat=readjson(RAW/'datasets_E.json')['Datasets']['Dataset']
 for x in cat:
  if x['DatasetCode'] in ['RP','QV','RL','QI']:fetch(x['FileLocation'],RAW/(x['DatasetCode']+'.zip'))
 fetch('https://api.worldbank.org/v2/country?format=json&per_page=400',RAW/'wb_countries.json')
 for c in ['NY.GDP.PCAP.KD','NV.AGR.TOTL.KD','AG.PRD.CROP.XD']:
  fetch(f'https://api.worldbank.org/v2/country/all/indicator/{c}?format=json&date=1990:2024&per_page=20000',RAW/f'wb_{c}.json')
cat={x['DatasetCode']:x for x in readjson(RAW/'datasets_E.json')['Datasets']['Dataset']}
register=[];retrieval=[];selected={};areas=[]
filters={'RP':'Element Code=5157; Item Code in 1357,1320,1309,1331,1316,1354,1355,1341,1345', 'QV':'Element Code=58; Item Code in 2041,1717,1730,1739', 'RL':'Element Code=5110; Item Code in 6621,6650','QI':'Element Code=432; Item Code=2041'}
items={'RP':[1357,1320,1309,1331,1316,1354,1355,1341,1345],'QV':[2041,1717,1730,1739],'RL':[6621,6650],'QI':[2041]};el={'RP':5157,'QV':58,'RL':5110,'QI':432}
for code in items:
 p=RAW/(code+'.zip');z=zipfile.ZipFile(p)
 for n in z.namelist():
  if 'Normalized' not in n:
   d=pd.read_csv(z.open(n),encoding='utf-8-sig');save(d,code+'_'+n,M)
   if 'AreaCodes' in n:areas.append(d.assign(domain=code))
 parts=[];nraw=0
 n=next(n for n in z.namelist() if 'Normalized' in n)
 for d in pd.read_csv(z.open(n),encoding='utf-8-sig',chunksize=180000,low_memory=False):
  nraw+=len(d);s=d[d['Item Code'].isin(items[code]) & d['Element Code'].eq(el[code]) & d.Year.between(1990,2024)].copy();parts.append(s)
 d=pd.concat(parts,ignore_index=True);selected[code]=d;save(d,code+'_selected_original_rows.csv',R/'data/processed')
 retrieval.append(dict(source=code,catalog_rows=cat[code]['FileRows'],retrieved_rows=nraw,selected_rows=len(d),complete=nraw==int(cat[code]['FileRows'])))
 register.append(dict(variable=code,source='FAOSTAT',url=cat[code]['FileLocation'],units='see selected_units.csv',code=code,filter=filters[code]+'; years1990:2024',extraction_date='2026-11-06',source_update=cat[code]['DateUpdate'],raw_file='data/raw/'+p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 print('read',code,nraw,len(d),flush=True)
save(pd.DataFrame(retrieval),'retrieval_reconciliation.csv',M)
save(pd.concat([d[['Item Code','Item','Element Code','Element','Unit']].drop_duplicates().assign(domain=c) for c,d in selected.items()]),'selected_units.csv',M)
# Explicit M49 to ISO crosswalk; historical territories never bridged.
a=pd.concat(areas).drop_duplicates(['Area Code','M49 Code','Area']);cw=[]
for _,r in a.iterrows():
 m=str(r['M49 Code']).lstrip("'").zfill(3);c=pycountry.countries.get(numeric=m);status='matched_iso' if c else 'historical_or_non_iso'
 if int(r['Area Code'])>=5000:status='aggregate';c=None
 if m=='159':status='China_aggregate_excluded';c=None
 cw.append({'fao_area_code':r['Area Code'],'m49':m,'fao_area':r['Area'],'iso3':c.alpha_3 if c else None,'status':status})
cross=pd.DataFrame(cw).drop_duplicates();save(cross,'country_crosswalk_name_variants.csv',M)
cross=cross.drop_duplicates(['fao_area_code','m49','iso3','status']);save(cross,'country_crosswalk.csv',M)
assert not cross.dropna(subset=['iso3']).duplicated('iso3').any(),'ambiguous crosswalk'
wb=readjson(RAW/'wb_countries.json');assert len(wb[1])==int(wb[0]['total'])
countries=pd.DataFrame([{'iso3':x['id'],'wb_name':x['name'],'region':x['region']['value'],'income':x['incomeLevel']['value']} for x in wb[1] if x['region']['value']!='Aggregates']);save(countries,'wb_country_metadata.csv',M)
iso=cross.iso3.dropna().unique();panel=pd.MultiIndex.from_product([iso,range(1990,2025)],names=['iso3','year']).to_frame(index=False).merge(cross.dropna(subset=['iso3']),on='iso3',validate='many_to_one').merge(countries,on='iso3',how='left',validate='many_to_one')
panel['region']=panel.region.fillna('Not classified in WDI');panel['income']=panel.income.fillna('Not classified in WDI')
joins=[];invalid=[];unitlog=[]
varmap={'RP':{1357:'p_kg',1320:'herbicide_kg',1309:'insecticide_kg',1331:'fungicide_kg',1316:'mineral_legacy_kg',1354:'mineral_kg',1355:'other_kg',1341:'growthreg_kg',1345:'rodent_kg'},'QV':{2041:'crop_usd2015',1717:'cereal_usd2015',1730:'oil_usd2015',1739:'vegfruit_usd2015'},'RL':{6621:'arable_ha',6650:'permanent_ha'},'QI':{2041:'qi_crop'}}
for code,d in selected.items():
 dd=d.merge(cross[['fao_area_code','iso3','status']],left_on='Area Code',right_on='fao_area_code',how='left',validate='many_to_one')
 joins.append(dict(source=code,selected_rows=len(d),matched_rows=int(dd.iso3.notna().sum()),excluded_rows=int(dd.iso3.isna().sum())))
 save(dd[dd.iso3.isna()].groupby(['Area Code','Area','status'],dropna=False).size().reset_index(name='rows'),code+'_unmatched_records.csv',M)
 for item,var in varmap[code].items():
  v=dd[dd['Item Code'].eq(item)&dd.iso3.notna()].copy();assert not v.duplicated(['iso3','Year']).any(),var
  orig=v.Value.copy();v.loc[v.Value<0,'Value']=np.nan
  invalid.append(dict(variable=var,negative=int((orig<0).sum()),zero=int(orig.eq(0).sum()),missing=int(orig.isna().sum())))
  factor=1000 if code in ['RP','QV','RL'] else 1
  expected={'RP':'t','QV':'1000 USD','RL':'1000 ha'}.get(code)
  if expected is not None:assert set(v.Unit.dropna())<={expected},(var,set(v.Unit))
  v[var]=v.Value*factor
  v=v.rename(columns={'Year':'year','Flag':var+'_flag'})
  panel=panel.merge(v[['iso3','year',var,var+'_flag']],on=['iso3','year'],how='left',validate='one_to_one')
for code,var in [('NY.GDP.PCAP.KD','gdp_pc_const'),('NV.AGR.TOTL.KD','ag_value_added_const'),('AG.PRD.CROP.XD','wb_crop_index')]:
 p=RAW/f'wb_{code}.json';j=readjson(p);assert j[0]['pages']==1 and len(j[1])==int(j[0]['total'])
 d=pd.DataFrame([dict(iso3=x['countryiso3code'],year=int(x['date']),**{var:x['value']}) for x in j[1]])
 valid=d.iso3.isin(countries.iso3);joins.append(dict(source=code,selected_rows=len(d),matched_rows=int(valid.sum()),excluded_rows=int((~valid).sum())))
 panel=panel.merge(d[valid],on=['iso3','year'],how='left',validate='one_to_one')
 register.append(dict(variable=var,source='World Bank WDI',url=f'https://api.worldbank.org/v2/country/all/indicator/{code}?format=json&date=1990:2024&per_page=20000',units='constant 2015 US$' if code!='AG.PRD.CROP.XD' else '2014-2016=100',code=code,filter='non aggregate country; 1990:2024; null preserved',extraction_date='2026-11-06',source_update=j[0]['lastupdated'],raw_file='data/raw/'+p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
save(pd.DataFrame(register),'source_register.csv',M);save(pd.DataFrame(joins),'join_log.csv',M);save(pd.DataFrame(invalid),'invalid_values.csv',M)
panel['cropland_ha']=panel[['arable_ha','permanent_ha']].sum(axis=1,min_count=2)
panel['p_kg_ha']=panel.p_kg/panel.cropland_ha.where(panel.cropland_ha>0)
panel['p_kg_per_1000usd']=panel.p_kg/(panel.crop_usd2015.where(panel.crop_usd2015>0)/1000)
for c in ['cereal','oil','vegfruit']:panel[c+'_share']=panel[c+'_usd2015']/panel.crop_usd2015.where(panel.crop_usd2015>0)
panel['joint_py']=panel.p_kg.notna()&panel.crop_usd2015.gt(0)
major=['herbicide_kg','insecticide_kg','fungicide_kg']
panel['major_complete']=panel[major].notna().all(axis=1)
panel['major_sum']=panel[major].sum(axis=1,min_count=3)
panel['residual_kg']=panel.p_kg-panel.major_sum
save(panel,'panel_country_year.csv',R/'data/processed')
save(panel[['iso3','fao_area','region','income','year','p_kg_flag','joint_py','major_complete']].assign(**{c+'_available':panel[c].notna() for c in ['p_kg','crop_usd2015','cropland_ha','gdp_pc_const','qi_crop','wb_crop_index']+major}),'country_year_coverage.csv')
coverage=[]
for v in ['p_kg','crop_usd2015','cropland_ha','gdp_pc_const','qi_crop','wb_crop_index']+major:
 d=panel[panel[v].notna()];coverage.append(dict(variable=v,countries=d.iso3.nunique(),country_years=len(d),first_year=d.year.min(),last_year=d.year.max()))
save(pd.DataFrame(coverage),'variable_coverage.csv')
flag=panel.groupby(['region','p_kg_flag'],dropna=False).size().reset_index(name='country_years');save(flag,'pesticide_flags_by_region.csv')
quality=panel[panel.p_kg.notna()].assign(official=lambda d:d.p_kg_flag.eq('A'),imputed=lambda d:d.p_kg_flag.eq('I')).groupby('iso3').agg(p_years=('p_kg','size'),official_share=('official','mean'),imputed_share=('imputed','mean')).reset_index();save(quality,'country_quality.csv')
panel=panel.merge(quality,on='iso3',how='left',validate='many_to_one')
sample=[];runs=[]
for iso,g in panel.groupby('iso3'):
 ys=g[g.joint_py].year.sort_values().tolist();best=cur=0;last=None
 for y in ys:cur=cur+1 if last is not None and y==last+1 else 1;best=max(best,cur);last=y
 runs.append(dict(iso3=iso,max_consecutive_joint_years=best,joint_years=len(ys)))
run=pd.DataFrame(runs);save(run,'country_consecutive_coverage.csv')
for start in [1990,2000,2010]:
 d=panel[panel.year.ge(start)&panel.joint_py];cnt=d.groupby('iso3').size()
 sample.append(dict(case=f'{start}-2024_any',countries=len(cnt),country_years=len(d)))
 ids=cnt[cnt.eq(2025-start)].index;sample.append(dict(case=f'{start}-2024_complete',countries=len(ids),country_years=len(d[d.iso3.isin(ids)])))
ids=run[run.max_consecutive_joint_years.ge(10)].iso3;sample.append(dict(case='at_least_10_consecutive_joint_years',countries=len(ids),country_years=int(panel[panel.iso3.isin(ids)].joint_py.sum())))
save(pd.DataFrame(sample),'candidate_samples.csv')
# Explicit mutually exclusive quadrants. Growth fractions, not percentages.
def classify(gp,gy):
 if not np.isfinite(gp) or not np.isfinite(gy):return 'undefined'
 if gy>0:
  if gp<0:return 'absolute'
  if gp<gy:return 'relative'
  return 'growth_coupling'
 if gy<0:return 'output_decline_p_decrease' if gp<0 else 'output_decline_p_nondecrease'
 return 'output_stagnant'
def comparisons(data,yvar='crop_usd2015',width=5,minobs=3,start=1990,rolling=False):
 x=data[data.year.ge(start)&data.p_kg.notna()&data[yvar].gt(0)].copy();res=[]
 for iso,g in x.groupby('iso3'):
  periods={}
  for b in (range(start,2025-width+1) if rolling else range(start,2025,width)):
   if b+width>2025:continue
   gg=g[g.year.between(b,b+width-1)]
   if len(gg)>=minobs:periods[b]=gg
  for b,gg in periods.items():
   prev=periods.get(b-width)
   if prev is None:continue
   p0=prev.p_kg.mean();p1=gg.p_kg.mean();y0=prev[yvar].mean();y1=gg[yvar].mean()
   if p0<=0 or y0<=0:continue
   gp=p1/p0-1;gy=y1/y0-1;cl=classify(gp,gy)
   row=dict(iso3=iso,previous_start=b-width,period_start=b,period=f'{b-width}-{b-1} to {b}-{b+width-1}',n_previous=len(prev),n_current=len(gg),p_previous=p0,p_current=p1,y_previous=y0,y_current=y1,g_p=gp,g_y=gy,classification=cl,decoupled=cl in ['absolute','relative'],tapio_unstable=abs(gy)<.01,tapio_elasticity=gp/gy if abs(gy)>=.01 else np.nan)
   for c in major:
    aa=prev[prev[c].notna()];bb=gg[gg[c].notna()];ok=len(aa)>=minobs and len(bb)>=minobs
    row[c+'_n_previous']=len(aa);row[c+'_n_current']=len(bb)
    row[c+'_direction']='unclassified' if not ok else ('increase' if bb[c].mean()>aa[c].mean() else ('decrease' if bb[c].mean()<aa[c].mean() else 'unchanged'))
    row[c+'_delta']=bb[c].mean()-aa[c].mean() if ok else np.nan
    row[c+'_share_previous']=aa[c].mean()/aa.p_kg.mean() if ok and aa.p_kg.mean()>0 else np.nan
    row[c+'_share_current']=bb[c].mean()/bb.p_kg.mean() if ok and bb.p_kg.mean()>0 else np.nan
    row[c+'_common_support_total_class']=classify(bb.p_kg.mean()/aa.p_kg.mean()-1,bb[yvar].mean()/aa[yvar].mean()-1) if ok and aa.p_kg.mean()>0 else 'unclassified'
   # Strict common support for all three categories and total/output.
   aa=prev[prev.major_complete];bb=gg[gg.major_complete];ok=len(aa)>=minobs and len(bb)>=minobs and aa.p_kg.mean()>0
   row['all_classes_comparable']=bool(ok)
   row['common_support_classification']=classify(bb.p_kg.mean()/aa.p_kg.mean()-1,bb[yvar].mean()/aa[yvar].mean()-1) if ok else 'unclassified'
   row['any_class_increase_common']=bool((bb[major].mean()>aa[major].mean()).any()) if ok else None
   res.append(row)
 return pd.DataFrame(res)
base=comparisons(panel);save(base,'decoupling_5yr_country.csv')
summary=base.groupby(['period_start','period','classification']).size().reset_index(name='country_comparisons');save(summary,'decoupling_5yr_summary.csv')
per=[]
for iso,g in base.groupby('iso3'):
 best=cur=absbest=abscur=0;last=None
 for _,r in g.sort_values('period_start').iterrows():
  adjacent=last is not None and r.period_start==last+5
  cur=(cur+1 if adjacent else 1) if r.decoupled else 0
  abscur=(abscur+1 if adjacent else 1) if r.classification=='absolute' else 0
  best=max(best,cur);absbest=max(absbest,abscur);last=r.period_start
 per.append(dict(iso3=iso,comparisons=len(g),max_consecutive_decoupled=best,max_consecutive_absolute=absbest,persistent_ge2=best>=2,persistent_ge3=best>=3))
persist=pd.DataFrame(per);save(persist,'persistence_country.csv')
trans=[]
for iso,g in base.groupby('iso3'):
 g=g.sort_values('period_start');prev=None
 for _,r in g.iterrows():
  if prev is not None and r.period_start==prev.period_start+5:trans.append(dict(iso3=iso,from_class=prev.classification,to_class=r.classification,to_period=r.period_start))
  prev=r
transition=pd.DataFrame(trans);save(transition,'decoupling_transitions_country.csv');save(transition.groupby(['from_class','to_class']).size().reset_index(name='n'),'transition_matrix.csv')
dec=base[base.decoupled];classrows=[]
for c in major:
 for direction,n in dec[c+'_direction'].value_counts().items():classrows.append(dict(category=c,direction=direction,n=int(n),total_decoupled=len(dec)))
save(pd.DataFrame(classrows),'class_changes_in_total_decoupling.csv')
strict=base[base.all_classes_comparable & base.common_support_classification.isin(['absolute','relative'])];save(strict,'common_support_decoupled_class_changes.csv')
# Category mass audit: residual is not toxicity and not necessarily an official other category.
audit=panel[panel.major_complete & panel.p_kg.notna()].copy();audit['major_exceeds_total']=audit.residual_kg < -np.maximum(.01,abs(audit.p_kg)*1e-5);save(audit[['iso3','year','p_kg']+major+['residual_kg','major_exceeds_total']],'category_additivity_audit.csv')
breaks=panel.sort_values(['iso3','year']).copy();prev=breaks.groupby('iso3').p_kg.shift();dy=breaks.groupby('iso3').year.diff();breaks['p_growth_annual']=breaks.p_kg/prev.where(prev>0)-1
save(breaks[dy.eq(1)&(breaks.p_growth_annual.abs()>1)][['iso3','year','p_kg','p_kg_flag','p_growth_annual']],'large_annual_breaks.csv')
# Balanced samples: no changing-country-composition global trend.
trendrows=[]
for start in [1990,2000,2010]:
 d=panel[panel.year.ge(start)];ct=d[d.joint_py].groupby('iso3').size();keep=ct[ct.eq(2025-start)].index
 for y,g in d[d.iso3.isin(keep)].groupby('year'):
  trendrows.append(dict(base_year=start,year=y,countries=len(keep),p_kg=g.p_kg.sum(),crop_usd2015=g.crop_usd2015.sum(),cropland_ha=g.cropland_ha.sum(min_count=len(g))))
trends=pd.DataFrame(trendrows)
for start in trends.base_year.unique():
 ix=trends.base_year.eq(start);ref=trends[ix & trends.year.eq(start)].iloc[0]
 for c in ['p_kg','crop_usd2015']:trends.loc[ix,c+'_index']=trends.loc[ix,c]/ref[c]*100
save(trends,'balanced_global_trends.csv')
regional=panel.groupby(['region','year']).agg(p_kg=('p_kg',lambda x:x.sum(min_count=1)),crop_usd2015=('crop_usd2015',lambda x:x.sum(min_count=1)),p_countries=('p_kg','count'),y_countries=('crop_usd2015','count')).reset_index();save(regional,'regional_unbalanced_totals.csv')
des=panel[['p_kg','crop_usd2015','cropland_ha','p_kg_ha','p_kg_per_1000usd','gdp_pc_const']+major].describe(percentiles=[.01,.25,.5,.75,.99]).T.reset_index(names='variable');save(des,'descriptive_statistics.csv')
# Sensitivities: each row reports its actual denominator.
sens=[]
def add_sens(name,data,**kw):
 b=comparisons(data,**kw)
 if len(b):sens.append(dict(specification=name,countries=b.iso3.nunique(),comparisons=len(b),absolute=int(b.classification.eq('absolute').sum()),relative=int(b.classification.eq('relative').sum()),decoupled_share=b.decoupled.mean()))
for start in [1990,2000,2010]:add_sens('start_'+str(start),panel,start=start)
add_sens('min4of5',panel,minobs=4)
add_sens('p_official_A_only',panel[panel.p_kg_flag.eq('A')])
add_sens('official_rich_A_ge50pct',panel[panel.official_share.ge(.5)])
add_sens('imputation_heavy_I_ge50pct',panel[panel.imputed_share.ge(.5)])
add_sens('QI_crop_output',panel,yvar='qi_crop');add_sens('WDI_crop_output',panel,yvar='wb_crop_index')
add_sens('3year_min2',panel,width=3,minobs=2)
add_sens('rolling5_min3',panel,rolling=True)
med=panel.groupby('iso3').cropland_ha.median();keep=med[med>=med.quantile(.05)].index;add_sens('exclude_smallest5pct_cropland',panel[panel.iso3.isin(keep)])
add_sens('all_major_classes_observed',panel[panel.major_complete])
save(pd.DataFrame(sens),'sensitivity_summary.csv')
# Country and year FE, cluster by country; associations only.
models=[];model_samples=[]
def fit(name,data,extra=False):
 cols=['p_kg','crop_usd2015','cropland_ha','gdp_pc_const'];d=data.copy()
 mixvars=[v for v in ['cereal_share','oil_share','vegfruit_share'] if d[v].notna().sum()>0] if extra else []
 d=d[d[cols].gt(0).all(axis=1)].dropna(subset=cols+mixvars)
 for c in cols:d['ln_'+c]=np.log(d[c])
 d=d[d.groupby('iso3').iso3.transform('size')>=5]
 if d.iso3.nunique()<20:return
 form='ln_p_kg ~ ln_crop_usd2015 + ln_cropland_ha + ln_gdp_pc_const'+(' + '+' + '.join(mixvars) if mixvars else '')+' + C(iso3) + C(year)'
 m=smf.ols(form,d).fit(cov_type='cluster',cov_kwds={'groups':d.iso3})
 for term in ['ln_crop_usd2015','ln_cropland_ha','ln_gdp_pc_const']+mixvars:
  ci=m.conf_int().loc[term];models.append(dict(model=name,term=term,coef=m.params[term],cluster_se=m.bse[term],pvalue=m.pvalues[term],ci_low=ci.iloc[0],ci_high=ci.iloc[1],nobs=int(m.nobs),countries=d.iso3.nunique(),first_year=d.year.min(),last_year=d.year.max(),country_FE=True,year_FE=True,cluster='iso3',r_squared=m.rsquared))
 model_samples.append(d[['iso3','year']].assign(model=name));print('fit',name,len(d),d.iso3.nunique(),flush=True)
fit('baseline_1990',panel);fit('crop_mix_1990',panel,True);fit('baseline_2000',panel[panel.year>=2000]);fit('baseline_2010',panel[panel.year>=2010]);fit('official_A_only',panel[panel.p_kg_flag.eq('A')]);fit('official_rich',panel[panel.official_share.ge(.5)])
save(pd.DataFrame(models),'fixed_effects_models.csv');save(pd.concat(model_samples),'model_sample_keys.csv')
# Figures.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':240})
colors=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#666666']
fig,axes=plt.subplots(1,2,figsize=(11,4))
for ax,start in zip(axes,[1990,2000]):
 g=trends[trends.base_year.eq(start)];ax.plot(g.year,g.p_kg_index,label='Pesticide active ingredient',color=colors[0]);ax.plot(g.year,g.crop_usd2015_index,label='Real crop gross output',color=colors[1]);ax.set(xlabel='Year',ylabel=f'Index ({start} = 100)',title=f'Balanced country sample, n = {int(g.countries.iloc[0])}');ax.legend(frameon=False,fontsize=8)
fig.tight_layout();fig.savefig(F/'Fig1_global_trends.png');plt.close(fig)
classes=['absolute','relative','growth_coupling','output_decline_p_decrease','output_decline_p_nondecrease','output_stagnant'];labels=['Absolute','Relative','Growth coupling','Output down / P down','Output down / P not down','Output stagnant']
fig,axes=plt.subplots(1,2,figsize=(12,4.5),gridspec_kw={'width_ratios':[1.3,1]})
pt=summary.pivot_table(index='period_start',columns='classification',values='country_comparisons',fill_value=0).reindex(columns=classes,fill_value=0)
bottom=np.zeros(len(pt))
for c,l,col in zip(classes,labels,colors):axes[0].bar(pt.index.astype(str),pt[c],bottom=bottom,label=l,color=col);bottom+=pt[c].values
axes[0].set(xlabel='Start year of current five-year period',ylabel='Country comparisons');axes[0].legend(frameon=False,fontsize=7,ncol=2,loc='lower left',bbox_to_anchor=(0,1.01))
mx=transition.groupby(['from_class','to_class']).size().unstack(fill_value=0).reindex(index=classes,columns=classes,fill_value=0);im=axes[1].imshow(mx,cmap='Blues');axes[1].set_xticks(range(6),['Abs','Rel','Growth','Y-/P-','Y-/P+','Y='],rotation=45);axes[1].set_yticks(range(6),['Abs','Rel','Growth','Y-/P-','Y-/P+','Y=']);axes[1].set(xlabel='Next classification',ylabel='Previous classification')
for i in range(6):
 for j in range(6):axes[1].text(j,i,str(mx.iloc[i,j]),ha='center',va='center',fontsize=8,color='white' if mx.iloc[i,j]>mx.max().max()/2 else 'black')
fig.colorbar(im,ax=axes[1],shrink=.7);fig.tight_layout();fig.savefig(F/'Fig2_decoupling_transitions.png');plt.close(fig)
# Balanced category sample, requiring all 35 annual values and nonnegative residual.
good=panel.major_complete & panel.p_kg.gt(0) & panel.residual_kg.ge(-.01)
ct=panel[good].groupby('iso3').size();classids=ct[ct.eq(35)].index
mix=panel[panel.iso3.isin(classids)].groupby('year')[['p_kg']+major].sum();mix['Residual']=mix.p_kg-mix[major].sum(axis=1);save(mix.reset_index().assign(countries=len(classids)),'balanced_class_mix.csv')
fig,axes=plt.subplots(1,2,figsize=(11,4))
axes[0].stackplot(mix.index,*[(mix[c]/mix.p_kg*100).values for c in major+['Residual']],labels=['Herbicides','Insecticides','Fungicides/bactericides','Residual'],colors=colors[:4]);axes[0].set(xlabel='Year',ylabel='Share of total mass (%)',title=f'Balanced category sample, n = {len(classids)}');axes[0].legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.22),ncol=2)
counts=[]
for c in major:counts.append([int((dec[c+'_direction']==s).sum()) for s in ['increase','decrease','unchanged','unclassified']])
vals=np.array(counts);bottom=np.zeros(3)
for i,l in enumerate(['Increase','Decrease','Unchanged','Not comparable']):axes[1].bar(range(3),vals[:,i],bottom=bottom,label=l,color=colors[i]);bottom+=vals[:,i]
axes[1].set_xticks(range(3),['Herbicides','Insecticides','Fungicides /\nbactericides']);axes[1].set(ylabel='Total-decoupled country comparisons');axes[1].legend(frameon=False,fontsize=7,loc='upper center',bbox_to_anchor=(.5,-.22),ncol=2)
fig.tight_layout();fig.savefig(F/'Fig3_pesticide_class_mix.png');plt.close(fig)
heat=panel.pivot(index='iso3',columns='year',values='joint_py').sort_index();fig,ax=plt.subplots(figsize=(10,9));ax.imshow(heat,aspect='auto',cmap='Blues',interpolation='none',vmin=0,vmax=1);ax.set_xticks(range(0,35,5),range(1990,2025,5));ax.set_yticks(range(0,len(heat),10),heat.index[::10]);ax.set(xlabel='Year',ylabel='Country / territory (ISO3)');fig.tight_layout();fig.savefig(F/'FigS1_country_year_coverage.png');plt.close(fig)
# Compact machine-readable run results for the Chinese report.
metrics=dict(panel_grid_countries=panel.iso3.nunique(),panel_grid_rows=len(panel),joint_countries=panel.loc[panel.joint_py,'iso3'].nunique(),joint_country_years=int(panel.joint_py.sum()),comparisons=len(base),comparison_countries=base.iso3.nunique(),absolute=int(base.classification.eq('absolute').sum()),relative=int(base.classification.eq('relative').sum()),persistent2=int(persist.persistent_ge2.sum()),persistent3=int(persist.persistent_ge3.sum()),strict_common_decoupled=len(strict),strict_common_any_increase=int(strict.any_class_increase_common.fillna(False).sum()),additivity_rows=len(audit),additivity_violations=int(audit.major_exceeds_total.sum()),class_balanced_countries=len(classids),large_breaks=int((breaks.p_growth_annual.abs()>1).sum()),flag_counts=panel.loc[panel.p_kg.notna(),'p_kg_flag'].value_counts().to_dict(),runtime=platform.python_version())
(M/'run_metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(metrics,ensure_ascii=False,indent=2),flush=True)
dictionary=[]
for c in panel.columns:
 unit='kg active ingredient' if c.endswith('_kg') or c=='p_kg' else ('ha' if c.endswith('_ha') else ('constant 2014-2016 US$' if c.endswith('usd2015') else 'see source metadata'))
 dictionary.append(dict(variable=c,unit=unit,definition={'joint_py':'P nonmissing and crop output >0','cropland_ha':'arable + permanent, both required','p_kg_per_1000usd':'P kg per 1000 constant crop gross output US$','region':'Current WDI region, descriptive classification only','income':'Current WDI income, not historical income class'}.get(c,c),missing_rule='No interpolation; source missing preserved; negatives invalid',provenance='FAO selected original rows + country_crosswalk.csv; WDI raw json; see source_register.csv'))
save(pd.DataFrame(dictionary),'data_dictionary.csv',M)
(M/'environment.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'platform':platform.platform()},indent=2),encoding='utf-8')
