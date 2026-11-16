"""Clustered descriptive inference, absolute persistence, prespecified panel LPMs."""
import pathlib,json,time
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
import sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'decoupling'))
from build import HERE,P1,TABLE,MAJOR,save,comparisons
SEED=20261106;B=4000;rng=np.random.default_rng(SEED)
def boolval(s):return pd.to_numeric(s,errors='coerce')
def cluster_ratio(data,event='any_increase',weight=None):
 d=data.dropna(subset=[event]).copy()
 if len(d)==0:return dict(n=0,countries=0,events=0,estimate=np.nan,ci_low=np.nan,ci_high=np.nan)
 d['_w']=1. if weight is None else d[weight];d=d[d._w.gt(0)&d._w.notna()];d['_num']=boolval(d[event])*d._w
 g=d.groupby('iso3').agg(num=('_num','sum'),den=('_w','sum'));est=g.num.sum()/g.den.sum();ci=[np.nan,np.nan]
 if len(g)>=20:
  counts=rng.multinomial(len(g),np.ones(len(g))/len(g),size=B);den=counts@g.den.to_numpy();rep=(counts@g.num.to_numpy())/den;ci=np.quantile(rep,[.025,.975])
 return dict(n=len(d),countries=len(g),events=float(boolval(d[event]).sum()),estimate=float(est),ci_low=float(ci[0]),ci_high=float(ci[1]))
def paired_contrast(a,b,event='any_increase'):
 a=a.dropna(subset=[event]);b=b.dropna(subset=[event]);ids=sorted(set(a.iso3)|set(b.iso3));vals=[]
 for d in [a,b]:
  g=d.assign(_num=boolval(d[event])).groupby('iso3').agg(num=('_num','sum'),den=('_num','size')).reindex(ids,fill_value=0);vals.append(g)
 if len(a)==0 or len(b)==0:return {}
 est=vals[0].num.sum()/len(a)-vals[1].num.sum()/len(b);cnt=rng.multinomial(len(ids),np.ones(len(ids))/len(ids),size=B);ra=(cnt@vals[0].num)/(cnt@vals[0].den);rb=(cnt@vals[1].num)/(cnt@vals[1].den);rep=ra-rb;rep=rep[np.isfinite(rep)];lo,hi=np.quantile(rep,[.025,.975]);return dict(difference=est,ci_low=lo,ci_high=hi,n_target=len(a),n_comparator=len(b),cluster_countries=len(ids))
def incidence(d,label):
 out=[]
 for event in ['any_increase','any_gt5','any_gt10','substantial_gt1pctP0','substantial_gt5pctP0','any_zero_entry']:
  for weight in [None,'cp0']:
   out.append(dict(sample=label,event=event,estimand='country_period_unweighted' if weight is None else 'baseline_P_mass_weighted',**cluster_ratio(d,event,weight)))
 g=d.dropna(subset=['any_increase']).copy();agg=g.groupby('iso3').any_increase.agg(['mean','max']).reset_index()
 for v,l in [('mean','country_average_period_rate'),('max','country_ever_fraction')]:out.append(dict(sample=label,event='any_increase',estimand=l,**cluster_ratio(agg.rename(columns={v:'value'}),'value')))
 return out
def main():
 panel=pd.read_csv(HERE/'data/annual_reverified.csv');base=pd.read_csv(TABLE/'country_period_all.csv');grid=pd.read_csv(TABLE/'period_grid_including_undetermined.csv');common=base[base.common_complete].copy();absolute=common[common.common_state.eq('absolute')].copy()
 print('Phase2 analysis starts',flush=True)
 allrows=[]
 groups={'absolute':absolute,'all_output_growth':common[common.cgy.gt(0)],'all_total_decline':common[common.cgp.lt(0)],'decline_without_output_growth':common[common.cgp.lt(0)&common.cgy.le(0)]}
 for lab,d in groups.items():allrows.extend(incidence(d,lab))
 save(pd.DataFrame(allrows),'divergence_incidence_bootstrap.csv')
 contrasts=[]
 for lab in list(groups)[1:]:
  for event in ['any_increase','any_gt10','substantial_gt1pctP0','substantial_gt5pctP0']:contrasts.append(dict(target='absolute',comparator=lab,event=event,**paired_contrast(absolute,groups[lab],event)))
 save(pd.DataFrame(contrasts),'comparable_period_contrasts.csv')
 mags=[]
 for label,d in [('absolute_all',absolute),('absolute_any_increase',absolute[boolval(absolute.any_increase).eq(1)])]:
  for col in ['positive_increment_over_P0','max_increment_over_P0','class_coverage0','class_coverage1','delta_total_t','delta_residual_t']:
   x=d[col].dropna();row=dict(sample=label,variable=col,n=len(x),mean=x.mean(),**{'q'+str(int(q*100)):x.quantile(q) for q in [.05,.25,.5,.75,.95]});mags.append(row)
 save(pd.DataFrame(mags),'divergence_magnitude_distribution.csv')
 classrows=[];long=[]
 for c in MAJOR:
  d=absolute.copy();d['class_increase']=d[c+'_delta_t'].gt(0);d['class_gt5']=np.where(d[c+'0'].gt(0),d[c+'_growth'].gt(.05),np.nan);d['class_gt10']=np.where(d[c+'0'].gt(0),d[c+'_growth'].gt(.1),np.nan);d['class_gt1pctP0']=d[c+'_delta_over_P0'].gt(.01)
  for event in ['class_increase','class_gt5','class_gt10','class_gt1pctP0']:classrows.append(dict(category=c,event=event,**cluster_ratio(d,event)))
  for _,r in d.iterrows():long.append(dict(id=r.id,iso3=r.iso3,period_start=r.period_start,category=c,start_t=r[c+'0']/1000,end_t=r[c+'1']/1000,delta_t=r[c+'_delta_t'],growth_pct=r[c+'_growth']*100,share_change_pp=r[c+'_share_pp'],increment_pct_P0=r[c+'_delta_over_P0']*100,known_break=r.known_break,national_A_share=r.national_A_share))
 save(pd.DataFrame(classrows),'category_incidence_bootstrap.csv');save(pd.DataFrame(long),'common_absolute_category_magnitudes.csv')
 strata=[]
 for key in ['region','income','period_start']:
  for level,d in absolute.groupby(key):
   for ev in ['any_increase','any_gt10','substantial_gt1pctP0']:strata.append(dict(stratum=key,level=level,event=ev,small_sample_descriptive_only=d.iso3.nunique()<20,**cluster_ratio(d,ev)))
 save(pd.DataFrame(strata),'divergence_by_region_income_period.csv')
 # Quality distributions of all yearly source flags, region and income retained.
 for key in ['region','income']:
  d=panel[panel.p_kg.notna()].groupby([key,'p_kg_flag']).size().reset_index(name='country_years');d['within_group_share']=d.country_years/d.groupby(key).country_years.transform('sum');save(d,'source_flags_by_'+key+'.csv')
 proxy=base.groupby(['region','income','proxy_note'],dropna=False).agg(comparisons=('id','size'),countries=('iso3','nunique')).reset_index();save(proxy,'proxy_screen_by_region_income.csv')
 # Exact arithmetic decomposition of six reported, mutually exclusive major headings.
 arithmetic=[]
 for _,r in absolute.iterrows():
  d0=panel[panel.iso3.eq(r.iso3)&panel.year.isin([int(x) for x in r.common_support0.split(';')])];d1=panel[panel.iso3.eq(r.iso3)&panel.year.isin([int(x) for x in r.common_support1.split(';')])]
  rr=dict(id=r.id,iso3=r.iso3,delta_total_t=r.delta_total_t,delta_HIF_t=r.delta_HIF_t,delta_residual_t=r.delta_residual_t,coverage0=r.class_coverage0,coverage1=r.class_coverage1)
  extra=['other_kg','growthreg_kg','rodent_kg'];all6=all(d[v].notna().all() for d in [d0,d1] for v in extra)
  rr['six_headings_complete']=all6;rr['mineral_oils_observed']=all(d.mineral_kg.notna().all() for d in [d0,d1]);delta_extra=sum(d1[v].mean()-d0[v].mean() for v in extra)/1000 if all6 else np.nan;rr['delta_3_other_reported_headings_t']=delta_extra;rr['unallocated_balance_delta_t']=r.delta_total_t-r.delta_HIF_t-delta_extra if all6 else np.nan
  arithmetic.append(rr)
 save(pd.DataFrame(arithmetic),'reported_heading_arithmetic.csv')
 # Nonoverlapping levels, overlapping adjacent changes: explicit full-grid state transition.
 trans=[]
 for iso,g in grid.sort_values(['iso3','period_start']).groupby('iso3'):
  prev=None
  for _,r in g.iterrows():
   if prev is not None and r.period_start==prev.period_start+5:trans.append(dict(iso3=iso,from_period=prev.period_start,to_period=r.period_start,from_state=prev.state,to_state=r.state,shared_middle_mean=True))
   prev=r
 transition=pd.DataFrame(trans);save(transition,'state_transition_opportunities.csv');mat=transition.groupby(['from_state','to_state']).size().reset_index(name='n');mat['row_probability']=mat.n/mat.groupby('from_state').n.transform('sum');save(mat,'transition_matrix_including_undetermined.csv')
 # Restrict matrix to country universe with >=1 determined comparison; no empty-country inflation.
 universe=set(base.iso3);t=transition[transition.iso3.isin(universe)];mat=t.groupby(['from_state','to_state']).size().reset_index(name='n');mat['row_probability']=mat.n/mat.groupby('from_state').n.transform('sum');save(mat,'transition_matrix_study_countries.csv')
 cont=[]
 for state in ['absolute','relative','coupled_or_other']:
  d=t[t.from_state.eq(state)&t.to_state.ne('undetermined')].copy();d['continued']=d.to_state.eq(state);cont.append(dict(from_state=state,**cluster_ratio(d,'continued')))
 save(pd.DataFrame(cont),'state_continuation_bootstrap.csv')
 pers=[]
 for iso,g in grid[grid.iso3.isin(universe)].sort_values(['iso3','period_start']).groupby('iso3'):
  best=cur=0
  for _,r in g.iterrows():cur=cur+1 if r.state=='absolute' else 0;best=max(best,cur)
  observed=g.state.ne('undetermined').to_numpy();op2=int((observed[1:]&observed[:-1]).sum());op3=int((observed[2:]&observed[1:-1]&observed[:-2]).sum())
  # Separated changes compare (1990–94->1995–99) then (2000–04->2005–09), no shared level mean.
  for offset in [0,1]:
   h=g.iloc[offset::2];a=h.state.eq('absolute').to_numpy();obs=h.state.ne('undetermined').to_numpy();opp=int((obs[1:]&obs[:-1]).sum());pairs=int((a[1:]&a[:-1]).sum());pers.append(dict(iso3=iso,design='separated_pairs_offset'+str(offset),max_consecutive_absolute=np.nan,observed_comparisons=int(obs.sum()),eligible_two_positions=opp,absolute_two_positions=pairs,has_two=pairs>0,eligible_three_positions=np.nan,has_three=np.nan))
  pers.append(dict(iso3=iso,design='adjacent_changes',max_consecutive_absolute=best,observed_comparisons=int(observed.sum()),eligible_two_positions=op2,absolute_two_positions=int((g.state.eq('absolute').to_numpy()[1:]&g.state.eq('absolute').to_numpy()[:-1]).sum()),has_two=best>=2,eligible_three_positions=op3,has_three=best>=3))
 per=pd.DataFrame(pers);save(per,'absolute_persistence_country.csv');pr=[]
 for label,g in per.groupby('design'):
  for k in [2,3]:
   field='eligible_two_positions' if k==2 else 'eligible_three_positions';outcome='has_two' if k==2 else 'has_three';d=g[g[field].gt(0)]
   if len(d):pr.append(dict(design=label,min_run=k,eligible_countries=len(d),countries_with_run=int(boolval(d[outcome]).sum()),rate=boolval(d[outcome]).mean(),eligible_windows=int(d[field].sum())))
 save(pd.DataFrame(pr),'absolute_persistence_summary.csv')
 # Sensitivity sampling precedes inference; all window definitions explicit and nonrolling.
 specs=[('primary',base),('official_rich_A_ge50pct',base[base.national_A_share.ge(.5)]),('period_A_majority_both',base[base.p_A_share0.ge(.5)&base.p_A_share1.ge(.5)]),('no_known_break',base[~base.known_break]),('official_rich_no_known_break',base[base.national_A_share.ge(.5)&~base.known_break]),('exclude_negative_residual_gt2t',base[~base.common_known_negative_residual.fillna(False).astype(bool)]),('through2019',base[base.period_end<=2019]),('exclude_BRA_CHN_USA',base[~base.iso3.isin(['BRA','CHN','USA'])])]
 complete=panel[panel.joint_py].groupby('iso3').size();balanced=set(complete[complete.eq(35)].index);specs.append(('balanced1990_2024',base[base.iso3.isin(balanced)]))
 for name,kw in [('min4of5',dict(minobs=4)),('3year_min2',dict(width=3,minobs=2)),('shift1992_5year',dict(start=1992)),('QI_output',dict(yvar='qi_crop')),('WDI_crop_index',dict(yvar='wb_crop_index')),('endpoints1994_1999etc',dict(endpoint=True)),('A_X_records_only',dict(allowed_flags=['A','X']))]:
  print('sensitivity',name,flush=True);d=comparisons(panel,**kw);save(d,'sensitivity_country_period_'+name+'.csv');specs.append((name,d[d.state.ne('undetermined')]))
 sensitivity=[]
 for label,d in specs:
  a=d[d.common_state.eq('absolute')];row=dict(specification=label,all_comparisons=len(d),all_countries=d.iso3.nunique(),absolute=int(d.state.eq('absolute').sum()),common_absolute=len(a),common_countries=a.iso3.nunique())
  for ev in ['any_increase','any_gt10','substantial_gt1pctP0','substantial_gt5pctP0']:
   result=cluster_ratio(a,ev);sensitivity.append(dict(**row,event=ev,weighting='unweighted',**result));sensitivity.append(dict(**row,event=ev,weighting='baseline_P_mass',**cluster_ratio(a,ev,'cp0')))
 save(pd.DataFrame(sensitivity),'core_sensitivity_bootstrap.csv')
 # Influence uses countries as units and exposes size-weighting concentration.
 influence=[]
 for iso in sorted(absolute.iso3.unique()):
  d=absolute[~absolute.iso3.eq(iso)];influence.append(dict(dropped_iso3=iso,n=len(d),unweighted=boolval(d.any_increase).mean(),mass_weighted=np.average(boolval(d.any_increase),weights=d.cp0),removed_mass_share=absolute[absolute.iso3.eq(iso)].cp0.sum()/absolute.cp0.sum()))
 save(pd.DataFrame(influence),'leave_one_country_out.csv')
 run_models(base)
 config=dict(seed=SEED,bootstrap_replicates=B,ci='country-cluster percentile95; minimum20countries',note='No ecological risk or causal mechanism inference');(HERE/'metadata/inference_config.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
 print('Analysis complete',flush=True)

def run_models(base):
 d=base.copy();d['absolute_outcome']=d.state.eq('absolute').astype(int);d['divergence_outcome']=boolval(d.substantial_gt1pctP0)
 d['ln_gdppc0']=np.log(d.gdp_pc_const0.where(d.gdp_pc_const0>0));d['ln_cropland0']=np.log(d.cropland_ha0.where(d.cropland_ha0>0));d['ln_crop_intensity0']=np.log((d.y0/d.cropland_ha0).where((d.y0>0)&(d.cropland_ha0>0)))
 d['cereal_share0_10pp']=d.cereal_share0/.1;d['vegfruit_share0_10pp']=d.vegfruit_share0/.1
 basic=['ln_gdppc0','ln_cropland0','ln_crop_intensity0'];extra=basic+['cereal_share0_10pp','vegfruit_share0_10pp'];out=[];samples=[];diagn=[];vifs=[];corr=[];gaps=[];skips=[]
 for v in basic+extra[3:]:gaps.append(dict(variable=v,comparisons=len(d),nonmissing=int(d[v].notna().sum()),countries=d[d[v].notna()].iso3.nunique()))
 save(pd.DataFrame(gaps),'model_covariate_availability.csv')
 def fit(label,frame,y,xvars):
  before=len(frame);g=frame.dropna(subset=[y]+xvars).copy();g=g[g.groupby('iso3').iso3.transform('size').ge(2)];G=g.iso3.nunique();switch=g.groupby('iso3')[y].nunique();events=g[y].sum()
  if G<20 or min(events,len(g)-events)<20:
   skips.append(dict(model=label,reason='fewer20clusters or20events/nonevents',n=len(g),clusters=G,events=events));return
  formula=y+' ~ '+' + '.join(xvars)+' + C(iso3) + C(period_start)';model=smf.ols(formula,g).fit(cov_type='cluster',cov_kwds={'groups':g.iso3},use_t=True)
  for term in xvars:
   ci=model.conf_int().loc[term];out.append(dict(model=label,outcome=y,term=term,coefficient=model.params[term],cluster_se=model.bse[term],ci_low=ci.iloc[0],ci_high=ci.iloc[1],pvalue=model.pvalues[term],nobs=len(g),clusters=G,countries_switching_outcome=int(switch.ge(2).sum()),events=events,fe='country+5year_period',reference='baseline prior5yr means',estimator='LPM',se='country cluster finite correction t(G-1)'))
  # Partial out country and period effects for within-design collinearity checks.
  fe=smf.ols(y+' ~ C(iso3) + C(period_start)',g).exog;X=g[xvars].to_numpy();within=X-fe@np.linalg.lstsq(fe,X,rcond=None)[0];sd=within.std(axis=0);Z=within/np.where(sd>1e-12,sd,1);rank=np.linalg.matrix_rank(Z);condition=np.linalg.cond(Z)
  for i,v in enumerate(xvars):vifs.append(dict(model=label,variable=v,within_sd=sd[i],within_vif=variance_inflation_factor(Z,i) if len(xvars)>1 else 1))
  cm=pd.DataFrame(Z,columns=xvars).corr()
  for v1 in xvars:
   for v2 in xvars:corr.append(dict(model=label,variable1=v1,variable2=v2,within_correlation=cm.loc[v1,v2]))
  pred=model.predict(g);diagn.append(dict(model=label,candidate_n=before,nobs=len(g),clusters=G,events=events,switching_countries=int(switch.ge(2).sum()),within_rank=rank,regressors=len(xvars),standardized_within_condition=condition,prediction_min=pred.min(),prediction_max=pred.max(),prediction_outside01=int(((pred<0)|(pred>1)).sum()),r2=model.rsquared,model_rank=np.linalg.matrix_rank(model.model.exog)))
  samples.append(g[['id','iso3','period_start',y]+xvars].assign(model=label));print('model',label,len(g),G,'switchers',int(switch.ge(2).sum()),flush=True)
 for label,frame,y in [('absolute_all',d,'absolute_outcome'),('absolute_officialrich',d[d.national_A_share.ge(.5)&~d.known_break],'absolute_outcome'),('class_divergence_all_common',d[d.common_complete],'divergence_outcome'),('class_divergence_officialrich',d[d.common_complete&d.national_A_share.ge(.5)&~d.known_break],'divergence_outcome'),('class_divergence_conditional_absolute',d[d.common_state.eq('absolute')],'divergence_outcome')]:
  fit(label+'_baseline',frame,y,basic);enhanced=frame.dropna(subset=extra);fit(label+'_baseline_same_enhanced_sample',enhanced,y,basic);fit(label+'_enhanced',frame,y,extra)
 save(pd.DataFrame(out),'economic_models_coefficients.csv');save(pd.concat(samples,ignore_index=True),'economic_model_sample_keys.csv');save(pd.DataFrame(diagn),'economic_model_diagnostics.csv');save(pd.DataFrame(vifs),'model_within_collinearity.csv');save(pd.DataFrame(corr),'model_within_correlations.csv');save(pd.DataFrame(skips,columns=['model','reason','n','clusters','events']),'models_not_estimable.csv')
if __name__=='__main__':main()
