"""Focused checks after primary estimates: quality persistence and event provenance."""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import sys,pathlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'decoupling'))
from build import HERE,P1,TABLE,MAJOR,save
from analyze import cluster_ratio
g=pd.read_csv(TABLE/'period_grid_including_undetermined.csv');b=pd.read_csv(TABLE/'country_period_all.csv');universe=set(b.iso3);g=g[g.iso3.isin(universe)].copy();p=pd.read_csv(HERE/'data/annual_reverified.csv');a=b[b.common_state.eq('absolute')].copy()
out=[];contin=[]
for label,d in [('all',g),('officialrich',g[g.national_A_share.ge(.5)]),('no_known_break',g.copy()),('officialrich_no_known_break',g[g.national_A_share.ge(.5)].copy()),('period_A_majority_no_known_break',g.copy())]:
 d=d.copy()
 if 'no_known_break' in label:d.loc[d.known_break,'state']='undetermined'
 if 'period_A_majority' in label:d.loc[~(d.p_A_share0.ge(.5)&d.p_A_share1.ge(.5)),'state']='undetermined'
 rows=[];trs=[]
 for iso,c in d.sort_values(['iso3','period_start']).groupby('iso3'):
  obs=c.state.ne('undetermined').to_numpy();ab=c.state.eq('absolute').to_numpy();cur=best=0
  for flag in ab:cur=cur+1 if flag else 0;best=max(best,cur)
  eligible2=bool((obs[:-1]&obs[1:]).any());eligible3=bool((obs[:-2]&obs[1:-1]&obs[2:]).any());rows.append(dict(iso3=iso,run=best,eligible2=eligible2,eligible3=eligible3))
  prev=None
  for _,r in c.iterrows():
   if prev is not None and prev.state=='absolute' and r.state!='undetermined':trs.append(dict(iso3=iso,continued=r.state=='absolute'))
   prev=r
 rows=pd.DataFrame(rows)
 for k in [2,3]:
  el=rows[rows['eligible'+str(k)]];out.append(dict(sample=label,run_length=k,eligible_countries=len(el),countries_with_run=int(el.run.ge(k).sum()),rate=el.run.ge(k).mean()))
 if trs:contin.append(dict(sample=label,**cluster_ratio(pd.DataFrame(trs),'continued')))
save(pd.DataFrame(out),'absolute_persistence_by_quality.csv');save(pd.DataFrame(contin),'absolute_continuation_by_quality.csv')
# All main triples have own-end support; enrich record-level flag composition for categories.
rows=[]
for _,r in a.iterrows():
 rr=dict(id=r.id,iso3=r.iso3,period_start=r.period_start)
 for side,years in [('0',r.common_support0),('1',r.common_support1)]:
  h=p[p.iso3.eq(r.iso3)&p.year.isin([int(x) for x in years.split(';')])]
  flags=h[[c+'_flag' for c in MAJOR]].to_numpy().flatten();rr['class_A_share'+side]=sum(x=='A' for x in flags)/len(flags);rr['class_AX_share'+side]=sum(x in ['A','X'] for x in flags)/len(flags)
 rows.append(rr)
flags=pd.DataFrame(rows);save(flags,'absolute_class_flag_composition.csv')
q=a.merge(flags,on=['id','iso3','period_start'],validate='one_to_one');quality=[]
for label,d in [('P_and_HIF_A_majority_both',q[q.p_A_share0.ge(.5)&q.p_A_share1.ge(.5)&q.class_A_share0.ge(.5)&q.class_A_share1.ge(.5)&~q.known_break]),('P_and_HIF_AX_majority_both',q[(q.p_A_share0+q.p_X_share0).ge(.5)&(q.p_A_share1+q.p_X_share1).ge(.5)&q.class_AX_share0.ge(.5)&q.class_AX_share1.ge(.5)&~q.known_break])]:
 for ev in ['any_increase','any_gt10','substantial_gt1pctP0']:quality.append(dict(sample=label,event=ev,**cluster_ratio(d,ev)))
save(pd.DataFrame(quality),'class_and_total_flag_quality.csv')
# Which five excluded countries/periods explain undefined relative category growth?
save(a[a.positive_class_baselines.eq(False)],'nonpositive_category_baselines.csv')
# Same sample outlier-drop indicators and magnitude deadband: ratio denominator is total P0.
for c in MAJOR:
 a[c+'_delta_over_P0']=pd.to_numeric(a[c+'_delta_over_P0'],errors='coerce')
rare=a[a.max_increment_over_P0.gt(.05)].sort_values('max_increment_over_P0',ascending=False);save(rare,'large_class_divergence_events.csv')
dead=[]
for margin in [0,.005,.01,.02]:
 h=a[(a.cgp < -margin)&(a.cgy > margin)]
 for ev in ['any_increase','any_gt10','substantial_gt1pctP0']:dead.append(dict(growth_margin=margin,event=ev,**cluster_ratio(h,ev)))
save(pd.DataFrame(dead),'total_output_growth_deadband_sensitivity.csv')
# Exploratory formal interactions: differing significance in two strata is not evidence of different slopes.
inter=[];keys=[]
d=b[~b.known_break].copy();d['officialrich']=d.national_A_share.ge(.5).astype(int);d['abs_y']=d.state.eq('absolute').astype(float);d['div_y']=pd.to_numeric(d.substantial_gt1pctP0,errors='coerce')
d['ln_income']=np.log(d.gdp_pc_const0.where(d.gdp_pc_const0>0));d['ln_area']=np.log(d.cropland_ha0.where(d.cropland_ha0>0));d['ln_intensity']=np.log((d.y0/d.cropland_ha0).where(d.y0.gt(0)&d.cropland_ha0.gt(0)))
for outcome in ['abs_y','div_y']:
 h=d.dropna(subset=[outcome,'ln_income','ln_area','ln_intensity']).copy();h=h[h.groupby('iso3').iso3.transform('size').ge(2)]
 form=outcome+' ~ ln_income + ln_area + ln_intensity + ln_income:officialrich + ln_area:officialrich + ln_intensity:officialrich + C(iso3) + C(period_start)'
 m=smf.ols(form,h).fit(cov_type='cluster',cov_kwds={'groups':h.iso3},use_t=True)
 for term in ['ln_income:officialrich','ln_area:officialrich','ln_intensity:officialrich']:
  ci=m.conf_int().loc[term];inter.append(dict(outcome=outcome,term=term,coefficient=m.params[term],cluster_se=m.bse[term],ci_low=ci.iloc[0],ci_high=ci.iloc[1],pvalue=m.pvalues[term],nobs=len(h),countries=h.iso3.nunique(),status='exploratory_after_primary_results'))
 keys.append(h[['id','iso3','period_start']].assign(outcome=outcome))
save(pd.DataFrame(inter),'quality_tier_interaction_tests.csv');save(pd.concat(keys),'quality_interaction_sample_keys.csv')
shares=[]
for c in MAJOR:
 h=a.copy();h['share_up']=h[c+'_share_pp'].gt(0);h['share_up_mass_not_up']=h[c+'_share_pp'].gt(0)&h[c+'_delta_t'].le(0)
 for ev in ['share_up','share_up_mass_not_up']:shares.append(dict(category=c,event=ev,**cluster_ratio(h,ev)))
shares=pd.DataFrame(shares);save(shares,'class_share_quantity_crossclassification.csv')
scale=[]
for label,h in [('all_common_absolute',a),('officialrich_common_absolute',a[a.national_A_share.ge(.5)])]:
 for weight in ['cp0','cy0','cropland_ha0']:
  for ev in ['any_increase','any_gt10','substantial_gt1pctP0','substantial_gt5pctP0']:scale.append(dict(sample=label,weight=weight,event=ev,**cluster_ratio(h,ev,weight)))
save(pd.DataFrame(scale),'agricultural_scale_weighted_incidence.csv')
print('Additional quality and provenance checks complete.')
