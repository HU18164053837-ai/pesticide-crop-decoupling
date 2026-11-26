"""One bounded, same-year land-intensity extension; frozen EDS outputs are read-only."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1];BASE=R/'00_BASELINE/repository' if (R/'00_BASELINE/repository').exists() else R.parents[1];D=BASE/'data/processed';T=BASE/'results/tables'
SEED=20261126;B=4000;MAJOR=['herbicide_kg','insecticide_kg','fungicide_kg'];VARIABLES=['p_kg','crop_usd2015','cropland_ha']
def classify(gp,gy):
    if not np.isfinite(gp) or not np.isfinite(gy):return 'undetermined'
    if gp<0 and gy>0:return 'absolute'
    if gy>0 and 0<=gp<gy:return 'relative'
    return 'coupled_or_other'
def ratio(d,event,seed):
    if not len(d):return dict(n=0,countries=0,events=0,estimate=None,ci_low=None,ci_high=None,ci_status='empty')
    g=d.groupby('iso3')[event].agg(['sum','size']);point=g['sum'].sum()/g['size'].sum();lo=hi=np.nan
    if len(g)>=20:
        rng=np.random.default_rng(seed);counts=rng.multinomial(len(g),np.ones(len(g))/len(g),size=B);v=(counts@g['sum'].to_numpy())/(counts@g['size'].to_numpy());lo,hi=np.quantile(v,[.025,.975])
    return dict(n=len(d),countries=len(g),events=int(g['sum'].sum()),estimate=point,ci_low=lo,ci_high=hi,ci_status='country-cluster percentile95' if len(g)>=20 else 'not estimated:fewer than20clusters')
def main():
    p=pd.read_csv(D/'annual_reverified.csv');b=pd.read_csv(T/'country_period_all.csv');a=pd.read_csv(T/'absolute_common_support_189.csv');records=[];audits=[]
    assert len(b)==893 and b.state.eq('absolute').sum()==200 and len(a)==189 and a.substantial_gt1pctP0.sum()==80
    by={iso:g for iso,g in p.groupby('iso3')}
    for _,frozen in b.iterrows():
        g=by[frozen.iso3]
        for support,cols in [('PYA',VARIABLES),('PYAHIF',VARIABLES+MAJOR)]:
            row=dict(id=frozen.id,iso3=frozen.iso3,previous_start=int(frozen.previous_start),period_start=int(frozen.period_start),frozen_state=frozen.state,support=support,known_break=bool(frozen.known_break),national_A_share=frozen.national_A_share,in_frozen189=frozen.id in set(a.id))
            valid=True
            for side,st in [('0',int(frozen.previous_start)),('1',int(frozen.period_start))]:
                raw=g[g.year.between(st,st+4)];joint=raw[raw.p_kg.notna()&np.isfinite(raw.p_kg)&raw.crop_usd2015.gt(0)&np.isfinite(raw.crop_usd2015)]
                land=joint[joint.cropland_ha.gt(0)&np.isfinite(joint.cropland_ha)]
                use=land[land[cols].notna().all(axis=1)&np.isfinite(land[cols]).all(axis=1)&land.p_kg.ge(0)]
                row['joint_n'+side]=len(joint);row['land_n'+side]=len(land);row['n'+side]=len(use);row['years'+side]=';'.join(use.year.astype(str));row['lost_area_years'+side]=';'.join(joint[~joint.year.isin(land.year)].year.astype(str))
                valid=valid and len(use)>=3
                for c in cols:row[c+side]=use[c].mean() if len(use) else np.nan
                row['P_A_share'+side]=use.p_kg_flag.eq('A').mean() if len(use) else np.nan
                row['area_A_share'+side]=use[['arable_ha_flag','permanent_ha_flag']].eq('A').to_numpy().mean() if len(use) else np.nan
                for _,v in use.iterrows():audits.append(dict(id=frozen.id,support=support,side=side,iso3=frozen.iso3,year=int(v.year),**{c:v[c] for c in cols}))
            row['support_available']=bool(valid and row['p_kg0']>0)
            if row['support_available']:
                row['gp']=row['p_kg1']/row['p_kg0']-1;row['gy']=row['crop_usd20151']/row['crop_usd20150']-1;row['ga']=row['cropland_ha1']/row['cropland_ha0']-1
                row['I0']=row['p_kg0']/row['cropland_ha0'];row['I1']=row['p_kg1']/row['cropland_ha1'];row['gI']=row['I1']/row['I0']-1;row['common_year_state']=classify(row['gp'],row['gy']);row['state_retained']=row['common_year_state']==row['frozen_state']
                row['intensity_up']=bool(row['I1']>row['I0']);row['joint_divergence']=bool(row['gp']<0 and row['gy']>0 and row['intensity_up']);row['mass_up']=bool(row['gp']>0)
                if support=='PYAHIF':
                    anyup=False
                    for c in MAJOR:
                        row[c+'_I0']=row[c+'0']/row['cropland_ha0'];row[c+'_I1']=row[c+'1']/row['cropland_ha1'];row[c+'_intensity_up']=bool(row[c+'_I1']>row[c+'_I0']);row[c+'_mass_up']=bool(row[c+'1']>row[c+'0']);row[c+'_mass_nonincrease_intensity_up']=bool(row[c+'1']<=row[c+'0'] and row[c+'_intensity_up']);anyup=anyup or row[c+'_intensity_up']
                    row['any_class_intensity_up']=anyup
            else:row['common_year_state']='undetermined';row['state_retained']=False
            records.append(row)
    d=pd.DataFrame(records);d.to_csv(R/'01_ANALYSIS/comparison_records.csv',index=False,encoding='utf-8-sig');pd.DataFrame(audits).to_csv(R/'01_ANALYSIS/same_year_observation_audit.csv',index=False,encoding='utf-8-sig')
    primary=d[d.support.eq('PYA')];avail=primary[primary.support_available].copy();summaries=[];counter=0
    def add(label,h,event,definition):
        nonlocal counter
        counter+=1
        # Identical supporting records and event vectors share one random stream.
        canonical=h[['id','iso3',event]].sort_values('id').rename(columns={event:'event'}).to_csv(index=False,lineterminator='\r\n')
        seed=(SEED+int(hashlib.sha256(canonical.encode()).hexdigest()[:8],16))%(2**32)
        summaries.append(dict(cohort=label,event=event,definition=definition,seed=seed,**ratio(h,event,seed)))
    for st in ['absolute','relative','coupled_or_other']:
        h=avail[avail.frozen_state.eq(st)];add('frozen_'+st,h,'intensity_up','Any intensity increase in frozen-state cohort,three-variable support')
        hstable=h[h.state_retained];add('stable_'+st,hstable,'intensity_up','Same-year state agrees with frozen state;descriptive comparable cohort')
        h2=avail[avail.common_year_state.eq(st)];add('common_year_'+st,h2,'intensity_up','Reclassified common-year state:sensitivity only')
    fa=avail[avail.frozen_state.eq('absolute')];stable=fa[fa.common_year_state.eq('absolute')]
    add('frozen_absolute_joint',fa,'joint_divergence','Same-year P declines,Y grows,I rises,denominator all land-supported frozen absolute records')
    add('stable_absolute_joint',stable,'joint_divergence','Same-year P declines,Y grows,I rises,among sign-retaining land-supported frozen absolute records')
    for lab,h in [('official_rich',stable[stable.national_A_share.ge(.5)]),('no_known_break',stable[~stable.known_break]),('official_rich_no_break',stable[stable.national_A_share.ge(.5)&~stable.known_break]),('P_and_area_A_majority_no_break',stable[stable.P_A_share0.ge(.5)&stable.P_A_share1.ge(.5)&stable.area_A_share0.ge(.5)&stable.area_A_share1.ge(.5)&~stable.known_break])]:add(lab,h,'joint_divergence','Reporting-quality restriction on unchanged same-year rows')
    cl=d[d.support.eq('PYAHIF')&d.support_available&d.in_frozen189].copy();cla=cl[cl.common_year_state.eq('absolute')]
    for lab,h in [('class_land_supported_frozen189',cl),('class_land_common_absolute',cla)]:
        add(lab,h,'any_class_intensity_up','Six-variable common support;any class kg/ha increases')
        for c in MAJOR:
            for ev in ['intensity_up','mass_up','mass_nonincrease_intensity_up']:add(lab,h,c+'_'+ev,'Six-variable common-year class comparison;zero baselines retained for level changes')
    out=pd.DataFrame(summaries);out.to_csv(R/'02_TABLES/intensity_cluster_estimates.csv',index=False,encoding='utf-8-sig')
    flow=[]
    for support,original in [('PYA',b),('PYAHIF',a)]:
        for state in (['absolute','relative','coupled_or_other'] if support=='PYA' else ['absolute']):
            ids=set(original[original.state.eq(state)].id) if support=='PYA' else set(a.id);h=d[d.support.eq(support)&d.id.isin(ids)];ok=h[h.support_available]
            flow.append(dict(support=support,frozen_state=state,original_n=len(h),land_supported_n=len(ok),excluded_n=len(h)-len(ok),retained_state_n=int(ok.state_retained.sum()),reclassified_n=int((~ok.state_retained).sum()),countries=int(ok.iso3.nunique()),missing_land_years=int(sum(h['joint_n0']-h['land_n0'])+sum(h['joint_n1']-h['land_n1']))))
    pd.DataFrame(flow).to_csv(R/'02_TABLES/support_flow.csv',index=False,encoding='utf-8-sig');pd.crosstab(avail.frozen_state,avail.common_year_state).to_csv(R/'02_TABLES/frozen_to_common_year_state.csv',encoding='utf-8-sig')
    events=stable[stable.joint_divergence];events.to_csv(R/'02_TABLES/joint_divergence_events.csv',index=False,encoding='utf-8-sig')
    classflow=dict(original_n=189,land_supported_n=len(cl),same_year_absolute_n=len(cla),state_changes=len(cl)-len(cla))
    result=dict(seed=SEED,bootstrap_replicates=B,minimum_clusters_for_interval=20,primary_support='same observed P,Y,positive cropland years;at least3 per5yr group;ratio of means',original_893=len(b),original_200=200,original_189=189,area_complete_893=len(avail),area_complete_frozen_absolute=len(fa),same_year_sign_retaining_absolute=len(stable),joint_events=len(events),event_countries=events.iso3.nunique(),class_support=classflow,new_models=0,threshold_search=False)
    (R/'01_ANALYSIS/SUMMARY.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));print(out.iloc[:15].to_string(index=False))
if __name__=='__main__':main()
