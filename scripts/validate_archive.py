"""Numerical replay comparisons, source hashes and frozen-figure integrity."""
from pathlib import Path
import json,hashlib,platform,sys,importlib.metadata
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1];T=R/'results/tables';FRESH=R/'.reproduction/phase2/tables';rows=[]
def ck(name,observed,expected,kind='recomputed',rtol=1e-9,atol=1e-10):
    ok=bool(np.isclose(observed,expected,rtol=rtol,atol=atol)) if isinstance(observed,(int,float,np.number)) else observed==expected
    rows.append(dict(check=name,observed=observed,expected=expected,kind=kind,status='PASS' if ok else 'FAIL'))
    assert ok,(name,observed,expected)
p=pd.read_csv(R/'.reproduction/phase2/data/annual_reverified.csv');b=pd.read_csv(FRESH/'country_period_all.csv');a=b[b.common_state.eq('absolute')].copy()
for n,v,e in [('reporting_entities',p[p.joint_py].iso3.nunique(),157),('joint_country_years',p.joint_py.sum(),5251),('comparisons',len(b),893),('absolute',b.state.eq('absolute').sum(),200),('relative',b.state.eq('relative').sum(),133),('common_absolute',len(a),189),('any_class_increase',a.any_increase.sum(),102),('class_gt1pctP0',a.substantial_gt1pctP0.sum(),80)]:ck(n,int(v),e)
h=a[a.cgp.lt(-.01)&a.cgy.gt(.01)];ck('bilateral1pct_n',len(h),172);ck('bilateral1pct_k',int(h.substantial_gt1pctP0.sum()),72)
c=b[b.common_complete&b.cgp.lt(0)&b.cgy.le(0)];ck('comparator_n',len(c),122);ck('comparator_k',int(c.substantial_gt1pctP0.sum()),35)
for label,h,w,expected in [('crop_value',a,'cy0',.15630935391918524),('crop_value_exclude_CHN',a[~a.iso3.eq('CHN')],'cy0',.31158242987177565),('pesticide_mass',a,'cp0',.31307237832474316)]:ck(label,float(np.average(h.substantial_gt1pctP0.astype(float),weights=h[w])),expected)
# Independently recalculate full study-country transition matrix, including unknown states.
grid=pd.read_csv(FRESH/'period_grid_including_undetermined.csv');pairs=[]
for iso,g in grid[grid.iso3.isin(set(b.iso3))].sort_values(['iso3','period_start']).groupby('iso3'):
    for i in range(1,len(g)):
        x,y=g.iloc[i-1],g.iloc[i]
        if y.period_start==x.period_start+5:pairs.append(dict(from_state=x.state,to_state=y.state))
trans=pd.DataFrame(pairs);m=trans.groupby(['from_state','to_state']).size().reset_index(name='n');m['row_probability']=m.n/m.groupby('from_state').n.transform('sum');f=pd.read_csv(T/'transition_matrix_study_countries.csv');j=m.merge(f,on=['from_state','to_state'],suffixes=('_new','_frozen'),validate='one_to_one');ck('transition_cells',len(j),len(f))
for _,x in j.iterrows():ck('matrix_count '+x.from_state+' '+x.to_state,int(x.n_new),int(x.n_frozen));ck('matrix_probability '+x.from_state+' '+x.to_state,float(x.row_probability_new),float(x.row_probability_frozen))
d=trans[trans.from_state.eq('absolute')&trans.to_state.ne('undetermined')];ck('continuation_denominator',len(d),171);ck('continuation_events',int(d.to_state.eq('absolute').sum()),55)
# Validate exact historical replay tables if --with-inference ran.
inference=FRESH/'divergence_incidence_bootstrap.csv';replayed=[]
for new in sorted(FRESH.glob('*.csv')):
    if new.name in ['country_period_all.csv','period_grid_including_undetermined.csv','absolute_common_support_189.csv'] or not (T/new.name).exists():continue
    x=pd.read_csv(new);y=pd.read_csv(T/new.name);ck(new.name+' shape',str(x.shape),str(y.shape),'replayed_table');ck(new.name+' columns',list(x.columns),list(y.columns),'replayed_table')
    for col in x:
        if pd.api.types.is_numeric_dtype(x[col]) and pd.api.types.is_numeric_dtype(y[col]):
            ok=np.allclose(x[col],y[col],equal_nan=True,rtol=1e-7,atol=1e-9)
            ck(new.name+' '+col+' numeric_match',bool(ok),True,'replayed_table')
        else:ck(new.name+' '+col+' exact_match',x[col].fillna('<NA>').astype(str).tolist()==y[col].fillna('<NA>').astype(str).tolist(),True,'replayed_table')
    replayed.append(new.name)
# Never merge CI objects from distinct frozen runs.
q=pd.read_csv(T/'TableS15_absolute_continuation_quality.csv');cont=pd.read_csv(T/'state_continuation_bootstrap.csv');u=q[q['sample'].eq('all')].iloc[0];v=cont[cont.from_state.eq('absolute')].iloc[0]
ci=pd.DataFrame([dict(source='Table S15 / absolute_continuation_by_quality.csv',estimate=u.estimate,ci_low=u.ci_low,ci_high=u.ci_high),dict(source='main-text archived state_continuation_bootstrap.csv',estimate=v.estimate,ci_low=v.ci_low,ci_high=v.ci_high)])
ci.to_csv(R/'reproducibility/ci_endpoint_provenance.csv',index=False,encoding='utf-8-sig');ck('S15 exact copy',q.equals(pd.read_csv(T/'absolute_continuation_by_quality.csv')),True,'frozen_provenance')
for file in T.glob('*.csv'):
    tab=pd.read_csv(file)
    if {'ci_low','ci_high'}.issubset(tab.columns):
        valid=tab[['ci_low','ci_high']].notna().all(axis=1);ck(file.name+' CI_order',bool(tab.loc[valid,'ci_low'].le(tab.loc[valid,'ci_high']).all()),True,'frozen_interval_integrity')
source=pd.read_csv(R/'reproducibility/source_manifest.csv')
for _,x in source.iterrows():ck('source_hash '+x.repository_path,hashlib.sha256((R/x.repository_path).read_bytes()).hexdigest(),x.archived_sha256,'frozen_provenance')
figs=list((R/'results/figures').glob('*.svg'))+list((R/'results/figures').glob('*.png'));ck('v1.4_figure_files',len(figs),18,'frozen_figure_integrity')
pd.DataFrame(rows).to_csv(R/'reproducibility/validation_results.csv',index=False,encoding='utf-8-sig')
result=dict(status='REPRODUCIBILITY_VERIFIED' if inference.exists() else 'CORE_DETERMINISTIC_REPRODUCIBILITY_VERIFIED',checks=len(rows),failed=0,python=platform.python_version(),isolated_venv=sys.prefix!=sys.base_prefix,packages={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','matplotlib','pycountry','Pillow']},replayed_tables=replayed,figures='Frozen v1.4 originals checked by hash;not redrawn',scope='Offline RP/QV measured-row reconstruction,comparison grid,original bootstrap/model replay if present,boundary/CHN audit,frozen table/figure integrity',raw_bulk_database_full_rerun='Not run:no database retrieval;retained original-row extracts used instead')
(R/'reproducibility/validation_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

