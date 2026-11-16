"""Offline staging plus original-row validation; never downloads a live database."""
from pathlib import Path
import shutil,sys,json
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[2];P1=R/'.reproduction/phase1';P2=R/'.reproduction/phase2'
for p in [P1/'metadata',P1/'results/tables',P1/'data/processed',P2/'data',P2/'tables',P2/'metadata',P2/'logs'] :p.mkdir(parents=True,exist_ok=True)
for n in ['country_notes_coded_changes.csv','country_proxy_notes.csv','country_crosswalk.csv']:shutil.copy2(R/'data/metadata'/n,P1/'metadata'/n)
shutil.copy2(R/'data/metadata/country_quality.csv',P1/'results/tables/country_quality.csv')
for n in ['official_world_and_matched_country_totals.csv','balanced_global_trends.csv']:shutil.copy2(R/'results/tables'/n,P1/'results/tables'/n)
shutil.copy2(R/'data/processed/panel_country_year.csv',P1/'data/processed/panel_country_year.csv')
panel=pd.read_csv(R/'data/processed/annual_reverified.csv');source=pd.read_csv(R/'data/processed/original_source_rows.csv');assert not source.duplicated(['iso3','Year','variable']).any()
for variable,d in source.groupby('variable'):
    expected_unit='1000 USD' if variable=='crop_usd2015' else 't';assert d.Unit.eq(expected_unit).all()
    rebuilt=d[['iso3','Year','Value']].rename(columns={'Year':'year','Value':'reconstructed'})
    rebuilt['reconstructed']=rebuilt.reconstructed*1000
    j=panel[['iso3','year',variable]].merge(rebuilt,on=['iso3','year'],how='left',validate='one_to_one')
    assert np.allclose(j[variable],j.reconstructed,equal_nan=True),variable
    # Rebuild the analysis values from the retained original observations.
    panel[variable]=j.reconstructed.to_numpy()
assert panel.joint_py.eq(panel.p_kg.notna()&panel.crop_usd2015.gt(0)).all()
panel.to_csv(P2/'data/annual_reverified.csv',index=False,encoding='utf-8-sig')
shutil.copy2(R/'data/processed/original_source_rows.csv',P2/'data/original_source_rows.csv')
sys.path.insert(0,str(R/'scripts/decoupling'));from build import comparisons
grid=comparisons(panel);grid.to_csv(P2/'tables/period_grid_including_undetermined.csv',index=False,encoding='utf-8-sig')
base=grid[grid.state.ne('undetermined')].copy();base.to_csv(P2/'tables/country_period_all.csv',index=False,encoding='utf-8-sig')
common=base[base.common_state.eq('absolute')];common.to_csv(P2/'tables/absolute_common_support_189.csv',index=False,encoding='utf-8-sig')
# Export once for the subsequent original inference script; frozen files stay untouched.
frozen=pd.read_csv(R/'results/tables/country_period_all.csv');j=base.merge(frozen,on='id',suffixes=('_new','_frozen'),validate='one_to_one');assert len(j)==893
for n in ['gp','gy','cp0','cp1','cy0','cy1','cgp','cgy','max_increment_over_P0']:assert np.allclose(j[n+'_new'],j[n+'_frozen'],equal_nan=True),n
for n in ['state','common_state']:assert j[n+'_new'].eq(j[n+'_frozen']).all(),n
result=dict(countries=int(panel[panel.joint_py].iso3.nunique()),country_years=int(panel.joint_py.sum()),comparisons=len(base),absolute=int(base.state.eq('absolute').sum()),relative=int(base.state.eq('relative').sum()),common_absolute=len(common),any_increase=int(common.any_increase.sum()))
(P2/'metadata/rebuilt_counts.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
