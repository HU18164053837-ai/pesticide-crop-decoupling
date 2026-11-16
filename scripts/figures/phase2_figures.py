"""Publication-candidate English figures, all from real Phase2 tables."""
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'decoupling'))
from build import HERE,P1,TABLE,save
F=HERE/'figures';plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'svg.fonttype':'none','pdf.fonttype':42})
C=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#666666']
def export(fig,name):
 for ext in ['png','svg','pdf']:fig.savefig(F/(name+'.'+ext),bbox_inches='tight',pad_inches=.15)
 plt.close(fig)
def lab(ax,t):ax.text(0,1.10,t,transform=ax.transAxes,fontweight='bold',fontsize=12)
def forest(ax,data,labelcol='label',scale=100,color=C[0]):
 y=np.arange(len(data));x=data.estimate.to_numpy()*scale;lo=data.ci_low.to_numpy()*scale;hi=data.ci_high.to_numpy()*scale
 valid=np.isfinite(lo)&np.isfinite(hi);ax.errorbar(x[valid],y[valid],xerr=np.array([x[valid]-lo[valid],hi[valid]-x[valid]]),fmt='o',color=color,capsize=3,markersize=5);ax.plot(x[~valid],y[~valid],'o',color=color);ax.set_yticks(y,data[labelcol]);ax.invert_yaxis()
def main():
 p=pd.read_csv(HERE/'data/annual_reverified.csv');b=pd.read_csv(TABLE/'country_period_all.csv');a=pd.read_csv(TABLE/'absolute_common_support_189.csv');inc=pd.read_csv(TABLE/'divergence_incidence_bootstrap.csv');sen=pd.read_csv(TABLE/'core_sensitivity_bootstrap.csv')
 # Fig1: source World totals distinct from fixed study sample; coverage in its own panels.
 w=pd.read_csv(P1/'results/tables/official_world_and_matched_country_totals.csv');t=pd.read_csv(P1/'results/tables/balanced_global_trends.csv');fig,axs=plt.subplots(2,2,figsize=(12,8.3))
 ax=axs[0,0]
 for v,label,col in [('FAO_world_p_kg','Pesticide active ingredient',C[0]),('FAO_world_crop_usd2015','Real crop gross output',C[1])]:ax.plot(w.year,w[v]/w.loc[w.year.eq(1990),v].iloc[0]*100,label=label,color=col,lw=2)
 ax.set(xlabel='Year',ylabel='Index (1990 = 100)',title='FAOSTAT reported world series');ax.legend(frameon=False,fontsize=8);lab(ax,'a')
 ax=axs[0,1];d=t[t.base_year.eq(1990)]
 for v,label,col in [('p_kg_index','Pesticide active ingredient',C[0]),('crop_usd2015_index','Real crop gross output',C[1])]:ax.plot(d.year,d[v],label=label,color=col,lw=2)
 ax.set(xlabel='Year',ylabel='Index (1990 = 100)',title='Balanced study subset (116 countries/territories)');ax.legend(frameon=False,fontsize=8);lab(ax,'b')
 ax=axs[1,0];h=p.pivot(index='iso3',columns='year',values='joint_py').sort_index();ax.imshow(h,aspect='auto',cmap='Blues',interpolation='none',vmin=0,vmax=1);ax.set_xticks(range(0,35,5),range(1990,2025,5));ax.set_yticks(range(0,len(h),30),h.index[::30]);ax.set(xlabel='Year',ylabel='ISO3 (ordered)',title='Joint pesticide/crop-output availability');lab(ax,'c')
 ax=axs[1,1];flags=p[p.p_kg.notna()].groupby(['year','p_kg_flag']).size().unstack(fill_value=0).reindex(columns=['A','E','I','X'],fill_value=0);save(flags.reset_index(),'source_flags_by_year.csv');ax.stackplot(flags.index,*[flags[v].to_numpy() for v in flags],labels=['Official (A)','Estimated (E)','Imputed (I)','External (X)'],colors=[C[2],C[0],C[1],C[3]]);ax.set(xlabel='Year',ylabel='Country/territory records',ylim=(0,250),title='Pesticide total: original source flags');ax.legend(frameon=False,fontsize=8,ncol=2,loc='upper left');lab(ax,'d');fig.tight_layout();export(fig,'Fig1_global_changes_and_coverage')
 # Fig2: retain undetermined state and absolute persistence opportunities.
 fig,axs=plt.subplots(1,3,figsize=(16,5.0),gridspec_kw={'width_ratios':[1.1,1.15,1]});m=pd.read_csv(TABLE/'transition_matrix_study_countries.csv');states=['absolute','relative','coupled_or_other','undetermined'];short=['Absolute','Relative','Other','Unknown'];pt=m.pivot(index='from_state',columns='to_state',values='row_probability').reindex(index=states,columns=states).fillna(0);ct=m.pivot(index='from_state',columns='to_state',values='n').reindex(index=states,columns=states).fillna(0)
 ax=axs[0];im=ax.imshow(pt,cmap='Blues',vmin=0,vmax=1);ax.set_xticks(range(4),short,rotation=35,ha='right');ax.set_yticks(range(4),short);ax.set(xlabel='Next comparison',ylabel='Current comparison',title='Adjacent-change transitions')
 for i in range(4):
  for j in range(4):ax.text(j,i,f'{pt.iloc[i,j]*100:.0f}%\n(n={ct.iloc[i,j]:.0f})',ha='center',va='center',fontsize=8,color='white' if pt.iloc[i,j]>.6 else 'black')
 lab(ax,'a');ax=axs[1];d=pd.read_csv(TABLE/'absolute_continuation_by_quality.csv');names={'all':'All records','officialrich':'Official-rich countries','no_known_break':'Selected breaks removed','officialrich_no_known_break':'Official-rich, breaks removed','period_A_majority_no_known_break':'A-majority periods, breaks removed'};d['label']=d['sample'].map(names);forest(ax,d);ax.set(xlim=(0,100),xlabel='Remaining absolute in next comparison (%)',title='Absolute-state continuation');lab(ax,'b')
 ax=axs[2];d=pd.read_csv(TABLE/'absolute_persistence_summary.csv');d=d[d.min_run.eq(2)].copy();names={'adjacent_changes':'Adjacent changes\n(shared middle mean)','separated_pairs_offset0':'Separated changes\n(offset 0)','separated_pairs_offset1':'Separated changes\n(offset 1)'};ax.barh(range(len(d)),d.rate*100,color=[C[0],C[2],C[3]])
 for i,(_,r) in enumerate(d.iterrows()):ax.text(r.rate*100+1,i,f'{r.countries_with_run}/{r.eligible_countries}',va='center',fontsize=9)
 ax.set_yticks(range(len(d)),d.design.map(names));ax.invert_yaxis();ax.set(xlim=(0,40),xlabel='Countries with two absolute states (%)',title='Opportunity-dependent persistence');lab(ax,'c');fig.tight_layout();export(fig,'Fig2_absolute_decoupling_persistence')
 # Fig3: original incidence, actual magnitude, and period comparator.
 fig,axs=plt.subplots(1,3,figsize=(16,5.1));ax=axs[0];d=inc[inc['sample'].eq('absolute')&inc.estimand.eq('country_period_unweighted')&inc.event.isin(['any_increase','any_gt5','any_gt10','substantial_gt1pctP0','substantial_gt5pctP0'])].copy();names={'any_increase':'Any increase (>0)','any_gt5':'Any class >5% growth*','any_gt10':'Any class >10% growth*','substantial_gt1pctP0':'Any increment >1% of total P0','substantial_gt5pctP0':'Any increment >5% of total P0'};d['label']=[names[r.event]+f'\n({int(r.events)}/{int(r.n)})' for _,r in d.iterrows()];forest(ax,d);ax.set(xlim=(0,75),xlabel='Absolute-decoupled comparisons (%)',title='Incidence and clustered 95% intervals');lab(ax,'a')
 ax=axs[1];cats=['herbicide_kg','insecticide_kg','fungicide_kg'];names=['Herbicides','Insecticides','Fungicides /\nbactericides'];vals=[a.loc[a[c+'_delta_over_P0'].gt(0),c+'_delta_over_P0'].to_numpy()*100 for c in cats];ax.boxplot(vals,tick_labels=names,showfliers=False,patch_artist=True,boxprops={'facecolor':'#d6e8f3'},medianprops={'color':C[1]});ax.set(ylabel='Positive class increment / total P0 (%)',title='Magnitude among class-increase events');ax.tick_params(axis='x',labelsize=8);lab(ax,'b')
 ax=axs[2];d=inc[inc.estimand.eq('country_period_unweighted')&inc.event.eq('any_increase')].copy();names={'absolute':'Absolute decoupling','all_output_growth':'All output-growth periods','all_total_decline':'All total-decline periods','decline_without_output_growth':'Total decline, output not growing'};d['label']=d['sample'].map(names);forest(ax,d);ax.set(xlim=(0,100),xlabel='Comparisons with any class increase (%)',title='Common-year period comparators');lab(ax,'c');fig.tight_layout();export(fig,'Fig3_category_divergence_and_magnitude')
 # Fig4: differences in significance are not interpreted as statistically different slopes.
 co=pd.read_csv(TABLE/'economic_models_coefficients.csv');terms=['ln_gdppc0','ln_cropland0','ln_crop_intensity0','cereal_share0_10pp','vegfruit_share0_10pp'];labels=['Real GDP per capita (+10%)','Cropland area (+10%)','Real crop output / ha (+10%)','Cereal value share (+10 pp)','Fruit/vegetable value share (+10 pp)'];fig,axs=plt.subplots(1,2,figsize=(13,5.2))
 for ax,stem,title in zip(axs,['absolute','class_divergence'],['Absolute-decoupling outcome','Class increment >1% of total P0 outcome']):
  for j,(name,label,col) in enumerate([(stem+'_all_enhanced' if stem=='absolute' else 'class_divergence_all_common_enhanced','All available',C[0]),(stem+'_officialrich_enhanced','Official-rich, selected breaks removed',C[1])]):
   ds=co[co.model.eq(name)].set_index('term').reindex(terms);scale=np.array([np.log(1.1)]*3+[1,1])*100;x=ds.coefficient.to_numpy()*scale;lo=ds.ci_low.to_numpy()*scale;hi=ds.ci_high.to_numpy()*scale;y=np.arange(5)+(-.12 if j==0 else .12);ax.errorbar(x,y,xerr=np.array([x-lo,hi-x]),fmt='o' if j==0 else 's',capsize=3,color=col,label=f'{label}\nN={int(ds.nobs.iloc[0])}, countries={int(ds.clusters.iloc[0])}')
  ax.axvline(0,color='#777777',lw=1);ax.set_yticks(range(5),labels);ax.invert_yaxis();ax.set(xlabel='LPM association (percentage points; 95% interval)',title=title);ax.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.18));lab(ax,'a' if stem=='absolute' else 'b')
 fig.tight_layout();export(fig,'Fig4_economic_associations')
 # Supplements: robustness and class-specific patterns, no small-cell significance claims.
 d=sen[sen.weighting.eq('unweighted')&sen.event.eq('substantial_gt1pctP0')].copy();d['label']=d.specification.str.replace('_',' ');fig,ax=plt.subplots(figsize=(10,8));forest(ax,d);ax.set(xlim=(0,75),xlabel='Any class increment >1% of total P0 (%)');fig.tight_layout();export(fig,'FigS1_quality_window_output_sensitivity')
 long=pd.read_csv(TABLE/'common_absolute_category_magnitudes.csv');fig,axs=plt.subplots(1,2,figsize=(13,5))
 for ax,key,labtxt in zip(axs,['region','income'],['Region','Current WDI income group']):
  data=long.merge(a[['id',key]],on='id',validate='many_to_one');rows=[]
  for (grp,c),h in data.groupby([key,'category']):rows.append(dict(group=grp,category=c,estimate=h.delta_t.gt(0).mean(),n=len(h)))
  f=pd.DataFrame(rows);groups=list(f.group.unique());width=.22
  for i,(c,l,col) in enumerate(zip(cats,['Herbicides','Insecticides','Fungicides/bactericides'],C)):
   h=f[f.category.eq(c)].set_index('group').reindex(groups);ax.barh(np.arange(len(groups))+(i-1)*width,h.estimate*100,height=width,label=l,color=col)
  ax.set_yticks(range(len(groups)),groups,fontsize=8);ax.set(xlabel='Class increase in absolute-decoupled periods (%)',ylabel=labtxt,xlim=(0,100));ax.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2)
 fig.tight_layout();export(fig,'FigS2_category_patterns_by_region_income')
 # Simple methods diagram with actual sample counts.
 fig,ax=plt.subplots(figsize=(12,3.5));ax.axis('off');boxes=[('Official Phase1 snapshot','157 countries/territories\n5,251 joint country-years'),('Adjacent five-year means','893 determined comparisons\n200 absolute / 133 relative'),('Common-year H/I/F support','189 absolute comparisons\n11 excluded for missing classes'),('Magnitude and durability','Country-cluster bootstrap\nQuality/window checks + FE LPM')]
 for i,(title,desc) in enumerate(boxes):
  x=.02+i*.25;ax.text(x+.105,.55,title+'\n\n'+desc,ha='center',va='center',fontsize=10,bbox=dict(boxstyle='round,pad=.65',facecolor='#edf4f8',edgecolor=C[0]),transform=ax.transAxes)
  if i<3:ax.annotate('',xy=(x+.245,.55),xytext=(x+.22,.55),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#555555'))
 export(fig,'FigS3_methods_and_sample_flow');print('4 main and3 supplementary figures exported to PNG/SVG/PDF.',flush=True)
if __name__=='__main__':main()
