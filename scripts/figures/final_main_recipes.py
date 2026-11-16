"""Original final Fig1/Fig3/S1 recipe; optional only, not executed by validation."""
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
REPO=Path(__file__).resolve().parents[2];HERE=REPO/'.reproduction/figure_recipes';P1=REPO/'.reproduction/phase1';T=REPO/'results/tables'
for d in ['FIGURES_MAIN','FIGURES_SI']:(HERE/d).mkdir(parents=True,exist_ok=True)
def read(name):return pd.read_csv(T/name)
annual=pd.read_csv(REPO/'data/processed/annual_reverified.csv');inc=read('divergence_incidence_bootstrap.csv');weights=read('agricultural_scale_weighted_incidence.csv');per=read('absolute_persistence_summary.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'svg.fonttype':'none','pdf.fonttype':42})
C=['#0072B2','#D55E00','#009E73','#CC79A7']
fig_sources=[]
def export(fig,name,folder,sources,operation):
 for ext in ['png','svg','pdf']:fig.savefig(HERE/folder/(name+'.'+ext),bbox_inches='tight',pad_inches=.15)
 plt.close(fig);fig_sources.append(dict(figure=name,folder=folder,source_tables=sources,operation=operation,generator='src/build_phase3.py'))
def label(ax,t):ax.text(-.03,1.08,t,transform=ax.transAxes,fontweight='bold',fontsize=12)
def forest(ax,d,labels):
 x=d.estimate.to_numpy()*100;lo=d.ci_low.to_numpy()*100;hi=d.ci_high.to_numpy()*100;y=np.arange(len(d))
 ax.errorbar(x,y,xerr=np.array([x-lo,hi-x]),fmt='o',color=C[0],capsize=3);ax.set_yticks(y,labels);ax.invert_yaxis()
w=pd.read_csv(P1/'results/tables/official_world_and_matched_country_totals.csv');bal=pd.read_csv(P1/'results/tables/balanced_global_trends.csv')
fig=plt.figure(figsize=(9.6,7.6));gs=fig.add_gridspec(2,2);axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
for v,l,c in [('FAO_world_p_kg','Pesticide mass',C[0]),('FAO_world_crop_usd2015','Real crop gross output',C[1])]:axs[0].plot(w.year,w[v]/w.loc[w.year.eq(1990),v].iloc[0]*100,label=l,color=c,lw=2)
axs[0].set(title='FAOSTAT World aggregates',xlabel='Year',ylabel='Index (1990 = 100)');axs[0].legend(frameon=False,fontsize=8)
d=bal[bal.base_year.eq(1990)]
for v,l,c in [('p_kg_index','Pesticide mass',C[0]),('crop_usd2015_index','Real crop gross output',C[1])]:axs[1].plot(d.year,d[v],label=l,color=c,lw=2)
axs[1].set(title='Fixed 116-entity balanced subset',xlabel='Year',ylabel='Index (1990 = 100)');axs[1].legend(frameon=False,fontsize=8)
f=read('source_flags_by_year.csv').set_index('year');axs[2].stackplot(f.index,*[f[v] for v in ['A','E','I','X']],labels=['Official (A)','Estimated (E)','Imputed (I)','External (X)'],colors=[C[2],C[0],C[1],C[3]])
axs[2].set(title='Total pesticide source flags',xlabel='Year',ylabel='Reporting entities',ylim=(0,250));axs[2].legend(frameon=False,fontsize=8,ncol=2,loc='upper left')
for ax,l in zip(axs,'abc'):label(ax,l)
fig.tight_layout();export(fig,'Fig1_trends_and_reporting','FIGURES_MAIN','official_world_and_matched_country_totals.csv; balanced_global_trends.csv; source_flags_by_year.csv','World indices rescaled to 1990=100; fixed-subset stored indices; flags plotted as stored')

fig=plt.figure(figsize=(9.6,8));gs=fig.add_gridspec(2,2);axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[1,:]),fig.add_subplot(gs[0,1])]
events=['any_increase','substantial_gt1pctP0','substantial_gt5pctP0','any_gt10'];d=inc[(inc['sample']=='absolute')&(inc.estimand=='country_period_unweighted')].set_index('event').loc[events].reset_index()
forest(axs[0],d,[f'Any class increase\n102/189',f'Increment >1% of total P₀\n80/189',f'Increment >5% of total P₀\n37/189',f'Class growth >10%*\n68/184']);axs[0].set(title='Event definitions and denominators',xlabel='Absolute comparisons (%)',xlim=(0,75))
u=d[d.event.eq('substantial_gt1pctP0')].iloc[0];wd=weights[(weights['sample']=='all_common_absolute')&(weights.event=='substantial_gt1pctP0')].set_index('weight').loc[['cp0','cy0','cropland_ha0']].reset_index();dd=pd.concat([pd.DataFrame([u]),wd],ignore_index=True)
forest(axs[1],dd,['Equal comparison weights (n=189)','Baseline pesticide mass (n=189)','Baseline real crop value (n=189)','Baseline cropland (n=188)']);axs[1].set(title='Same >1% event, different estimands',xlabel='Weighted incidence (%)',xlim=(0,65))
ca=read('category_incidence_bootstrap.csv');dd=ca[ca.event.eq('class_gt1pctP0')].set_index('category').loc[['herbicide_kg','insecticide_kg','fungicide_kg']].reset_index();forest(axs[2],dd,[f'Herbicides\n42/189',f'Insecticides\n22/189',f'Fungicides / bactericides\n30/189']);axs[2].set(title='Which class exceeds 1% of total P₀?',xlabel='Absolute comparisons (%)',xlim=(0,45))
for ax,l in zip(axs,'abc'):label(ax,l)
fig.tight_layout(w_pad=2);export(fig,'Fig3_class_thresholds_and_weights','FIGURES_MAIN','divergence_incidence_bootstrap.csv; agricultural_scale_weighted_incidence.csv; category_incidence_bootstrap.csv','Stored canonical point estimates and clustered 95% intervals; no bootstrap rerun')

fig,ax=plt.subplots(figsize=(10,8));h=annual.pivot(index='iso3',columns='year',values='joint_py').sort_index();ax.imshow(h,aspect='auto',cmap='Blues',interpolation='none',vmin=0,vmax=1);ax.set_xticks(range(0,35,5),range(1990,2025,5));ax.set_yticks(range(0,len(h),10),h.index[::10],fontsize=7);ax.set(xlabel='Year',ylabel='ISO3 in full mapped grid',title='Joint pesticide and crop-output availability: 236 mapped entities');fig.tight_layout();export(fig,'FigS1_full_coverage_grid','FIGURES_SI','annual_reverified.csv','Binary joint_py availability; mapped universe differs from 157 observed study entities')