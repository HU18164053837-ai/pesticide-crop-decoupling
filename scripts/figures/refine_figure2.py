"""Readability-only revision, using frozen Phase2 points and confidence intervals."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def draw(here,p2):
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'svg.fonttype':'none','pdf.fonttype':42})
 fig=plt.figure(figsize=(9.6,8.1));axes=[fig.add_axes([.09,.57,.32,.34]),fig.add_axes([.76,.57,.22,.34]),fig.add_axes([.35,.12,.63,.29])]
 tr=pd.read_csv(p2/'tables/transition_matrix_study_countries.csv');states=['absolute','relative','coupled_or_other','undetermined'];names=['Absolute','Relative','Other','Unknown']
 prob=tr.pivot(index='from_state',columns='to_state',values='row_probability').reindex(index=states,columns=states).fillna(0);ct=tr.pivot(index='from_state',columns='to_state',values='n').reindex(index=states,columns=states).fillna(0)
 ax=axes[0];ax.imshow(prob,cmap='Blues',vmin=0,vmax=1);ax.set_xticks(range(4),names,rotation=35,ha='right');ax.set_yticks(range(4),names);ax.set(title='Complete transition grid',xlabel='Next comparison',ylabel='Current comparison')
 for i in range(4):
  for j in range(4):ax.text(j,i,f'{prob.iloc[i,j]*100:.0f}%\nn={int(ct.iloc[i,j])}',ha='center',va='center',fontsize=9.5,color='white' if prob.iloc[i,j]>.6 else 'black')
 d=pd.read_csv(p2/'tables/absolute_continuation_by_quality.csv');labels={'all':'All records','officialrich':'Official-rich countries','no_known_break':'Selected breaks removed','officialrich_no_known_break':'Official-rich,\nbreaks removed','period_A_majority_no_known_break':'Period A-majority,\nbreaks removed'};ax=axes[1];y=np.arange(len(d));x=d.estimate.to_numpy()*100;lo=d.ci_low.to_numpy()*100;hi=d.ci_high.to_numpy()*100
 ax.errorbar(x,y,xerr=np.array([x-lo,hi-x]),fmt='o',color='#0072B2',capsize=3);ax.set_yticks(y,[labels[r['sample']]+f'\n{int(r.events)}/{int(r.n)}' for _,r in d.iterrows()],fontsize=10.5);ax.invert_yaxis();ax.set(xlim=(0,100),title='Next determined state\nremains absolute',xlabel='Continuation (%)');ax.set_xticks([0,25,50,75,100])
 d=pd.read_csv(p2/'tables/absolute_persistence_summary.csv');ax=axes[2];ax.barh(range(len(d)),d.rate*100,color=['#0072B2','#0072B2','#009E73','#CC79A7']);labs=[]
 for _,r in d.iterrows():
  nm={'adjacent_changes':'Adjacent','separated_pairs_offset0':'Separated offset 0','separated_pairs_offset1':'Separated offset 1'}[r.design]
  labs.append(f'{nm}, run ≥{int(r.min_run)}\n{int(r.countries_with_run)}/{int(r.eligible_countries)} countries; {int(r.eligible_windows)} windows')
 ax.set_yticks(range(len(d)),labs,fontsize=10.5);ax.invert_yaxis();ax.set(xlim=(0,34),title='Different opportunity grids',xlabel='Eligible countries ever meeting run (%)')
 for ax,l in zip(axes,'abc'):ax.text(-.03,1.08,l,transform=ax.transAxes,fontweight='bold',fontsize=12)
 name='Fig2_continuation_and_opportunities'
 for ext in ['png','svg','pdf']:fig.savefig(here/'FIGURES_MAIN'/(name+'.'+ext),bbox_inches='tight',pad_inches=.15)
 plt.close(fig)
 fm=pd.read_csv(here/'metadata/FIGURE_SOURCE_MAP.csv');ix=fm.figure.eq(name);fm.loc[ix,'operation']='Phase3.1 readability-only layout; frozen counts, points and intervals unchanged';fm.loc[ix,'generator']='src/refine_figure2.py';fm.to_csv(here/'metadata/FIGURE_SOURCE_MAP.csv',index=False,encoding='utf-8-sig')
 return here/'FIGURES_MAIN'/(name+'.png')
