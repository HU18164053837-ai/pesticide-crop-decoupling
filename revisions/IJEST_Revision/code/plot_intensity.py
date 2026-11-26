"""One restrained two-panel main figure,from verified same-year extension tables."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
R=Path(__file__).resolve().parents[1];D=pd.read_csv(R/'02_TABLES/intensity_cluster_estimates.csv');BLUE='#23679C';GRAY='#A5ADB4';INK='#30363B'
plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.labelsize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#AAB2B9','axes.linewidth':.65,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
def one(cohort,event):
    d=D[D.cohort.eq(cohort)&D.event.eq(event)];assert len(d)==1;return d.iloc[0]
f,axes=plt.subplots(1,2,figsize=(174/25.4,96/25.4));f.subplots_adjust(left=.085,right=.98,bottom=.24,top=.89,wspace=.40)
source=[]
a=axes[0];rows=[one('frozen_'+s,'intensity_up') for s in ['absolute','relative','coupled_or_other']];xx=np.arange(3)
for i,r in enumerate(rows):
    a.bar(i,r.estimate*100,width=.55,color=BLUE,zorder=2);a.errorbar(i,r.estimate*100,yerr=[[100*(r.estimate-r.ci_low)],[100*(r.ci_high-r.estimate)]],fmt='none',ecolor=INK,capsize=3,elinewidth=.85)
    a.text(i,r.ci_high*100+3,f'{r.estimate*100:.1f}',ha='center',va='bottom',fontsize=8.5)
    source.append(dict(panel='a',category=['Absolute','Relative','Other'][i],**r.to_dict()))
a.set_xticks(xx,[l+f'\n{int(r.events)}/{int(r.n)}' for l,r in zip(['Absolute','Relative','Other'],rows)]);a.set(ylabel='Comparisons with rising kg/ha (%)',ylim=(0,100));a.set_yticks([0,25,50,75,100]);a.text(-.20,1.08,'a',transform=a.transAxes,fontweight='bold',fontsize=10)
b=axes[1];cats=['herbicide_kg','insecticide_kg','fungicide_kg']
for j,(event,label,color) in enumerate([('mass_up','Mass',GRAY),('intensity_up','kg/ha',BLUE)]):
    rr=[one('class_land_common_absolute',c+'_'+event) for c in cats];y=np.array([r.estimate*100 for r in rr]);lo=np.array([r.ci_low*100 for r in rr]);hi=np.array([r.ci_high*100 for r in rr]);pos=xx+(j-.5)*.32
    b.bar(pos,y,width=.29,color=color,label=label,zorder=2);b.errorbar(pos,y,yerr=[y-lo,hi-y],fmt='none',ecolor=INK,capsize=2,elinewidth=.8)
    for c,r in zip(cats,rr):source.append(dict(panel='b',category=c,**r.to_dict()))
b.set_xticks(xx,['Herbicides','Insecticides','Fungicides /\nbactericides']);b.set(ylabel='Class increase frequency (%)',ylim=(0,45));b.set_yticks([0,15,30,45]);b.legend(loc='upper center',bbox_to_anchor=(.5,1.16),ncol=2,fontsize=8);b.text(-.20,1.08,'b',transform=b.transAxes,fontweight='bold',fontsize=10)
for ax in axes:ax.tick_params(length=3,width=.6);ax.set_axisbelow(True)
f.canvas.draw();renderer=f.canvas.get_renderer();bad=[]
for t in f.findobj(Text):
    if not t.get_visible() or not t.get_text().strip():continue
    z=t.get_window_extent(renderer)
    if z.x0<0 or z.y0<0 or z.x1>f.bbox.width or z.y1>f.bbox.height:bad.append(t.get_text())
assert not bad,bad
f.savefig(R/'03_FIGURES/Fig4_land_intensity.svg');f.savefig(R/'03_FIGURES/Fig4_land_intensity.pdf');f.savefig(R/'03_FIGURES/Fig4_land_intensity.png',dpi=600)
f.savefig(R/'03_FIGURES/Fig4_land_intensity.eps');f.savefig(R/'03_FIGURES/Fig4_land_intensity.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
pd.DataFrame(source).to_csv(R/'03_FIGURES/Fig4_SOURCE_VALUES.csv',index=False,encoding='utf-8-sig');(R/'07_QA/FIG4_LAYOUT.json').write_text(json.dumps(dict(width_mm=174,height_mm=96,min_font=8,no_clipping=True)))
print('Fig4 exported:SVG/vectorPDF/600dpiPNG;two panels,no annotation blocks.')
