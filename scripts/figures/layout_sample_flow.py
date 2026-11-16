"""Necessary layout repair of frozen FigureS6; static original labels and counts."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
def draw(here):
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none','pdf.fonttype':42,'savefig.dpi':300})
 fig,ax=plt.subplots(figsize=(8.8,4.4));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
 boxes=[(.03,.55,'1  Official Phase1 snapshot','157 countries/territories\n5,251 joint country-years'),(.54,.55,'2  Adjacent five-year means','893 determined comparisons\n200 absolute / 133 relative'),(.54,.05,'3  Common-year H/I/F support','189 absolute comparisons\n11 excluded for missing classes'),(.03,.05,'4  Magnitude and durability','Country-cluster bootstrap\nQuality/window checks + FE LPM')]
 for x,y,title,desc in boxes:
  ax.add_patch(FancyBboxPatch((x,y),.43,.37,boxstyle='round,pad=0.008',facecolor='#edf4f8',edgecolor='#0072B2',linewidth=1))
  ax.text(x+.215,y+.27,title,ha='center',va='center',fontsize=10.5,fontweight='bold');ax.text(x+.215,y+.13,desc,ha='center',va='center',fontsize=11,linespacing=1.4)
 for start,end in [((.47,.735),(.53,.735)),((.755,.54),(.755,.43)),((.53,.235),(.47,.235))]:ax.annotate('',xy=end,xytext=start,arrowprops=dict(arrowstyle='->',color='#555555',lw=1.3))
 name='FigS6_methods_and_flow'
 for ext in ['png','svg','pdf']:fig.savefig(here/'FIGURES_SI'/(name+'.'+ext),bbox_inches='tight',pad_inches=.12)
 plt.close(fig)
 return here/'FIGURES_SI'/(name+'.png')
