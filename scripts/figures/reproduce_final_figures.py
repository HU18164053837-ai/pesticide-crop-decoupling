"""Optional original figure replay;not run during archive preparation."""
from pathlib import Path
import sys,subprocess,shutil
R=Path(__file__).resolve().parents[2];HERE=R/'.reproduction/figure_recipes';P2=R/'.reproduction/phase2'
if not (P2/'tables/economic_models_coefficients.csv').exists():raise SystemExit('Run scripts/run_all.py --with-inference first.')
for d in ['FIGURES_MAIN','FIGURES_SI']:(HERE/d).mkdir(parents=True,exist_ok=True)
subprocess.run([sys.executable,'-I',str(Path(__file__).with_name('phase2_figures.py'))],check=True)
subprocess.run([sys.executable,'-I',str(Path(__file__).with_name('final_main_recipes.py'))],check=True)
sys.path.insert(0,str(Path(__file__).resolve().parent))
from refine_figure2 import draw
draw(HERE,P2)
from layout_sample_flow import draw as flow
flow(HERE)
for original,final in [('Fig4_economic_associations','FigS2_economic_associations'),('FigS1_quality_window_output_sensitivity','FigS3_quality_and_window_sensitivity'),('FigS2_category_patterns_by_region_income','FigS4_region_and_income'),('Fig3_category_divergence_and_magnitude','FigS5_class_magnitudes_and_comparators')]:
    for ext in ['svg','png','pdf']:shutil.copy2(P2/'figures'/f'{original}.{ext}',HERE/'FIGURES_SI'/f'{final}.{ext}')
print('Original recipes replayed in .reproduction/figure_recipes;frozen figures untouched.')
