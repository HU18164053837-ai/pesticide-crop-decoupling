# Reproduction acceptance

**REPRODUCIBILITY_VERIFIED** for the final analytic results,with the scope below.

The existing analysis was replayed in a new Python3.12.14 venv using isolated mode. A second run was completed from a clean directory containing only the curated archive files. The second run required neither the original project paths nor the historical ZIP packages.

The736 initial checks expanded to748 after including exact high-resolution v1.4PDF assets and upstream indicator metadata. The final record is validation_results.csv:748 checks,zero failures. Checks include reconstructed sample counts,all requested principal event/weight results,the study-country state-transition matrix,36 replayed inference/diagnostic tables,and source/figure integrity. Floating-point fields were compared at rtol1e−7/atol1e−9;categorical values and sample keys were compared exactly. These are executed assertions,not independent statistical tests.

The replay used the original4,000-draw bootstrap seed and call order,the original15 fixed-effects model specifications,and existing quality-interaction checks. It did not add models,choose specifications or modify frozen outputs. All original CI endpoints remain unchanged. Table S15 and the state-continuation table retain their different archived upper endpoints,documented in ci_endpoint_provenance.csv.

Separate source-document checks passed211 assertions:all26 existing v1.4paper/SI tables were exported cell-for-cell;all9embedded PNGs matched the archived originals byte-for-byte;main Table2 conditional rates,the15 model-table coefficients/intervals,and strict-quality39/17 strata were checked against source outputs. See document_consistency_checks.csv.

## Explicit limits

- The full historical Phase1 live bulk-database retrieval and exploratory modelling routine was not rerun. Its code,source URLs,version dates,selected original rows and raw checksums are preserved;the offline final analysis starts from frozen observations and covariates. No fresh database is substituted.
- A new package-index installation was not tested. The independent venv contains physically copied historical dependency distributions;versions and this provisioning method are recorded.
- Figure-generating recipes were not executed in this archive task,as requested. Frozen original v1.4PNG/SVG/PDF assets were verified by hash and the PNGs against the final Word documents. Optional original-recipe replay is documented below;pixel-level cross-platform rendering reproducibility is not claimed.
- The original manuscript files are not uploaded and have not been edited. An anonymous attachment is prepared,but no journal submission or hosted anonymous-view URL is created by this task.

## Optional figure replay after review

After the statistical replay,`python scripts/figures/reproduce_final_figures.py` runs the existing recipes into an ignored output directory. It does not overwrite results/figures/. This optional command was not run during archive preparation. Different font/rendering/library builds can affect pixels,SVG IDs and PDF timestamps even when the underlying data agree.
