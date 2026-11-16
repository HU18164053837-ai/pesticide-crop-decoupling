# Actual reproduction procedure

Tested with Python3.12.14 and the pinned requirements. A new venv was created without access to user site-packages;the already frozen dependency distributions were physically copied into its site-packages. Python isolated mode(-I) was used and confirmed that scientific imports resolved within the venv. No external PYTHONPATH or live data download was needed. This tests an independent environment using the historical dependency versions;it does not claim a fresh package-index installation was tested.

Run `python scripts/run_all.py --with-inference` from any working directory. Order:

1. prepare_frozen.py validates RP/QV units and original-row values,reconstructs measured fields,and runs the original comparisons() algorithm.
2. recompute_targeted.py independently enumerates original-row common supports and repeats the original boundary,CHN-weight and concentration checks.
3. analyze.py replays original primary bootstrap,transition/persistence,sensitivities and all15 model specifications.
4. additional_checks.py replays the original quality-tier bootstrap,exploratory interaction checks and weighting tables.
5. validate_archive.py compares recomputed/frozen tables,checks all requested counts/rates,verifies CI ordering and confirms source/figure hashes.

Numerical table comparison tolerances are rtol1e−7/atol1e−9 for replayed floating-point fields;IDs,categories,columns and counts must match exactly. No frozen interval endpoint is modified. Actual per-check results and environment metadata are in validation_results.csv and validation_summary.json.

The original full Phase1 bulk-database retrieval/cleaning/model/figure routine is archived as scripts/data_preparation/phase1_pipeline.py but **not executed** by the offline entry point. It needs the historical RP/QV/RL/QI bulk ZIPs and World Bank JSONs under its --root path;these large source databases are deliberately not uploaded. Their exact URLs,dates and SHA-256 values are in data/metadata/source_register.csv. Current live versions may differ. Source extracts and the annual panel allow the principal analysis to be reproduced without recovering those full bulk databases. The Phase1 historical raw-bulk replay is therefore not claimed verified in this archive exercise.

Figure recipes are supplied for transparency,but no plotting script is run during acceptance. Original v1.4SVG/PNG bytes are compared with provenance hashes. Visual/source-data traceability follows figure_source_map.csv;the map's historical archive paths describe lineage,not missing runtime dependencies. Current executable locations are scripts/figures/.
