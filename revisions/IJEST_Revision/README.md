# IJEST land-intensity extension candidate

Status: **CANDIDATE FOR AUTHOR REVIEW — NOT FORMALLY FROZEN**. No release/tag or merge to the EDS baseline is performed.

The unchanged EDS data and original results are retained at commit d959d7eea039290d312663395907108edf4bf800. This branch adds a bounded, common-year land-normalized analysis:18/199 absolute comparisons have falling total pesticide mass,growing crop output and rising kg/ha;country-cluster95% interval5.2–13.3%. Reporting quality matters:the official-rich/no-known-break subset gives2/67. The original200/893,55/171 and80/189 definitions/estimates remain unchanged.

Data use is limited to the existing root data/processed and data/source_extracts. P/Y/positive cropland must share at least3observed years in each5-year group. Intensity is ratio of period means,not mean of annual ratios. Class intensities additionally require all three class records on identical years. All887retained comparisons preserve their frozen states;188/189class comparisons remain. No new model or threshold search is performed.

Using the repository's pinned Python environment,run from the repository root:

```
python revisions/IJEST_Revision/code/land_intensity.py
python revisions/IJEST_Revision/code/validate_intensity.py
python revisions/IJEST_Revision/code/plot_intensity.py
```

The scripts write only inside this revision directory;original data/results are read-only. Bootstrap uses4000whole-country draws and deterministic cohort/event seeds. New intervals must not replace any frozen EDS intervals. The pooled arable/permanent-crop A-share screen with16countries has noCI. Results are national input descriptors,not application doses,concentrations,exposure or risk.

01_ANALYSIS stores all record-level support and annual-observation audits;02_TABLES stores estimates,support loss and event IDs;03_FIGURES contains only the newFig.4;06_LITERATURE documents ordinary established indicators and limited empirical originality;07_QA stores executed validation. The original9figures are not redrawn. The candidate manuscript,Title Page,Cover Letter and personal author details are deliberately not included in this public branch. A separate anonymous reviewer ZIP is assembled for submission.

Code/data licences remain those of the root repository. Author approval is required before a final IJEST freeze.
