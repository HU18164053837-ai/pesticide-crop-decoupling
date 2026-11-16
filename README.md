# Pesticide–crop output decoupling

Research archive for **Beyond aggregate pesticide decoupling: state recurrence and class-specific input changes across 157 countries**.

This repository preserves the frozen analysis supporting the EDS submission. It studies whether declining national pesticide mass alongside increasing real crop output recurs across five-year comparisons, and whether individual pesticide classes can increase despite aggregate decline. The analysis is descriptive and observational; active-ingredient mass is not a measure of toxicity, exposure or environmental benefit.

## Frozen scope

- Data: 1990–2024; 236 mapped country/territory entities,157 with joint pesticide/output observations;5,251 joint country-years.
- Five-year means: nonrolling groups from1990–1994 through2020–2024,with at least three common valid annual observations per group.
- Final figure assets: original EDS v1.4 figures only. Supplementary tables: v1.4 plus the v1.5 Table S15 quality-tier continuation table. Later graphical experiments are excluded.
- No manuscript or historical ZIP bundle is distributed here. The paper-table representations in results/tables/manuscript/ preserve its tabulated values.

## Data sources and variables

[FAOSTAT](https://www.fao.org/faostat/en/#data) domains: Pesticides Use(RP),Value of Agricultural Production(QV),Land Use(RL) and Production Indices(QI). Total pesticide mass is reported item1357,not a sum of selected classes. Herbicides,insecticides and fungicides/bactericides are separate reported headings. Tonnes of active ingredient are converted to kg; QV item2041/element58 is converted from1000 USD at constant2014–2016 prices to USD. Historical variable names ending usd2015 are retained but do not change that QV price basis.

[World Bank World Development Indicators](https://databank.worldbank.org/source/world-development-indicators): NY.GDP.PCAP.KD,NV.AGR.TOTL.KD and AG.PRD.CROP.XD. WDI also supplies descriptive region/current income classifications. Cropland is arable plus permanent-crop area,requiring both components. Full definitions,units,flags,codes,source update dates,download URLs and original hashes are retained in data/metadata/ and docs/data_dictionary.md.

Data were retrieved on6November2026. Frozen extracts,not current API responses,are the reproducibility inputs. No interpolation or historical-territory bridging is performed.

## Definitions

For adjacent group means, gP=P1/P0−1 and gY=Y1/Y0−1. Absolute decoupling requires gP<0 and gY>0;relative decoupling requires gY>0 and0≤gP<gY. Other determined comparisons are coupled_or_other. Missing support yields undetermined.

Class comparisons require pesticide total,crop output and all three major pesticide classes observed on the same retained years in each group. The principal class event is at least one class-mass increment exceeding1% of baseline total pesticide mass P0. A class share increase alone is not sufficient. Own-class percentage-growth events require positive class baselines and therefore have a different denominator. Weighted estimates use baseline masses,real crop values or cropland areas;they answer different aggregation questions.

Adjacent changes share an intermediate five-year mean. Next-state continuation and country-ever run frequencies have different denominators and are not independent transition probabilities.

## Run the analysis

Python3.12 was used; exact tested library versions are in requirements.txt and reproducibility/validation_summary.json. Create an isolated environment,install dependencies,and run:

```bash
python -m venv .venv
# Activate .venv using the command appropriate to your operating system.
python -m pip install -r requirements.txt
python scripts/run_all.py --with-inference
```

The default `python scripts/run_all.py` rebuilds original RP/QV measurements,comparison support and deterministic boundary/CHN sensitivities. `--with-inference` additionally replays the **original** seed20261106,4,000-replicate country-cluster bootstrap,15 pre-existing fixed-effects model specifications,and exploratory quality interactions. This is verification of frozen methods,not specification search. All fresh computations go to ignored .reproduction/;frozen results/ is never overwritten.

See reproducibility/run_all.md for order,environment provenance and precise validation scope. No database retrieval,figure redraw or manuscript generation occurs in this entry point. Original figure-generating recipes are retained in scripts/figures/;the nine final figures are reused byte-for-byte from v1.4 and checked against provenance hashes.

## Main frozen findings

| Result | Value |
|---|---:|
| Reporting countries/territories |157|
| Joint country-years |5,251|
| Determined five-year comparisons |893|
| Absolute / relative decoupling |200 /133|
| Common-supported absolute comparisons |189|
| At least one class increases |102/189|
| A class increment exceeds1% of baseline total mass |80/189|
| Bilateral1% total-change sensitivity |72/172|
| Next determined state remains absolute |55/171|
| Class event after decline without output growth |35/122|
| Baseline crop-value weighted incidence |15.63%|
| Crop-value weighted incidence excluding CHN |31.16%|

Class-specific>1%-of-P0 events are42/189(herbicides),22/189(insecticides) and30/189(fungicides/bactericides);classes can overlap. Counts and country numbers refer to the relevant support,not one universal denominator.

## Uncertainty and limitations

Stored confidence intervals are retained with their exact originating table. Table S15's all-record continuation CI is24.683544–39.585129%,whereas the archived state-continuation table gives24.683544–39.325843%. The point estimate55/171 is identical. These separately generated frozen bootstrap endpoints are **not harmonised**;see reproducibility/ci_endpoint_provenance.csv.

National reporting may combine surveys,sales/import proxies,estimates and imputations. A/E/I/X flags are reporting categories,not environmental-risk grades. Source notes,known breaks,common-year restrictions,nonpositive class baselines and non-additivity are documented and retained. Present-day WDI income/region groups are not historical assignments. Size weights are concentrated in some countries. Overlapping comparisons,measurement heterogeneity and ecological aggregation limit interpretation. The analysis does not identify causal mechanisms or pesticide risk.

## Licence,attribution and anonymous review

Original code is MIT-licensed;this does **not** license third-party data. FAO and World Bank data retain their CC BY4.0 plus provider-specific additional terms,attribution requirements and applicable exceptions. See data/LICENSE.md and docs/DATA_REDISTRIBUTION_AUDIT.md. Do not infer provider endorsement.

The working GitHub archive is private. GitHub ownership and commit metadata are not an anonymous-review solution. A separate identity-scrubbed ZIP without .git is prepared for submission as a reviewer-accessible supplementary archive. No identifying repository link is inserted into the anonymous manuscript. See docs/ANONYMOUS_REVIEW.md.

## Repository map

data/ contains frozen original-row extracts,analytic panels and metadata;scripts/ preserves preparation,classification,inference and figure recipes;results/ contains frozen tables and original figures;reproducibility/ contains replay checks,hash manifests and instructions;docs/ explains methods,licensing and review access.
