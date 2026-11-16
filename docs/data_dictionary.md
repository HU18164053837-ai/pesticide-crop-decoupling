# Data dictionary

The exhaustive annual-field dictionary is data/metadata/data_dictionary.csv;comparison/outcome patterns are in PHASE2_DATA_DICTIONARY.csv. The source register gives domain,element,item,unit,update date,URL and raw-file hash for each measured variable.

| Field | Meaning / unit |
|---|---|
|iso3,year|Country/territory ISO3 and calendar year;crosswalk is explicit|
|p_kg|FAOSTAT total pesticide use,item1357/element5157;kg active ingredient,reported tonnes×1000|
|herbicide_kg,insecticide_kg,fungicide_kg|Three separate major pesticide headings;kg active ingredient|
|crop_usd2015|Historical variable name for QV2041/58 crop gross output at constant2014–2016 prices;source1000USD×1000|
|cropland_ha|Arable+permanent-crop area;both required|
|gdp_pc_const|WDI real GDP per capita,constant2015USD|
|qi_crop,wb_crop_index|Alternative crop-output indices;2014–2016=100|
|*_flag|Original FAOSTAT A/E/I/X reporting provenance,not risk grades|
|joint_py|Total pesticide observed and real crop output strictly positive|
|id|ISO3_previous-start_current-start_window-width comparison key|
|p_kg0,p_kg1,y0,y1|Baseline/current five-year means on joint P/output years|
|gp,gy,state|Relative mean changes and aggregate decoupling classification|
|common_support0,common_support1|Exact annual years supporting total/output/three classes jointly|
|cp0,cp1,cy0,cy1,cgp,cgy|Corresponding common-year mass/output means and changes|
|*_delta_over_P0|Class-mass increment divided by common-year baseline total pesticide mass|
|any_increase|At least one class increases strictly above zero|
|substantial_gt1pctP0|At least one class increment exceeds1% of baseline total mass|
|any_gt10|Own-class growth exceeds10%,requiring positive baselines for all classes|
|estimate,ci_low,ci_high|Frozen point/interval;fraction for incidences,but original coefficient units in model tables|
|events,n,countries|Event comparisons,eligible comparisons and country clusters for that row|

Missing values are not interpolated. Class shares are distinct from class quantities. Modern WDI groups describe the archived classifications,not historical income status.
