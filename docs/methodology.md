# Frozen methodology

1. Link FAO M49/area identifiers explicitly to ISO3;exclude aggregates and unbridged historical entities. Preserve original rows,units and flags. Retain1990–2024;convert RP tonnes and QV1000USD to analysis units. Preserve missingness.
2. Form adjacent,nonrolling five-year groups. At least three same-year P/output observations per group are required. Use arithmetic means and classify gP=P1/P0−1,gY=Y1/Y0−1 into absolute(gP<0,gY>0),relative(gY>0,0≤gP<gY),other,or undetermined.
3. For class events,restrict both groups to years jointly observing total P,output and all three classes. Recompute aggregate changes on that support. Compare each class increment with baseline **total** mass;the principal event is max positive class increment/P0>0.01. Own-class growth>10% uses184 positive-baseline comparisons,not189.
4. Compute equal-comparison and baseline P-mass,crop-value,cropland-weighted incidences. Cropland has188 eligible observations/79 events;other principal weights189/80. CHN exclusion and bilateral thresholds are descriptive frozen sensitivities;no new CIs were generated for them.
5. Enumerate the full study-country state matrix including unknown states. Continuation55/171 conditions on a determined next state after absolute. Adjacent runs share means;separated opportunities use different windows/denominators.
6. Replay the original country-cluster percentile bootstrap:seed20261106,4,000 draws,at least20clusters. The primary and additional-check scripts start their archived random streams separately. Do not combine or substitute intervals across tables.
7. Supporting economic models are the original country+period fixed-effects linear probability models,using baseline covariates and country-cluster finite-correction t(G−1) inference. Original model samples,coefficients,diagnostics,interactions and non-estimable specifications are preserved. No causal interpretation follows.

The archive reorganises file paths only. Frozen outputs are read-only comparators;all replays write .reproduction/. No new specification,model selection,parameter search or figure redesign is part of archive preparation.
