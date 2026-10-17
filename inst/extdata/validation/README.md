# Calibration results

Generated on 2026-10-06 by inst/scripts/10-calibration.R with 4 workers,
1000 null replicates per sample-label test, 100 per node-label test,
100 (node-label: 20) replicates per power cell, in 192 minutes.
Settings as in the paper's analyses: signal_mode logfc, logfc_k 0.7,
resolution 30, smoothing 40, n_perm 999 or exact enumeration, 100 TFCE steps.

* type1_error.csv: family-wise error rate of each test under the global
  null (no effect anywhere), by network, with its SE and 95% Wilson
  interval over `replicates` null datasets.
* null_pvalues_*.csv: smallest p-value per null replicate, for QQ plots.
* power.csv: probability (and SE) of at least one significant result
  when an 8-gene module, the hub module or a peripheral path, is
  shifted by `effect` in 4 vs 4 samples; three tests on the same data.
* layout_sensitivity*.csv: the same data on twenty Fruchterman-Reingold
  layouts of the 300-node graph; pairwise Jaccard of the genes in
  significant regions, versus the layout-free TFCE result.

levi 1.99.0, R 4.6.1.

Denominators corrected on 2026-10-17: datasets with no detected region are
non-rejections and remain in the FWER and power denominators. The original
simulation draws were retained. Type I summaries were recomputed from the
saved null p-values (1000 datasets per sample-label test). Power rejection
counts were recovered from the saved rates and denominators, using the
original 100 datasets per sample-label/TFCE cell and 20 per node-label cell.
The ba300 path/effect-2 sample-label power is 4/100, rather than 4/99.
NA in null_pvalues_*.csv denotes no detected region, not a failed simulation.
