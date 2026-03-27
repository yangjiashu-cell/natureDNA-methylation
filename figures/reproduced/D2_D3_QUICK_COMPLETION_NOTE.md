## D2 / D3 Quick Completion (reuse-first)

This note records the fast-track completion strategy requested by the user:
reuse existing author/public data and existing computed outputs, avoid long reruns.

### D2 (PISCA methylation phylogeny)

- Kept primary Zenodo-based PISCA outputs in `figures/reproduced/pisca_zenodo_fig4/`:
  - `fig4a_case12_zenodo.xml/.log/.trees/_MCC.tree`
  - `fig4b_case19_zenodo.xml/.log/.trees/_MCC.tree`
  - `fig4a_case12_zenodo_MCC_ggtree.pdf`
  - `fig4b_case19_zenodo_MCC_ggtree.pdf`
- Removed extra/non-primary PISCA files to reduce clutter:
  - all files under `figures/reproduced/pisca_moesm8_fallback/`
  - redundant state files:
    - `fig4a_case12_zenodo_states_beta_em.csv`
    - `fig4a_case12_zenodo_states_beta_mixture.csv`

### D3 (SNV / BEAST 1.10)

- Reused existing prepared outputs:
  - `figures/reproduced/snv_beast10_prep/step13_preflight_report.json`
  - `figures/reproduced/snv_beast10_prep/case_*_binary_matrix.csv`
  - `figures/reproduced/snv_beast_from_moesm3/alignment_binary.nexus`
- BEAST 1.10.4 installed and detected:
  - `tools/BEAST v1.10.4/bin/beast.cmd`
  - reported in `step13_preflight_report.json` as `beast10_found`.

### Scope statement

- This quick-completion pass is designed for "run-through" reproducibility with available data/results.
- It does **not** claim full paper-scale recomputation (e.g., exhaustive long-chain reruns for every step).

