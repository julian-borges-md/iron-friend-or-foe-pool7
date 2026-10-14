# Iron Overload and Male Hypogonadism: Systematic Review and Meta-Analysis

**Iron: Friend or Foe? Pool 7 (Male Hypogonadism)**

Borges JYV.

PROSPERO registration: [CRD420261344877](https://www.crd.york.ac.uk/prospero/display_record.php?RecordID=1344877)

## Scope

Systematic review of the association between iron overload (hereditary haemochromatosis, transfusion-dependent thalassaemia, dysmetabolic iron overload syndrome) and male hypogonadism, with a random-effects meta-analysis of hypogonadism prevalence in transfusion-dependent thalassaemia (TDT) males.

## Repository layout

| Path | Content |
|---|---|
| `data/tdt_male_hypogonadism_prevalence.csv` | Study-level extraction for the six TDT studies entered into the meta-analysis |
| `analysis/meta_analysis.R` | Primary analysis (`meta::metaprop`, Freeman-Tukey double arcsine, DerSimonian-Laird), forest plot, leave-one-out sensitivity, `metafor` cross-check |
| `analysis/meta_analysis.py` | Independent Python implementation for numerical verification |
| `figures/` | PRISMA 2020 flow diagram and forest plot outputs |
| `protocol/` | Search strategy, eligibility criteria, risk of bias and GRADE procedures |
| `supplementary/` | Supplementary tables and checklists submitted with the manuscript |

## Reproduce

```bash
Rscript analysis/meta_analysis.R
python3 analysis/meta_analysis.py
```

Expected pooled random-effects prevalence: 36.4% [95% CI 21.0 to 53.5], I2 = 87.1%, Q = 38.87 (df = 5), p < 0.001.

## Methods summary

| Item | Specification |
|---|---|
| Databases | PubMed/MEDLINE, Embase-equivalent supplementary search, Google Scholar, reference lists |
| Eligibility | Human studies reporting iron overload and quantitative hormonal outcome data in males |
| Risk of bias | Adapted NIH Quality Assessment Tool; ROBINS-I principles |
| Certainty | GRADE |
| Effect model | DerSimonian-Laird random effects |
| Transformation | Freeman-Tukey double arcsine; Miller (1978) back-transform with harmonic mean n |
| Software | R >= 4.3 with `meta`, `metafor`; Python >= 3.10 with `numpy`, `scipy`, `pandas` |

## Citation

Borges JYV. Iron overload and male hypogonadism: a systematic review and meta-analysis. 2026. PROSPERO CRD420261344877.

## License

Code: MIT. Data extraction tables: CC BY 4.0.
