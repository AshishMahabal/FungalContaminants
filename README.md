# Fungal Contamination Checker

A Streamlit dashboard that flags likely fungal contaminants in a sequencing
or microbiome dataset by scoring observed species against a curated database
of risk-relevant properties.

Built on top of the upstream [FungalContaminants](https://github.com/AshishMahabal/FungalContaminants)
tool used at Caltech / TRISH / NASA-CSIF.

## What it does

You provide a CSV of species detected across one or more locations
(read counts per location). The app:

1. **Matches** each input species against a curated database of ~1,500
   fungi annotated with six contamination-relevant properties:
   antimicrobial resistance, biofilm formation, human pathogenicity,
   thermophily, radiation resistance, spore formation.
2. **Scores** each match using tri-valued evidence per property:
   - `0` — property not observed
   - `1` — observed at ≥35% sequence identity (weak evidence)
   - `2` — observed at ≥75% sequence identity (strong evidence)
3. **Filters** by two thresholds: minimum weighted score, and minimum
   read count in at least one location.
4. **Visualises** what's flagged: a sortable table, an interactive
   sunburst (Phyla → Species), and per-property breakdowns.

The score is a weighted sum: `score = Σ weight[prop] × evidence[prop]`.
Weights default to 1 per property, adjustable in the sidebar (0 ignores a
property; 2 double-weights it).

## Inputs

A CSV with the species name in the first column and read counts in the
remaining columns. Group entries (`Genus sp.`) are expanded against the
curated database; their score is aggregated across the genus by mean
(default), max, or sum.

```
#Datasets,loc1,loc2,loc3
Candida albicans,200,1240,0
Aspergillus sp.,300,4240,0
```

## Running locally

```bash
cd app
pip install -r requirements.txt
streamlit run main.py
```

Tests:

```bash
cd app
python -m pytest tests/
```

## Project layout

```
app/
├── main.py              # Streamlit UI (sidebar config + main results pane)
├── core/
│   ├── checker.py       # Pure scoring/filtering logic
│   └── viz_prep.py      # Pure dataframe prep for charts
├── data/
│   ├── curated_fungi_both.csv   # Merged tri-valued curated DB
│   ├── score_weights.txt
│   └── samples/         # Bundled sample inputs
└── tests/               # pytest suite
scripts/                 # One-off data-prep helpers (merge, query, score)
FungalContaminants/      # Upstream clone, read-only reference
```

## Credits

- Concept: Ashish Mahabal (Caltech), Nitin K. Singh
- Domain expertise: Swati Bijlani (COH), Nitin K. Singh
- Original Streamlit app: Vannsh Jani (Caltech VURP '25)
- Funded in part by TRISH (NASA-funded) via the Caltech Space-Health
  Innovation Fund (CSIF).
