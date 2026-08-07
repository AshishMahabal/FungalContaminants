"""Streamlit UI for the fungal contamination checker.

Sidebar holds all configuration; the main pane shows a dashboard.
Logic lives in :mod:`core.checker` and :mod:`core.viz_prep` for testability.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from core.checker import (
    filter_fungi,
    get_property_columns,
    load_synonyms,
    reconcile_input_names,
)
from core.properties import code_legend, short_code
from core.viz_prep import (
    build_heatmap_matrix,
    build_property_bar_records,
    build_upset_indicators,
    format_filtered_for_display,
)
from core.viz_upset import make_upset_figure

APP_DIR = Path(__file__).parent
DATA_DIR = APP_DIR / "data"
SAMPLES_DIR = DATA_DIR / "samples"
CURATED_PATH = DATA_DIR / "curated_fungi_both.csv"
WEIGHTS_PATH = DATA_DIR / "score_weights.txt"
SYNONYMS_PATH = DATA_DIR / "synonyms.csv"
CREDITS_PATH = APP_DIR / "CREDITS.md"

SAMPLES: dict[str, Path] = {
    "Lab cell culture": SAMPLES_DIR / "lab_culture.csv",
    "ISS surface microbiome": SAMPLES_DIR / "iss_environment.csv",
    "Hospital ICU air sampling": SAMPLES_DIR / "hospital_icu.csv",
}
UPLOAD_LABEL = "Upload my own…"

st.set_page_config(
    page_title="Fungal Contamination Checker",
    page_icon="🧫",
    layout="wide",
)


@st.cache_data
def load_curated() -> pd.DataFrame:
    return pd.read_csv(CURATED_PATH)


@st.cache_data
def load_default_weights() -> dict[str, int]:
    with open(WEIGHTS_PATH) as f:
        return json.load(f)


@st.cache_data
def load_synonym_map() -> dict[str, str]:
    return load_synonyms(SYNONYMS_PATH)


curated_df = load_curated()
property_cols = get_property_columns(curated_df)
default_weights = load_default_weights()
synonym_map = load_synonym_map()
seed_weights = {p: int(default_weights.get(p, 1)) for p in property_cols}

if "score_weights" not in st.session_state:
    st.session_state["score_weights"] = seed_weights.copy()


def _read_input_for(choice: str, uploaded) -> pd.DataFrame | None:
    if choice == UPLOAD_LABEL:
        if uploaded is None:
            st.error("Please upload a CSV file or pick a sample dataset.")
            return None
        try:
            return pd.read_csv(uploaded)
        except Exception as e:
            st.error(f"Could not read uploaded CSV: {e}")
            return None
    return pd.read_csv(SAMPLES[choice])


def run_analysis(
    choice: str,
    uploaded,
    weights: dict[str, int],
    score_threshold: float,
    reads_threshold: int,
    group_aggregator: str,
) -> None:
    input_df = _read_input_for(choice, uploaded)
    if input_df is None:
        return
    input_df, rename_map, notices = reconcile_input_names(
        input_df, curated_df, synonym_map
    )
    result = filter_fungi(
        input_df,
        curated_df,
        weights,
        score_threshold,
        reads_threshold,
        group_aggregator=group_aggregator,
    )
    st.session_state["result"] = result
    st.session_state["input_df"] = input_df
    st.session_state["rename_map"] = rename_map
    st.session_state["row_notices"] = notices
    st.session_state["dataset_label"] = choice


# ====================================================================
# Sidebar: all configuration
# ====================================================================
with st.sidebar:
    st.header("Configuration")

    sample_choice = st.selectbox(
        "Dataset",
        list(SAMPLES.keys()) + [UPLOAD_LABEL],
        index=0,
    )
    uploaded = (
        st.file_uploader("Upload CSV file", type="csv")
        if sample_choice == UPLOAD_LABEL
        else None
    )

    with st.expander("ℹ️  CSV format", expanded=False):
        st.markdown(
            """
- First column: species name (header is treated as a label).
- Remaining columns: read counts per location (any column names).
- Genus-level entries (`Candida sp.`, `Aspergillus spp.`) are
  expanded against the curated database.

Example:

```
#Datasets,loc1,loc2,loc3
Candida albicans,200,1240,0
Aspergillus sp.,300,4240,0
```

Upload limit: **10 MB**.
            """
        )

    st.caption(
        f"Curated DB: `{CURATED_PATH.name}` — {len(curated_df)} species, "
        f"{len(property_cols)} properties (tri-valued 0/1/2)."
    )

    st.divider()

    with st.expander("⚖️  Weights", expanded=False):
        st.caption("0 = ignore property; 2 = double-weight it.")
        # Apply any pending reset/upload *before* the sliders are created:
        # Streamlit forbids writing a widget's own key after it is instantiated.
        pending = st.session_state.pop("_pending_weights", None)
        if pending is not None:
            for prop in property_cols:
                if prop in pending:
                    st.session_state[f"w_{prop}"] = int(pending[prop])
        new_weights: dict[str, int] = {}
        for prop in property_cols:
            code = short_code(prop)
            label = f"{code}" if code != prop else prop
            new_weights[prop] = st.slider(
                label,
                0,
                2,
                int(st.session_state["score_weights"].get(prop, 1)),
                key=f"w_{prop}",
                help=prop,
            )
        st.session_state["score_weights"] = new_weights
        legend = code_legend(property_cols)
        if legend:
            st.caption(f"Codes: {legend}")

        if st.button("Restore defaults"):
            st.session_state["_pending_weights"] = dict(seed_weights)
            st.rerun()

        custom = st.file_uploader("Upload weights JSON", type="json")
        if custom is not None and st.session_state.get("_weights_file_id") != custom.file_id:
            st.session_state["_weights_file_id"] = custom.file_id
            loaded = json.load(custom)
            st.session_state["_pending_weights"] = {
                prop: int(loaded[prop]) for prop in property_cols if prop in loaded
            }
            st.rerun()

    st.markdown("**Thresholds**")
    group_aggregator = st.radio(
        "Group scoring (Genus sp.)",
        ["mean", "max", "sum"],
        index=0,
        horizontal=True,
        help=(
            "How to combine scores across the curated species expanded from a "
            "Genus sp. input. mean = average member; max = strongest single "
            "member; sum = total across members (grows with genus size)."
        ),
    )
    score_threshold = st.number_input(
        "A-score threshold",
        min_value=0.0,
        value=4.0,
        step=0.5,
        help=(
            "Minimum A-score (weighted evidence sum) to flag a species. "
            "S-score (count of contributing properties) is shown for "
            "reference but not currently used as a filter."
        ),
    )
    reads_threshold = st.slider(
        "Reads threshold",
        min_value=1,
        max_value=10000,
        value=10,
        help="Per-location minimum read count for a sample to count toward Num loc.",
    )

    st.divider()
    run_clicked = st.button("🚀 Run analysis", type="primary", width="stretch")

    with st.expander("ℹ️  About", expanded=False):
        if CREDITS_PATH.exists():
            st.markdown(CREDITS_PATH.read_text(encoding="utf-8"))
        else:
            st.caption("Credits file missing.")

# ====================================================================
# Auto-run on first load using whatever the sidebar currently shows
# ====================================================================
if "first_load_done" not in st.session_state:
    st.session_state["first_load_done"] = True
    if sample_choice != UPLOAD_LABEL:
        run_analysis(
            sample_choice,
            uploaded,
            st.session_state["score_weights"],
            score_threshold,
            reads_threshold,
            group_aggregator,
        )

if run_clicked:
    run_analysis(
        sample_choice,
        uploaded,
        st.session_state["score_weights"],
        score_threshold,
        reads_threshold,
        group_aggregator,
    )

# ====================================================================
# Main pane
# ====================================================================
st.title("🧫 Fungal Contamination Checker")
st.caption(
    "Score = Σ weight × evidence over a curated property set "
    "(tri-valued: 0 absent, 1 ≥35% identity, 2 ≥75% identity)."
)

with st.expander("ℹ️  How to read this dashboard", expanded=False):
    st.markdown(
        """
- **Pick a dataset** (or upload your own CSV) in the sidebar.
- The app matches each species against ~1,500 curated fungi annotated for
  six contamination-relevant properties: *antimicrobial resistance, biofilm
  formation, human pathogenicity, thermophily, radiation resistance,
  spore formation*.
- Each match gets two scores:
  - **A-score** (depth) = Σ weight × evidence — sensitive to evidence
    strength. Evidence per property is `0`, `1` (≥35% identity), or `2`
    (≥75% identity).
  - **S-score** (breadth) = count of contributing properties.
- A species is **flagged** if its A-score ≥ *A-score threshold* AND at
  least one location has reads ≥ *Reads threshold*. S-score is shown for
  reference.
- Genus-level inputs (`Candida sp.`) expand to all curated *Candida*
  species; their score is aggregated by mean / max / sum (sidebar choice).
- The **sunburst** shows where the contamination concentrates
  (Phyla → Species, sized by total reads). Click to drill in.
- Tweak weights or thresholds in the sidebar and click **Run analysis**.
        """
    )

if "result" not in st.session_state:
    st.info("Set parameters in the sidebar and click **Run analysis**.")
else:
    result = st.session_state["result"]
    label = st.session_state.get("dataset_label", "")
    if label:
        st.markdown(f"**Dataset:** {label}")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Input records", result.n_input)
    c2.metric("Matched to curated", result.n_matched_to_curated)
    c3.metric("Above threshold", result.n_above_thresholds)
    c4.metric("Unmatched", result.n_unmatched)

    notices = st.session_state.get("row_notices", [])
    if notices:
        lines = []
        for n in notices:
            if n["kind"] == "synonym":
                names = ", ".join(f"*{x}*" for x in n["inputs"])
                lines.append(
                    f"- **{n['organism']}** — {names} are the same organism; "
                    "their reads were **summed**."
                )
            else:  # duplicate
                detail = (
                    "read counts differed — kept the **larger** per location"
                    if n["reads_differed"]
                    else "identical rows, kept one"
                )
                lines.append(
                    f"- **{n['name']}** — {n['rows']} rows with the same name "
                    f"({detail})."
                )
        st.info(
            "🔗 **Combined input rows** (so nothing is counted twice):\n"
            + "\n".join(lines)
        )
        st.caption(
            "Different names for one organism are summed (see `synonyms.csv`); "
            "repeated identical names keep the larger read count."
        )

    if result.n_unmatched > 0:
        with st.expander(
            f"⚠️  {result.n_unmatched} unmatched species "
            "(not in the curated DB)",
            expanded=False,
        ):
            st.dataframe(result.unmatched, width="stretch", hide_index=True)
            st.download_button(
                "Download unmatched CSV",
                data=result.unmatched.to_csv(index=False),
                file_name="unmatched_species.csv",
                mime="text/csv",
                help=(
                    "Species in your input that don't appear in the curated "
                    "database (no scoring possible). Worth checking for typos, "
                    "alternate naming, or candidates to add to the curated DB."
                ),
            )

    if result.filtered.empty:
        st.warning(
            "No contamination detected at these thresholds. "
            "Try lowering Score threshold or Reads threshold in the sidebar."
        )
    else:
        st.subheader("Filtered contamination")
        species_col = result.filtered.columns[0]
        display_df = format_filtered_for_display(result.filtered)
        # Size to content; cap at 420 so the header stays sticky for long
        # tables. Streamlit row ≈ 35px, header ≈ 38px.
        table_height = min(420, 38 + 35 * len(display_df) + 8)
        st.dataframe(
            display_df,
            width="stretch",
            height=table_height,
            hide_index=True,
            column_config={
                species_col: st.column_config.TextColumn(
                    "Species", width=200
                ),
                "A-score": st.column_config.TextColumn(
                    "A", width=70,
                    help=(
                        "A-score = Σ weight × evidence (depth — sensitive to "
                        "evidence strength: value 1 vs 2). Filtered by "
                        "Score threshold in the sidebar."
                    ),
                ),
                "S-score": st.column_config.TextColumn(
                    "S", width=70,
                    help=(
                        "S-score = count of contributing properties (breadth). "
                        "Not filtered, shown for reference."
                    ),
                ),
                "Property Values": st.column_config.TextColumn(
                    "AMR BF HP RAD SF TH", width=160,
                    help=(
                        "Per-property strength fingerprint at fixed positions. "
                        "🟥 = strong (≥75% identity) · 🟨 = weak (≥35%) · "
                        "⬜ = none. For genus groups: max value across members."
                    ),
                ),
                "Num loc": st.column_config.NumberColumn(
                    "# loc", width=60
                ),
                "Locations": st.column_config.TextColumn(
                    "Reads per location", width="large"
                ),
            },
        )
        st.download_button(
            "Download filtered CSV",
            data=result.filtered.to_csv(index=False),
            file_name="filtered_fungi.csv",
            mime="text/csv",
        )

        with st.expander("📍  Where is the contamination? (heatmap)", expanded=True):
            heat = build_heatmap_matrix(
                result.filtered, st.session_state.get("input_df")
            )
            if heat.empty:
                st.caption("No reads to plot.")
            else:
                use_log = st.checkbox(
                    "Log colour scale",
                    value=True,
                    help="log(reads + 1). Recommended when reads span orders of magnitude.",
                    key="heatmap_log",
                )
                z = np.log1p(heat.values) if use_log else heat.values
                hover_text = [
                    [f"{sp}<br>{loc}<br>reads = {int(heat.values[i, j])}"
                     for j, loc in enumerate(heat.columns)]
                    for i, sp in enumerate(heat.index)
                ]
                fig = px.imshow(
                    z,
                    x=heat.columns,
                    y=heat.index,
                    color_continuous_scale="Tealgrn",
                    aspect="auto",
                    labels=dict(
                        x="Location",
                        y="Species",
                        color="log(reads+1)" if use_log else "reads",
                    ),
                )
                fig.update_traces(
                    customdata=hover_text,
                    hovertemplate="%{customdata}<extra></extra>",
                )
                fig.update_layout(
                    margin=dict(t=20, l=0, r=0, b=0),
                    height=min(700, max(220, 28 * len(heat.index) + 80)),
                )
                st.plotly_chart(fig, width="stretch")
                st.caption(
                    "Cell = raw reads (incl. below the Reads threshold). "
                    "Hover for exact counts."
                )

        with st.expander("🔗  Which properties co-occur? (UpSet)", expanded=False):
            ind = build_upset_indicators(result.filtered)
            if ind.empty or ind.shape[1] < 2:
                st.caption(
                    "Need at least two contributing properties to draw an UpSet plot."
                )
            else:
                fig = make_upset_figure(ind)
                st.plotly_chart(fig, width="stretch")
                st.caption(
                    "Top bars = species count for each property combination "
                    "(sorted by combo size). Left bars = total species per "
                    "property. Hover any bar or dot for member species."
                )

        with st.expander("Reads per property (bar chart)", expanded=False):
            bar_df = build_property_bar_records(result.filtered)
            if bar_df.empty:
                st.caption("No contributing properties.")
            else:
                bar_df = bar_df.assign(Code=bar_df["Property"].map(short_code))
                fig = px.bar(
                    bar_df,
                    x="Reads",
                    y="Code",
                    orientation="h",
                    text="Species count",
                    color="Reads",
                    color_continuous_scale="Tealgrn",
                    hover_data={"Property": True, "Code": False, "Reads": True,
                                "Species count": True},
                )
                fig.update_layout(
                    yaxis={"categoryorder": "total ascending", "title": "Property"},
                    margin=dict(t=20, l=0, r=0, b=0),
                    height=320,
                    coloraxis_showscale=False,
                )
                fig.update_traces(textposition="outside")
                st.plotly_chart(fig, width="stretch")
                st.caption(
                    "Sum of reads across flagged species exhibiting each "
                    "property; bar text = species count."
                )

        if not result.reverse_table.empty:
            with st.expander("Species per property (table)", expanded=False):
                rev = result.reverse_table.copy()
                rev["Code"] = rev["Property"].map(short_code)
                rev["Matching Species"] = rev["Matching Species"].apply(
                    lambda xs: ", ".join(xs) if isinstance(xs, list) else ""
                )
                rev = rev[["Code", "Property", "Number", "Matching Species"]]
                st.dataframe(
                    rev,
                    width="stretch",
                    hide_index=True,
                    column_config={
                        "Code": st.column_config.TextColumn("Code", width=60),
                        "Property": st.column_config.TextColumn(
                            "Property", width=180
                        ),
                        "Number": st.column_config.NumberColumn(
                            "Species", width=70
                        ),
                        "Matching Species": st.column_config.TextColumn(
                            "Matching species", width="large"
                        ),
                    },
                )

        if not result.group_stats.empty:
            with st.expander("Group match statistics", expanded=False):
                st.dataframe(result.group_stats, width="stretch")
