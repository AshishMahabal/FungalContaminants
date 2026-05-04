"""Pure data-prep for visualizations. No plotting library imports here.

Each function turns a :class:`FilterResult` (or its parts) into a flat
dataframe ready to hand to a chart constructor.
"""

from __future__ import annotations

import pandas as pd

from .checker import normalize_species_name
from .properties import fingerprint_string, short_code


def _phylum_for(
    species: str, phylum_lookup: dict[str, str], curated_df: pd.DataFrame
) -> str:
    """Best-effort phylum assignment, including for ``Genus sp.`` entries."""
    if species in phylum_lookup:
        return phylum_lookup[species]
    norm, is_group = normalize_species_name(species)
    if is_group and norm:
        prefix = norm.lower() + " "
        mask = curated_df["Species"].str.lower().str.startswith(prefix)
        phyla = curated_df.loc[mask, "Phyla"].unique().tolist()
        if len(phyla) == 1:
            return phyla[0]
        if len(phyla) > 1:
            return "Multiple"
    return "Unknown"


def _row_total_reads(row: pd.Series) -> int:
    """Sum the values in a filtered row's ``Locations`` dict."""
    locations = row.get("Locations") if hasattr(row, "get") else None
    return sum(locations.values()) if isinstance(locations, dict) else 0


def build_sunburst_records(
    filtered: pd.DataFrame, curated_df: pd.DataFrame
) -> pd.DataFrame:
    """Flat dataframe for a Phyla → Species sunburst.

    Sized by total reads above threshold. Reads are NOT split; each species
    contributes its full reads. Group entries (``Genus sp.``) inherit a
    phylum if all matched curated species share one, else ``Multiple``.

    Returns columns ``Phyla``, ``Species``, ``Reads``.
    """
    columns = ["Phyla", "Species", "Reads"]
    if filtered is None or filtered.empty:
        return pd.DataFrame(columns=columns)

    species_col = filtered.columns[0]
    phylum_lookup = dict(zip(curated_df["Species"], curated_df["Phyla"]))

    records: list[dict] = []
    for _, row in filtered.iterrows():
        total_reads = _row_total_reads(row)
        if total_reads <= 0:
            continue
        species = row[species_col]
        records.append(
            {
                "Phyla": _phylum_for(species, phylum_lookup, curated_df),
                "Species": species,
                "Reads": total_reads,
            }
        )
    return pd.DataFrame(records, columns=columns)


def build_property_bar_records(filtered: pd.DataFrame) -> pd.DataFrame:
    """Flat dataframe for a 'reads per property' bar chart.

    For each contributing property, sums the reads of every flagged species
    that exhibits the property and counts the species. Reads are NOT split
    across properties — a species exhibiting three properties contributes
    its full reads to each of those three property bars (the chart compares
    properties to each other; bars are not meant to sum to a grand total).

    Returns columns ``Property``, ``Reads``, ``Species count``,
    sorted by ``Reads`` descending.
    """
    columns = ["Property", "Reads", "Species count"]
    if filtered is None or filtered.empty:
        return pd.DataFrame(columns=columns)

    by_prop: dict[str, dict[str, int]] = {}
    for _, row in filtered.iterrows():
        total_reads = _row_total_reads(row)
        if total_reads <= 0:
            continue
        props = row.get("Contributing Properties") if hasattr(row, "get") else None
        if not isinstance(props, list):
            continue
        for prop in props:
            entry = by_prop.setdefault(prop, {"Reads": 0, "Species count": 0})
            entry["Reads"] += total_reads
            entry["Species count"] += 1

    df = pd.DataFrame(
        [
            {"Property": p, "Reads": v["Reads"], "Species count": v["Species count"]}
            for p, v in by_prop.items()
        ],
        columns=columns,
    )
    if not df.empty:
        df = df.sort_values("Reads", ascending=False).reset_index(drop=True)
    return df


def build_heatmap_matrix(
    filtered: pd.DataFrame, input_df: pd.DataFrame
) -> pd.DataFrame:
    """Species × location reads matrix for flagged species, from raw input.

    Index = species name (from the filtered table), columns = location names
    in input order, values = raw read counts (including below-threshold ones,
    so the user sees the full per-location distribution).
    """
    if filtered is None or filtered.empty or input_df is None or input_df.empty:
        return pd.DataFrame()
    species_col = filtered.columns[0]
    if species_col not in input_df.columns:
        return pd.DataFrame()
    flagged = filtered[species_col].tolist()
    location_cols = [c for c in input_df.columns if c != species_col]
    sub = input_df[input_df[species_col].isin(flagged)].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = sub.set_index(species_col)[location_cols]
    sub = sub.apply(pd.to_numeric, errors="coerce").fillna(0)
    # Preserve filtered-table species ordering when possible.
    sub = sub.reindex([s for s in flagged if s in sub.index])
    return sub


def build_upset_indicators(filtered: pd.DataFrame) -> pd.DataFrame:
    """Boolean species × property indicator matrix for upsetplot.

    Each row is a flagged species; each column is a contributing property
    that appears at least once across the filtered set; cell is True iff
    that species exhibits the property.
    """
    if filtered is None or filtered.empty:
        return pd.DataFrame()
    species_col = filtered.columns[0]
    species_props: list[tuple[str, list[str]]] = []
    all_props: list[str] = []
    seen_props: set[str] = set()
    for _, row in filtered.iterrows():
        species = row[species_col]
        props = row.get("Contributing Properties") if hasattr(row, "get") else None
        if not isinstance(props, list):
            continue
        species_props.append((species, props))
        for p in props:
            if p not in seen_props:
                seen_props.add(p)
                all_props.append(p)
    if not species_props:
        return pd.DataFrame()
    df = pd.DataFrame(
        False,
        index=[s for s, _ in species_props],
        columns=all_props,
        dtype=bool,
    )
    df.index.name = "Species"
    for species, props in species_props:
        for p in props:
            df.loc[species, p] = True
    return df


def format_filtered_for_display(filtered: pd.DataFrame) -> pd.DataFrame:
    """Compact list/dict/value columns into strings for table rendering.

    - ``Property Values`` (dict {prop: 0|1|2}) → 6 Unicode squares at fixed
      positions (red=2, yellow=1, white=0). Drops the column ``Contributing
      Properties`` since the squares already encode it (with strength).
    - ``Locations`` (dict) → ``"loc1: 200, loc2: 1240"``, sorted by key.
    """
    if filtered is None or filtered.empty:
        return pd.DataFrame() if filtered is None else filtered.copy()
    df = filtered.copy()
    if "Property Values" in df.columns:
        df["Property Values"] = df["Property Values"].apply(
            lambda v: fingerprint_string(v) if isinstance(v, dict) else ""
        )
    if "Contributing Properties" in df.columns:
        df = df.drop(columns=["Contributing Properties"])
    if "Locations" in df.columns:
        df["Locations"] = df["Locations"].apply(
            lambda d: ", ".join(f"{k}: {v}" for k, v in sorted(d.items()))
            if isinstance(d, dict)
            else ""
        )
    return df
