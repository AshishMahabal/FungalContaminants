"""Pure-logic contamination checker.

No Streamlit imports here so the module is testable in isolation.

Curated DB schema: columns ``Phyla``, ``Species``, then one column per property.
Property values are tri-valued: 0 (absent), 1 (35% identity evidence),
2 (75% identity evidence). Score per matched curated row is the weighted
sum ``Σ weight[prop] * value[prop]`` over property columns. Group entries
(``Genus sp.``) score the mean across all curated members of the genus.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
import pandas as pd

from .properties import format_score

GROUP_SUFFIXES = (" sp.", " spp.", " sp", " spp")
META_COLUMNS = ("Phyla", "Species")
GROUP_AGGREGATORS = ("mean", "max", "sum")
_AGG_LABEL = {"mean": "avg", "max": "max", "sum": "sum"}


@dataclass
class FilterResult:
    """Output of :func:`filter_fungi`."""

    filtered: pd.DataFrame  # one row per matched input species above thresholds
    group_stats: pd.DataFrame  # Group, Genus, Number of matches, Matched curated species
    reverse_table: pd.DataFrame  # one row per contributing property
    unmatched: pd.DataFrame  # input rows whose species didn't match the curated DB
    n_input: int
    n_matched_to_curated: int
    n_above_thresholds: int

    @property
    def n_unmatched(self) -> int:
        return len(self.unmatched)


def normalize_species_name(name: object) -> tuple[str | None, bool]:
    """Return ``(normalized_name, is_group)``.

    A "group" entry is something like ``Candida sp.`` — we strip the suffix
    and treat the first token as a genus query.
    """
    if name is None or (isinstance(name, float) and pd.isna(name)):
        return None, False
    text = str(name).strip()
    if not text:
        return None, False
    lower = text.lower()
    for suffix in GROUP_SUFFIXES:
        if lower.endswith(suffix):
            genus = text.split()[0]
            return genus, True
    return text, False


def get_property_columns(curated_df: pd.DataFrame) -> list[str]:
    """Property columns are everything except the metadata columns."""
    return [c for c in curated_df.columns if c not in META_COLUMNS]


def reconcile_input_names(
    input_df: pd.DataFrame, curated_df: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Best-effort rewrite of input species names to match curated names.

    Strategy: for each input species, take the first two tokens (genus + species)
    and find a curated species whose first two tokens match case-insensitively.
    This is stricter than the original "substring in" check, which would match
    ``Candida`` inside any longer name.

    Returns the (possibly renamed) input dataframe and a mapping
    ``{curated_name: original_input_name}`` for display.
    """
    df = input_df.copy()
    species_col = df.columns[0]
    if species_col != "#Datasets":
        df = df.rename(columns={species_col: "#Datasets"})
        species_col = "#Datasets"

    # Drop rows where the species cell is numeric (sample CSVs sometimes have stray rows).
    is_numeric = pd.to_numeric(df[species_col], errors="coerce").notna()
    df = df[~is_numeric].copy()

    curated_names = curated_df["Species"].tolist()
    curated_first2 = {
        " ".join(name.split()[:2]).lower(): name for name in curated_names
    }

    rename_map: dict[str, str] = {}
    new_names: list[str] = []
    for raw in df[species_col]:
        text = str(raw).strip()
        first2 = " ".join(text.split()[:2]).lower()
        match = curated_first2.get(first2)
        if match and match != text:
            rename_map[match] = text
            new_names.append(match)
        else:
            new_names.append(text)
    df[species_col] = new_names

    if len(df.columns) == 1:
        df["sample_loc1"] = 100
    return df, rename_map


def _expand_groups(
    input_df: pd.DataFrame, curated_species: list[str]
) -> tuple[list[int], dict[int, list[str]], list[dict]]:
    """For each input row, find curated species it represents.

    Returns matched row indices, ``{idx: [curated species, ...]}``, and a
    list of group-stat dicts (one per group entry).
    """
    species_col = input_df.columns[0]
    matched_indices: list[int] = []
    group_mapping: dict[int, list[str]] = {}
    group_stats: list[dict] = []
    curated_set = set(curated_species)

    for idx, raw in input_df[species_col].items():
        norm, is_group = normalize_species_name(raw)
        if norm is None:
            continue
        if is_group:
            prefix = norm.lower() + " "
            matches = [s for s in curated_species if s.lower().startswith(prefix)]
            if matches:
                matched_indices.append(idx)
                group_mapping[idx] = matches
                group_stats.append(
                    {
                        "Group": raw,
                        "Genus": norm,
                        "Number of matches": len(matches),
                        "Matched curated species": matches,
                    }
                )
        else:
            if norm in curated_set:
                matched_indices.append(idx)
                group_mapping[idx] = [norm]
    return matched_indices, group_mapping, group_stats


def _score_row(
    matched_curated: list[str],
    curated_df: pd.DataFrame,
    properties: list[str],
    weights: np.ndarray,
    aggregator: str = "mean",
) -> tuple[float, float, list[str], dict[str, int]]:
    """Compute (A-score, S-score, contributing property names).

    - **A-score** = ``Σ weight × value`` per matched curated species,
      then aggregated across the genus by ``aggregator``. Sensitive to
      evidence strength (value 1 vs 2).
    - **S-score** = count of properties with ``weight × value > 0`` per
      matched species, then aggregated. Counts breadth, not depth.
    - **Contributing properties** = union across matched species of
      properties whose ``weight × value > 0`` for at least one member.

    ``aggregator`` is ``mean | max | sum`` and only matters for group
    entries (``Genus sp.``) that match more than one curated species.
    """
    if aggregator not in GROUP_AGGREGATORS:
        raise ValueError(
            f"aggregator must be one of {GROUP_AGGREGATORS}, got {aggregator!r}"
        )
    rows = curated_df[curated_df["Species"].isin(matched_curated)]
    if rows.empty:
        return 0.0, 0.0, [], {}
    values = rows[properties].to_numpy(dtype=float)  # (n_matches, n_props)
    weighted = values * weights  # broadcast over rows
    a_scores = weighted.sum(axis=1)
    s_scores = (weighted > 0).sum(axis=1).astype(float)
    contributing_mask = weighted.sum(axis=0) > 0
    contributing = [p for p, hit in zip(properties, contributing_mask) if hit]
    # Per-property fingerprint: max raw value across matched species,
    # zeroed where the user weight is 0 (so the visual matches the score).
    raw_max = values.max(axis=0)
    masked = np.where(weights > 0, raw_max, 0)
    fingerprint = {p: int(v) for p, v in zip(properties, masked)}
    if aggregator == "mean":
        a_agg, s_agg = float(a_scores.mean()), float(s_scores.mean())
    elif aggregator == "max":
        a_agg, s_agg = float(a_scores.max()), float(s_scores.max())
    else:  # sum
        a_agg, s_agg = float(a_scores.sum()), float(s_scores.sum())
    return a_agg, s_agg, contributing, fingerprint


def filter_fungi(
    input_df: pd.DataFrame,
    curated_df: pd.DataFrame,
    score_weights: dict[str, int],
    score_threshold: float,
    reads_threshold: int,
    group_aggregator: str = "mean",
) -> FilterResult:
    """Run the full pipeline and return a :class:`FilterResult`.

    ``group_aggregator`` controls how scores are combined across the curated
    species expanded from a ``Genus sp.`` input row: ``mean`` (default),
    ``max``, or ``sum``. It has no effect on exact-match input rows.
    """
    if group_aggregator not in GROUP_AGGREGATORS:
        raise ValueError(
            f"group_aggregator must be one of {GROUP_AGGREGATORS}, "
            f"got {group_aggregator!r}"
        )
    n_input = len(input_df)
    species_col = input_df.columns[0]
    location_cols = list(input_df.columns[1:])
    curated_species = curated_df["Species"].tolist()

    matched_indices, group_mapping, group_stats = _expand_groups(
        input_df, curated_species
    )
    matching = input_df.loc[matched_indices].copy()
    unmatched_indices = [i for i in input_df.index if i not in set(matched_indices)]
    unmatched = input_df.loc[unmatched_indices].copy().reset_index(drop=True)

    if matching.empty:
        empty = pd.DataFrame()
        return FilterResult(
            filtered=empty,
            group_stats=pd.DataFrame(group_stats).drop_duplicates(subset=["Group"])
            if group_stats
            else empty,
            reverse_table=empty,
            unmatched=unmatched,
            n_input=n_input,
            n_matched_to_curated=0,
            n_above_thresholds=0,
        )

    reads = matching[location_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    above_reads = reads >= reads_threshold
    matching["Num loc"] = above_reads.sum(axis=1)
    matching["Locations"] = [
        {loc: int(reads.loc[i, loc]) for loc in location_cols if above_reads.loc[i, loc]}
        for i in matching.index
    ]
    matching = matching[matching["Num loc"] > 0].copy()
    n_matched = len(matching)

    if matching.empty:
        empty = pd.DataFrame()
        return FilterResult(
            filtered=empty,
            group_stats=pd.DataFrame(group_stats).drop_duplicates(subset=["Group"])
            if group_stats
            else empty,
            reverse_table=empty,
            unmatched=unmatched,
            n_input=n_input,
            n_matched_to_curated=n_matched,
            n_above_thresholds=0,
        )

    properties = [p for p in score_weights if p in curated_df.columns]
    weights = np.array([score_weights[p] for p in properties], dtype=float)

    a_scores, s_scores, contributing, fingerprints = [], [], [], []
    for idx in matching.index:
        a, s, props, fp = _score_row(
            group_mapping[idx], curated_df, properties, weights, group_aggregator
        )
        a_scores.append(a)
        s_scores.append(s)
        contributing.append(props)
        fingerprints.append(fp)
    matching["Weight Score"] = a_scores
    matching["S Score"] = s_scores
    matching["Contributing Properties"] = contributing
    matching["Property Values"] = fingerprints

    matching = matching[matching["Weight Score"] >= score_threshold].copy()
    n_above = len(matching)
    if matching.empty:
        empty = pd.DataFrame()
        return FilterResult(
            filtered=empty,
            group_stats=pd.DataFrame(group_stats).drop_duplicates(subset=["Group"])
            if group_stats
            else empty,
            reverse_table=empty,
            unmatched=unmatched,
            n_input=n_input,
            n_matched_to_curated=n_matched,
            n_above_thresholds=0,
        )

    suffix = _AGG_LABEL[group_aggregator]

    def _label(value: float, is_group: bool) -> str:
        text = format_score(value)
        return f"{text} ({suffix})" if is_group else text

    is_group_per_row = matching[species_col].apply(
        lambda n: normalize_species_name(n)[1]
    )
    matching["A-score"] = [
        _label(v, ig) for v, ig in zip(matching["Weight Score"], is_group_per_row)
    ]
    matching["S-score"] = [
        _label(v, ig) for v, ig in zip(matching["S Score"], is_group_per_row)
    ]
    filtered = matching[
        [species_col, "A-score", "S-score", "Property Values",
         "Contributing Properties", "Num loc", "Locations"]
    ].reset_index(drop=True)

    reverse_rows = []
    for prop in sorted({p for props in contributing for p in props}):
        species_for_prop = filtered[
            filtered["Contributing Properties"].apply(lambda x: prop in x)
        ][species_col].tolist()
        reverse_rows.append(
            {
                "Property": prop,
                "Number": len(species_for_prop),
                "Matching Species": species_for_prop,
            }
        )
    reverse_table = pd.DataFrame(reverse_rows)

    group_stats_df = (
        pd.DataFrame(group_stats).drop_duplicates(subset=["Group"])
        if group_stats
        else pd.DataFrame()
    )

    return FilterResult(
        filtered=filtered,
        group_stats=group_stats_df,
        reverse_table=reverse_table,
        unmatched=unmatched,
        n_input=n_input,
        n_matched_to_curated=n_matched,
        n_above_thresholds=n_above,
    )


def unique_contributing_properties(filtered: pd.DataFrame) -> list[str]:
    """Flatten the ``Contributing Properties`` column into a sorted unique list."""
    if filtered.empty or "Contributing Properties" not in filtered.columns:
        return []
    seen: set[str] = set()
    for props in filtered["Contributing Properties"]:
        if isinstance(props, Iterable):
            seen.update(props)
    return sorted(seen)
