"""Tests for app/core/viz_prep.py."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from core.checker import filter_fungi
from core.viz_prep import (
    build_heatmap_matrix,
    build_property_bar_records,
    build_sunburst_records,
    build_upset_combinations,
    build_upset_indicators,
    build_upset_set_sizes,
)


@pytest.fixture
def curated_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"Phyla": "Ascomycota", "Species": "Candida albicans",
             "Biofilm-formation": 2, "Human-pathogenicity": 2, "Thermophile": 1},
            {"Phyla": "Ascomycota", "Species": "Candida glabrata",
             "Biofilm-formation": 1, "Human-pathogenicity": 2, "Thermophile": 0},
            {"Phyla": "Ascomycota", "Species": "Aspergillus niger",
             "Biofilm-formation": 0, "Human-pathogenicity": 1, "Thermophile": 2},
            {"Phyla": "Basidiomycota", "Species": "Cryptococcus neoformans",
             "Biofilm-formation": 0, "Human-pathogenicity": 0, "Thermophile": 0},
        ]
    )


@pytest.fixture
def weights() -> dict[str, int]:
    return {"Biofilm-formation": 1, "Human-pathogenicity": 1, "Thermophile": 1}


# ---------- build_sunburst_records ----------

def test_sunburst_empty_filtered_returns_empty_with_correct_columns(curated_df):
    df = build_sunburst_records(pd.DataFrame(), curated_df)
    assert list(df.columns) == ["Phyla", "Species", "Reads"]
    assert df.empty


def test_sunburst_single_species_keeps_full_reads(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida albicans"], "loc1": [120]})
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    df = build_sunburst_records(result.filtered, curated_df)
    assert len(df) == 1
    assert df["Phyla"].iloc[0] == "Ascomycota"
    assert df["Species"].iloc[0] == "Candida albicans"
    assert df["Reads"].iloc[0] == 120


def test_sunburst_no_property_column_anymore(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida albicans"], "loc1": [120]})
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    df = build_sunburst_records(result.filtered, curated_df)
    assert "Property" not in df.columns


def test_sunburst_multi_species_each_gets_one_row(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Aspergillus niger"], "loc1": [100, 250]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    df = build_sunburst_records(result.filtered, curated_df)
    assert len(df) == 2
    by_species = dict(zip(df["Species"], df["Reads"]))
    assert by_species["Candida albicans"] == 100
    assert by_species["Aspergillus niger"] == 250


def test_sunburst_genus_group_single_phylum(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida sp."], "loc1": [60]})
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    df = build_sunburst_records(result.filtered, curated_df)
    assert (df["Phyla"] == "Ascomycota").all()
    assert (df["Species"] == "Candida sp.").all()
    assert df["Reads"].iloc[0] == 60


def test_sunburst_genus_group_multiple_phyla_labels_multiple():
    curated = pd.DataFrame(
        [
            {"Phyla": "Ascomycota", "Species": "Foo bar", "P": 1},
            {"Phyla": "Basidiomycota", "Species": "Foo baz", "P": 1},
        ]
    )
    inp = pd.DataFrame({"#Datasets": ["Foo sp."], "loc1": [10]})
    result = filter_fungi(inp, curated, {"P": 1}, score_threshold=0, reads_threshold=1)
    df = build_sunburst_records(result.filtered, curated)
    assert (df["Phyla"] == "Multiple").all()


def test_sunburst_skips_rows_with_zero_reads(curated_df):
    filtered = pd.DataFrame(
        [
            {
                "#Datasets": "Candida albicans",
                "Score": "5.00",
                "Contributing Properties": ["Biofilm-formation"],
                "Num loc": 1,
                "Locations": {},
            }
        ]
    )
    df = build_sunburst_records(filtered, curated_df)
    assert df.empty


def test_sunburst_unknown_species_falls_back_to_unknown_phylum(curated_df):
    filtered = pd.DataFrame(
        [
            {
                "#Datasets": "Mystery genus species",
                "Score": "2.00",
                "Contributing Properties": ["Biofilm-formation"],
                "Num loc": 1,
                "Locations": {"loc1": 50},
            }
        ]
    )
    df = build_sunburst_records(filtered, curated_df)
    assert (df["Phyla"] == "Unknown").all()


# ---------- build_property_bar_records ----------

def test_property_bar_empty_returns_empty_with_correct_columns():
    df = build_property_bar_records(pd.DataFrame())
    assert list(df.columns) == ["Property", "Reads", "Species count"]
    assert df.empty


def test_property_bar_single_species_three_properties(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida albicans"], "loc1": [120]})
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    df = build_property_bar_records(result.filtered)
    # Candida albicans contributes to all 3 properties; full 120 reads each.
    assert len(df) == 3
    assert (df["Reads"] == 120).all()
    assert (df["Species count"] == 1).all()


def test_property_bar_multi_species_sums_reads(curated_df, weights):
    # Candida albicans (100 reads, all 3 props) + Aspergillus niger (200, props HP+TH)
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Aspergillus niger"], "loc1": [100, 200]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    df = build_property_bar_records(result.filtered)
    by_prop = dict(zip(df["Property"], df["Reads"]))
    assert by_prop["Biofilm-formation"] == 100  # only C. albicans
    assert by_prop["Human-pathogenicity"] == 300  # both
    assert by_prop["Thermophile"] == 300  # both


def test_property_bar_sorted_by_reads_desc(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Aspergillus niger"], "loc1": [100, 200]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    df = build_property_bar_records(result.filtered)
    reads = df["Reads"].tolist()
    assert reads == sorted(reads, reverse=True)


def test_property_bar_species_count_per_property(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Aspergillus niger"], "loc1": [100, 200]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    df = build_property_bar_records(result.filtered)
    by_prop = dict(zip(df["Property"], df["Species count"]))
    assert by_prop["Biofilm-formation"] == 1
    assert by_prop["Human-pathogenicity"] == 2
    assert by_prop["Thermophile"] == 2


def test_property_bar_skips_zero_reads_rows():
    filtered = pd.DataFrame(
        [
            {
                "#Datasets": "Candida albicans",
                "Contributing Properties": ["Biofilm-formation"],
                "Num loc": 0,
                "Locations": {},
            }
        ]
    )
    df = build_property_bar_records(filtered)
    assert df.empty


# ---------- build_heatmap_matrix ----------

def test_heatmap_empty_inputs_return_empty(curated_df, weights):
    assert build_heatmap_matrix(pd.DataFrame(), pd.DataFrame()).empty
    assert build_heatmap_matrix(None, pd.DataFrame()).empty


def test_heatmap_uses_raw_input_including_below_threshold(curated_df, weights):
    # Reads threshold of 50 → only loc1 (100) qualifies for the filter,
    # but the heatmap should still show loc2 (5) — it's the raw input.
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans"], "loc1": [100], "loc2": [5]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=50)
    matrix = build_heatmap_matrix(result.filtered, inp)
    assert list(matrix.columns) == ["loc1", "loc2"]
    assert matrix.loc["Candida albicans", "loc1"] == 100
    assert matrix.loc["Candida albicans", "loc2"] == 5


def test_heatmap_only_includes_flagged_species(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Cryptococcus neoformans"],
         "loc1": [100, 100]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    # Cryptococcus has all-zero properties → score 0 → filtered out.
    matrix = build_heatmap_matrix(result.filtered, inp)
    assert matrix.index.tolist() == ["Candida albicans"]


def test_heatmap_preserves_location_column_order(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans"],
         "z_last": [10], "a_first": [20], "m_middle": [30]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=1)
    matrix = build_heatmap_matrix(result.filtered, inp)
    assert list(matrix.columns) == ["z_last", "a_first", "m_middle"]


# ---------- build_upset_indicators ----------

def test_upset_empty_returns_empty():
    assert build_upset_indicators(pd.DataFrame()).empty
    assert build_upset_indicators(None).empty


def test_upset_indicators_marks_contributing_properties(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Aspergillus niger"],
         "loc1": [100, 100]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    ind = build_upset_indicators(result.filtered)
    # Candida albicans → all 3; Aspergillus niger → HP + TH only.
    assert ind.loc["Candida albicans", "Biofilm-formation"]
    assert ind.loc["Candida albicans", "Human-pathogenicity"]
    assert ind.loc["Candida albicans", "Thermophile"]
    assert not ind.loc["Aspergillus niger", "Biofilm-formation"]
    assert ind.loc["Aspergillus niger", "Human-pathogenicity"]
    assert ind.loc["Aspergillus niger", "Thermophile"]


def test_upset_indicators_all_columns_are_bool(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans"], "loc1": [100]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    ind = build_upset_indicators(result.filtered)
    assert all(ind.dtypes == bool)


# ---------- build_upset_combinations ----------

def test_upset_combinations_empty_returns_empty():
    df = build_upset_combinations(pd.DataFrame())
    assert df.empty
    assert list(df.columns) == ["Properties", "Species count", "Members"]


def test_upset_combinations_groups_distinct_property_sets(curated_df, weights):
    # 3 species with distinct property combos in the fixture.
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Candida glabrata", "Aspergillus niger"],
         "loc1": [100, 100, 100]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    ind = build_upset_indicators(result.filtered)
    combos = build_upset_combinations(ind)
    assert len(combos) == 3
    assert all(combos["Species count"] == 1)
    for members in combos["Members"]:
        assert len(members) == 1


def test_upset_combinations_collapses_identical_property_sets():
    # Hand-built indicator matrix: two species with the same combo.
    ind = pd.DataFrame(
        {"A": [True, True, False], "B": [True, True, True], "C": [False, False, True]},
        index=["x", "y", "z"],
    )
    combos = build_upset_combinations(ind)
    # Two combos: {A,B} (x,y) and {B,C} (z).
    assert len(combos) == 2
    top = combos.iloc[0]
    assert top["Species count"] == 2
    assert set(top["Members"]) == {"x", "y"}
    assert set(top["Properties"]) == {"A", "B"}


def test_upset_combinations_sorted_by_species_count_desc():
    ind = pd.DataFrame(
        {"A": [True, True, True, False], "B": [True, True, True, True],
         "C": [False, False, False, True]},
        index=["s1", "s2", "s3", "lonely"],
    )
    combos = build_upset_combinations(ind)
    assert list(combos["Species count"]) == [3, 1]


def test_upset_combinations_skips_species_with_no_active_props():
    ind = pd.DataFrame(
        {"A": [True, False], "B": [False, False]},
        index=["has_one", "has_none"],
    )
    combos = build_upset_combinations(ind)
    assert len(combos) == 1
    assert combos.iloc[0]["Members"] == ["has_one"]


# ---------- build_upset_set_sizes ----------

def test_upset_set_sizes_empty_returns_empty_series():
    s = build_upset_set_sizes(pd.DataFrame())
    assert s.empty


def test_upset_set_sizes_counts_active_per_property_sorted_desc():
    ind = pd.DataFrame(
        {"A": [True, True, False], "B": [True, False, False],
         "C": [True, True, True]},
        index=["x", "y", "z"],
    )
    s = build_upset_set_sizes(ind)
    assert s["A"] == 2
    assert s["B"] == 1
    assert s["C"] == 3
    assert list(s.index) == ["C", "A", "B"]
