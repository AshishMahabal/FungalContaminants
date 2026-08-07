"""Tests for app/core/checker.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from core.checker import (
    FilterResult,
    filter_fungi,
    get_property_columns,
    load_synonyms,
    normalize_species_name,
    reconcile_input_names,
    unique_contributing_properties,
)


@pytest.fixture
def curated_df() -> pd.DataFrame:
    """Tiny curated table covering both binary and tri-valued cases."""
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


# ---------- normalize_species_name ----------

@pytest.mark.parametrize(
    "value,expected",
    [
        ("Candida albicans", ("Candida albicans", False)),
        ("Candida sp.", ("Candida", True)),
        ("Candida sp", ("Candida", True)),
        ("Candida spp.", ("Candida", True)),
        ("Candida spp", ("Candida", True)),
        ("  Candida albicans  ", ("Candida albicans", False)),
        ("", (None, False)),
        (None, (None, False)),
        (np.nan, (None, False)),
    ],
)
def test_normalize_species_name(value, expected):
    assert normalize_species_name(value) == expected


# ---------- get_property_columns ----------

def test_get_property_columns_excludes_metadata(curated_df):
    assert get_property_columns(curated_df) == [
        "Biofilm-formation",
        "Human-pathogenicity",
        "Thermophile",
    ]


# ---------- reconcile_input_names ----------

def test_reconcile_uses_first_two_tokens(curated_df):
    inp = pd.DataFrame(
        {"Sample": ["Candida albicans strain ABC", "Aspergillus niger ATCC 16404"],
         "loc1": [100, 200]}
    )
    out, mapping, merged = reconcile_input_names(inp, curated_df)
    assert out["#Datasets"].tolist() == ["Candida albicans", "Aspergillus niger"]
    assert mapping == {
        "Candida albicans": "Candida albicans strain ABC",
        "Aspergillus niger": "Aspergillus niger ATCC 16404",
    }
    assert merged == []


def test_reconcile_drops_numeric_rows(curated_df):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "42", "Aspergillus niger"],
         "loc1": [10, 20, 30]}
    )
    out, _, _ = reconcile_input_names(inp, curated_df)
    assert "42" not in out["#Datasets"].tolist()
    assert len(out) == 2


def test_reconcile_does_not_match_unrelated_substrings(curated_df):
    inp = pd.DataFrame({"Sample": ["Pseudo candida foo"], "loc1": [100]})
    out, mapping, _ = reconcile_input_names(inp, curated_df)
    assert out["#Datasets"].tolist() == ["Pseudo candida foo"]
    assert mapping == {}


def test_reconcile_adds_default_location_when_only_species_column(curated_df):
    inp = pd.DataFrame({"Sample": ["Candida albicans"]})
    out, _, _ = reconcile_input_names(inp, curated_df)
    assert "sample_loc1" in out.columns
    assert out["sample_loc1"].iloc[0] == 100


# ---------- synonym resolution & merging ----------

def _curated_with_candidozyma() -> pd.DataFrame:
    return pd.DataFrame(
        [{"Phyla": "Ascomycota", "Species": "Candidozyma auris CBS 10913",
          "Human-pathogenicity": 2}]
    )


def test_reconcile_resolves_synonym_alias_to_canonical():
    syn = {"candida auris": "candidozyma auris"}
    inp = pd.DataFrame({"Sample": ["Candida auris"], "loc1": [100]})
    out, mapping, merged = reconcile_input_names(inp, _curated_with_candidozyma(), synonyms=syn)
    assert out["#Datasets"].tolist() == ["Candidozyma auris CBS 10913"]
    assert mapping == {"Candidozyma auris CBS 10913": "Candida auris"}
    assert merged == []


def test_reconcile_merges_two_names_for_same_organism():
    syn = {"candida auris": "candidozyma auris"}
    inp = pd.DataFrame(
        {"Sample": ["Candida auris", "Candidozyma auris"], "loc1": [100, 50]}
    )
    out, mapping, merged = reconcile_input_names(inp, _curated_with_candidozyma(), synonyms=syn)
    assert out["#Datasets"].tolist() == ["Candidozyma auris CBS 10913"]  # one row
    assert out["loc1"].tolist() == [150]                                 # reads summed
    assert len(merged) == 1
    assert merged[0]["kind"] == "synonym"
    assert merged[0]["organism"] == "Candidozyma auris CBS 10913"
    assert set(merged[0]["inputs"]) == {"Candida auris", "Candidozyma auris"}


def test_reconcile_takes_larger_for_duplicate_same_name(curated_df):
    inp = pd.DataFrame(
        {"Sample": ["Candida albicans", "Candida albicans"],
         "loc1": [100, 80], "loc2": [5, 200]}
    )
    out, _, notices = reconcile_input_names(inp, curated_df)
    assert out["#Datasets"].tolist() == ["Candida albicans"]  # one row
    assert out["loc1"].tolist() == [100]                      # larger kept per location
    assert out["loc2"].tolist() == [200]
    assert len(notices) == 1
    assert notices[0]["kind"] == "duplicate"
    assert notices[0]["name"] == "Candida albicans"
    assert notices[0]["rows"] == 2
    assert notices[0]["reads_differed"] is True


def test_reconcile_underscore_variant_is_same_name_duplicate(curated_df):
    inp = pd.DataFrame(
        {"Sample": ["Candida albicans", "Candida_albicans"], "loc1": [30, 70]}
    )
    out, _, notices = reconcile_input_names(inp, curated_df)
    assert out["#Datasets"].tolist() == ["Candida albicans"]
    assert out["loc1"].tolist() == [70]                       # larger kept
    assert len(notices) == 1 and notices[0]["kind"] == "duplicate"


def test_reconcile_without_synonyms_is_backward_compatible(curated_df):
    inp = pd.DataFrame({"Sample": ["Candida auris"], "loc1": [100]})
    out, mapping, merged = reconcile_input_names(inp, curated_df)
    assert out["#Datasets"].tolist() == ["Candida auris"]  # not in tiny DB, unchanged
    assert mapping == {}
    assert merged == []


def test_load_synonyms(tmp_path):
    p = tmp_path / "syn.csv"
    p.write_text(
        "alias,canonical_name,current_accepted_name,relationship\n"
        "Candida_auris,Candidozyma_auris,Candidozyma_auris,synonym\n"
        "Botrytis_cinerea,Botryotinia_fuckeliana,Botrytis_cinerea,synonym\n"
    )
    assert load_synonyms(p) == {
        "candida auris": "candidozyma auris",
        "botrytis cinerea": "botryotinia fuckeliana",
    }


# ---------- filter_fungi ----------

def test_filter_basic_scoring(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Aspergillus niger"],
         "loc1": [100, 100]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    assert isinstance(result, FilterResult)
    assert result.n_above_thresholds == 2
    # Candida albicans: 2 + 2 + 1 = 5 (A); 3 properties contribute (S=3)
    # Aspergillus niger: 0 + 1 + 2 = 3 (A); 2 properties contribute (S=2)
    a_scores = dict(zip(result.filtered["#Datasets"], result.filtered["A-score"]))
    s_scores = dict(zip(result.filtered["#Datasets"], result.filtered["S-score"]))
    assert a_scores["Candida albicans"] == "5"
    assert a_scores["Aspergillus niger"] == "3"
    assert s_scores["Candida albicans"] == "3"
    assert s_scores["Aspergillus niger"] == "2"


def test_filter_score_threshold_excludes_low_scoring(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Cryptococcus neoformans"],
         "loc1": [100, 100]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    # Cryptococcus has all-zero properties → A-score 0 → below threshold.
    assert result.filtered["#Datasets"].tolist() == ["Candida albicans"]


def test_filter_reads_threshold_excludes_low_counts(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans"],
         "loc1": [5], "loc2": [3]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    assert result.filtered.empty


def test_filter_groups_average_across_genus(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida sp."], "loc1": [100]})
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    assert len(result.filtered) == 1
    # mean of Candida albicans (5) and Candida glabrata (3) = 4 (whole)
    assert result.filtered["A-score"].iloc[0] == "4 (avg)"
    # S-score: mean of (3, 2) = 2.5 (fractional)
    assert result.filtered["S-score"].iloc[0] == "2.50 (avg)"
    assert not result.group_stats.empty
    assert result.group_stats["Genus"].iloc[0] == "Candida"
    assert result.group_stats["Number of matches"].iloc[0] == 2


def test_filter_groups_max_aggregator(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida sp."], "loc1": [100]})
    result = filter_fungi(
        inp, curated_df, weights, score_threshold=1, reads_threshold=10,
        group_aggregator="max",
    )
    # max of (5, 3) = 5; S-score max of (3, 2) = 3
    assert result.filtered["A-score"].iloc[0] == "5 (max)"
    assert result.filtered["S-score"].iloc[0] == "3 (max)"


def test_filter_groups_sum_aggregator(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida sp."], "loc1": [100]})
    result = filter_fungi(
        inp, curated_df, weights, score_threshold=1, reads_threshold=10,
        group_aggregator="sum",
    )
    # sum of (5, 3) = 8; S-score sum of (3, 2) = 5
    assert result.filtered["A-score"].iloc[0] == "8 (sum)"
    assert result.filtered["S-score"].iloc[0] == "5 (sum)"


def test_filter_aggregator_does_not_affect_exact_match(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida albicans"], "loc1": [100]})
    for agg in ("mean", "max", "sum"):
        result = filter_fungi(
            inp, curated_df, weights, score_threshold=1, reads_threshold=10,
            group_aggregator=agg,
        )
        assert result.filtered["A-score"].iloc[0] == "5", f"failed for {agg}"
        assert result.filtered["S-score"].iloc[0] == "3", f"failed for {agg}"


def test_filter_invalid_aggregator_raises(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida albicans"], "loc1": [100]})
    with pytest.raises(ValueError):
        filter_fungi(
            inp, curated_df, weights, score_threshold=1, reads_threshold=10,
            group_aggregator="median",
        )


def test_filter_unmatched_input_returns_empty(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Made up species"], "loc1": [100]})
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    assert result.filtered.empty
    assert result.n_matched_to_curated == 0


def test_filter_records_unmatched_species(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Made up species", "Other unknown"],
         "loc1": [100, 200, 300]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    assert result.n_unmatched == 2
    assert set(result.unmatched["#Datasets"]) == {"Made up species", "Other unknown"}
    # Reads preserved
    assert result.unmatched.loc[
        result.unmatched["#Datasets"] == "Made up species", "loc1"
    ].iloc[0] == 200


def test_filter_unmatched_empty_when_all_match(curated_df, weights):
    inp = pd.DataFrame({"#Datasets": ["Candida albicans"], "loc1": [100]})
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    assert result.n_unmatched == 0
    assert result.unmatched.empty


def test_filter_handles_nan_reads(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans"], "loc1": [np.nan], "loc2": [100]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    assert len(result.filtered) == 1
    assert result.filtered["Locations"].iloc[0] == {"loc2": 100}


def test_filter_locations_dict_only_includes_above_threshold(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans"], "loc1": [5], "loc2": [200], "loc3": [50]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    assert result.filtered["Locations"].iloc[0] == {"loc2": 200, "loc3": 50}


def test_filter_reverse_table_lists_species_per_property(curated_df, weights):
    inp = pd.DataFrame(
        {"#Datasets": ["Candida albicans", "Aspergillus niger"], "loc1": [100, 100]}
    )
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    by_prop = dict(zip(result.reverse_table["Property"], result.reverse_table["Matching Species"]))
    assert "Candida albicans" in by_prop["Biofilm-formation"]
    assert "Aspergillus niger" not in by_prop["Biofilm-formation"]
    assert "Aspergillus niger" in by_prop["Thermophile"]


def test_filter_zero_weight_excludes_property(curated_df):
    inp = pd.DataFrame({"#Datasets": ["Candida albicans"], "loc1": [100]})
    weights = {"Biofilm-formation": 0, "Human-pathogenicity": 1, "Thermophile": 0}
    result = filter_fungi(inp, curated_df, weights, score_threshold=1, reads_threshold=10)
    # Only Human-pathogenicity contributes: A = 2*1 = 2, S = 1
    assert result.filtered["A-score"].iloc[0] == "2"
    assert result.filtered["S-score"].iloc[0] == "1"
    assert result.filtered["Contributing Properties"].iloc[0] == ["Human-pathogenicity"]


# ---------- unique_contributing_properties ----------

def test_unique_contributing_properties_empty():
    assert unique_contributing_properties(pd.DataFrame()) == []


def test_unique_contributing_properties_dedupes_and_sorts():
    df = pd.DataFrame(
        {"Contributing Properties": [["b", "a"], ["a", "c"], []]}
    )
    assert unique_contributing_properties(df) == ["a", "b", "c"]
