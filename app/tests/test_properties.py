"""Tests for app/core/properties.py and viz_prep formatting helper."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from core.properties import (
    PROPERTY_CODE,
    PROPERTY_ORDER,
    code_legend,
    fingerprint_string,
    format_score,
    short_code,
)
from core.viz_prep import format_filtered_for_display


def test_short_code_known_property():
    assert short_code("Biofilm-formation") == "BF"
    assert short_code("antimicrobial-resistance") == "AMR"


def test_short_code_unknown_falls_through():
    assert short_code("Made-up-property") == "Made-up-property"


def test_code_legend_skips_unknown():
    legend = code_legend(["Biofilm-formation", "Made-up", "Thermophile"])
    assert "BF = Biofilm-formation" in legend
    assert "TH = Thermophile" in legend
    assert "Made-up" not in legend


def test_format_filtered_empty_returns_empty():
    assert format_filtered_for_display(pd.DataFrame()).empty
    assert format_filtered_for_display(None).empty


def test_format_filtered_renders_fingerprint_and_drops_text_props():
    df = pd.DataFrame(
        [
            {
                "#Datasets": "Candida albicans",
                "Score": "5.00",
                "Contributing Properties": ["Biofilm-formation", "Human-pathogenicity"],
                "Property Values": {
                    "antimicrobial-resistance": 0,
                    "Biofilm-formation": 2,
                    "Human-pathogenicity": 2,
                    "Radiation-resistance": 0,
                    "Spore-formation": 0,
                    "Thermophile": 1,
                },
                "Num loc": 2,
                "Locations": {"loc2": 1240, "loc1": 200},
            }
        ]
    )
    out = format_filtered_for_display(df)
    # Squares in canonical order: AMR BF HP RAD SF TH → 0 2 2 0 0 1 → ⬜🟥🟥⬜⬜🟨
    assert out["Property Values"].iloc[0] == "⬜🟥🟥⬜⬜🟨"
    # Contributing Properties text column dropped (squares already encode it).
    assert "Contributing Properties" not in out.columns
    # Sorted by key: loc1 before loc2.
    assert out["Locations"].iloc[0] == "loc1: 200, loc2: 1240"


def test_format_filtered_preserves_other_columns():
    df = pd.DataFrame(
        [
            {
                "#Datasets": "Candida albicans",
                "Score": "5.00",
                "Contributing Properties": [],
                "Num loc": 0,
                "Locations": {},
            }
        ]
    )
    out = format_filtered_for_display(df)
    assert out["Score"].iloc[0] == "5.00"
    assert out["Num loc"].iloc[0] == 0


def test_format_score_int_when_whole():
    assert format_score(5) == "5"
    assert format_score(5.0) == "5"
    assert format_score(0) == "0"


def test_format_score_decimal_when_fractional():
    assert format_score(2.5) == "2.50"
    assert format_score(4.333) == "4.33"


def test_fingerprint_canonical_order_red_yellow_white():
    values = {
        "antimicrobial-resistance": 1,
        "Biofilm-formation": 2,
        "Human-pathogenicity": 0,
        "Radiation-resistance": 2,
        "Spore-formation": 0,
        "Thermophile": 1,
    }
    # Order: AMR BF HP RAD SF TH → 1 2 0 2 0 1 → 🟨🟥⬜🟥⬜🟨
    assert fingerprint_string(values) == "🟨🟥⬜🟥⬜🟨"


def test_fingerprint_missing_keys_default_to_zero():
    assert fingerprint_string({}) == "⬜⬜⬜⬜⬜⬜"


def test_fingerprint_clamps_out_of_range_values():
    # 3 → clamps to 2 (red); -1 → clamps to 0 (white)
    values = {p: 3 for p in PROPERTY_ORDER}
    assert fingerprint_string(values) == "🟥🟥🟥🟥🟥🟥"
    values = {p: -1 for p in PROPERTY_ORDER}
    assert fingerprint_string(values) == "⬜⬜⬜⬜⬜⬜"


def test_fingerprint_rounds_fractional_values():
    # 1.4 → 1 (yellow); 1.6 → 2 (red).
    values = dict.fromkeys(PROPERTY_ORDER, 1.6)
    assert fingerprint_string(values) == "🟥" * 6


def test_property_code_covers_all_known_six():
    assert set(PROPERTY_CODE) == {
        "antimicrobial-resistance",
        "Biofilm-formation",
        "Human-pathogenicity",
        "Radiation-resistance",
        "Spore-formation",
        "Thermophile",
    }
