"""Short-form codes and display helpers for property names."""

from __future__ import annotations

PROPERTY_CODE: dict[str, str] = {
    "antimicrobial-resistance": "AMR",
    "Biofilm-formation": "BF",
    "Human-pathogenicity": "HP",
    "Radiation-resistance": "RAD",
    "Spore-formation": "SF",
    "Thermophile": "TH",
}

# Canonical display order for the per-row property fingerprint.
PROPERTY_ORDER: list[str] = list(PROPERTY_CODE.keys())

# Unicode squares for the per-property fingerprint.
SQUARE_FOR_VALUE: dict[int, str] = {
    0: "⬜",  # absent
    1: "🟨",  # weak (≥35% identity)
    2: "🟥",  # strong (≥75% identity)
}


def short_code(name: str) -> str:
    """Return the canonical short code for a property name, or the name itself."""
    return PROPERTY_CODE.get(name, name)


def code_legend(properties: list[str]) -> str:
    """Render a one-line ``CODE = full name`` legend, omitting unknown properties."""
    parts = [f"{PROPERTY_CODE[p]} = {p}" for p in properties if p in PROPERTY_CODE]
    return "  ·  ".join(parts)


def fingerprint_string(values: dict[str, int], order: list[str] | None = None) -> str:
    """Render a per-property strength fingerprint as Unicode squares.

    ``values`` maps property name → ``0`` / ``1`` / ``2``. Missing keys
    default to ``0``. ``order`` controls left→right column order; defaults
    to :data:`PROPERTY_ORDER`. Unknown values are rounded to {0, 1, 2}.
    """
    if order is None:
        order = PROPERTY_ORDER
    squares = []
    for prop in order:
        v = int(round(float(values.get(prop, 0))))
        v = max(0, min(2, v))
        squares.append(SQUARE_FOR_VALUE[v])
    return "".join(squares)


def format_score(value: float) -> str:
    """Format a score as int when whole, two decimals otherwise."""
    f = float(value)
    return f"{int(f)}" if f.is_integer() else f"{f:.2f}"
