"""Interactive Plotly UpSet diagram.

Composes intersection bars (top), per-property set-size bars (left,
mirrored), and a dot matrix (bottom-right) into a single figure.

Pure plotting: data prep lives in :mod:`core.viz_prep`.
"""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .viz_prep import build_upset_combinations, build_upset_set_sizes

_MAX_MEMBERS_IN_HOVER = 12
_BAR_COLOR = "#2c7fb8"
_SET_BAR_COLOR = "#7fcdbb"
_DOT_ACTIVE = "#252525"
_DOT_INACTIVE = "#d9d9d9"


def _members_block(members: list[str]) -> str:
    if len(members) <= _MAX_MEMBERS_IN_HOVER:
        return "<br>".join(f"• {m}" for m in members)
    head = "<br>".join(f"• {m}" for m in members[:_MAX_MEMBERS_IN_HOVER])
    return f"{head}<br>… +{len(members) - _MAX_MEMBERS_IN_HOVER} more"


def make_upset_figure(indicators: pd.DataFrame) -> go.Figure:
    """Build a full UpSet figure (top bars + left set-size bars + dot matrix).

    ``indicators`` is the boolean Species × Property matrix from
    :func:`core.viz_prep.build_upset_indicators`. Properties are ordered
    largest set at top of the matrix.
    """
    combos = build_upset_combinations(indicators)
    set_sizes = build_upset_set_sizes(indicators)
    if combos.empty or set_sizes.empty:
        return go.Figure()

    # Largest set at top of matrix → reverse so plotly's bottom-up category
    # axis renders the largest at the top visually.
    y_order = list(reversed(list(set_sizes.index)))
    n_combos = len(combos)
    bar_x = list(range(n_combos))

    bar_hover = []
    for _, row in combos.iterrows():
        bar_hover.append(
            f"<b>{row['Species count']} species</b><br>"
            f"<b>Properties:</b> {', '.join(row['Properties'])}<br>"
            f"<b>Members:</b><br>{_members_block(row['Members'])}<extra></extra>"
        )

    inactive_x: list[int] = []
    inactive_y: list[str] = []
    active_x: list[int] = []
    active_y: list[str] = []
    active_hover: list[str] = []
    line_x: list = []
    line_y: list = []
    for combo_idx, combo_row in combos.iterrows():
        active_props = set(combo_row["Properties"])
        ys_in_combo = [p for p in y_order if p in active_props]
        for p in y_order:
            if p in active_props:
                active_x.append(combo_idx)
                active_y.append(p)
                active_hover.append(
                    f"<b>{p}</b><br>"
                    f"In a combination of {len(active_props)} "
                    f"propert{'y' if len(active_props) == 1 else 'ies'}<br>"
                    f"{combo_row['Species count']} species in this combo"
                    "<extra></extra>"
                )
            else:
                inactive_x.append(combo_idx)
                inactive_y.append(p)
        if len(ys_in_combo) > 1:
            line_x.extend([combo_idx, combo_idx, None])
            line_y.extend([ys_in_combo[0], ys_in_combo[-1], None])

    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.20, 0.80],
        row_heights=[0.45, 0.55],
        horizontal_spacing=0.02,
        vertical_spacing=0.04,
        specs=[[None, {"type": "xy"}], [{"type": "xy"}, {"type": "xy"}]],
    )

    fig.add_trace(
        go.Bar(
            x=bar_x,
            y=combos["Species count"].tolist(),
            marker_color=_BAR_COLOR,
            hovertemplate=bar_hover,
            showlegend=False,
            name="Combination size",
        ),
        row=1,
        col=2,
    )

    fig.add_trace(
        go.Bar(
            y=y_order,
            x=[int(set_sizes[p]) for p in y_order],
            orientation="h",
            marker_color=_SET_BAR_COLOR,
            hovertemplate="<b>%{y}</b><br>%{x} species<extra></extra>",
            showlegend=False,
            name="Set size",
        ),
        row=2,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=inactive_x,
            y=inactive_y,
            mode="markers",
            marker=dict(size=14, color=_DOT_INACTIVE),
            hoverinfo="skip",
            showlegend=False,
            name="Inactive",
        ),
        row=2,
        col=2,
    )
    if line_x:
        fig.add_trace(
            go.Scatter(
                x=line_x,
                y=line_y,
                mode="lines",
                line=dict(color=_DOT_ACTIVE, width=2),
                hoverinfo="skip",
                showlegend=False,
                name="Combo link",
            ),
            row=2,
            col=2,
        )
    fig.add_trace(
        go.Scatter(
            x=active_x,
            y=active_y,
            mode="markers",
            marker=dict(size=14, color=_DOT_ACTIVE),
            hovertemplate=active_hover,
            showlegend=False,
            name="Active",
        ),
        row=2,
        col=2,
    )

    fig.update_xaxes(visible=False, row=1, col=2,
                     range=[-0.5, n_combos - 0.5])
    fig.update_yaxes(title_text="Species count", row=1, col=2,
                     gridcolor="#eeeeee", zeroline=False)
    fig.update_xaxes(autorange="reversed", title_text="Set size",
                     gridcolor="#eeeeee", zeroline=False, row=2, col=1)
    fig.update_yaxes(visible=False, row=2, col=1,
                     categoryorder="array", categoryarray=y_order)
    fig.update_xaxes(visible=False, row=2, col=2,
                     range=[-0.5, n_combos - 0.5])
    fig.update_yaxes(showticklabels=False, row=2, col=2,
                     categoryorder="array", categoryarray=y_order)

    fig.update_layout(
        height=440,
        margin=dict(t=20, l=0, r=0, b=0),
        plot_bgcolor="white",
        hovermode="closest",
        bargap=0.25,
    )
    return fig
