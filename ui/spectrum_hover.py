"""Render a Plotly figure with isotopologue details computed on hover."""

import json
from pathlib import Path

from streamlit.components.v1 import html


_HOVER_SCRIPT = Path(__file__).with_suffix(".js").read_text(encoding="utf-8")


def render_hover_chart(fig, species, traces, mixture, resolving_power, y_axis_title):
    metadata = json.dumps({
        "species": species,
        "traces": traces,
        "mixture": mixture,
        "resolvingPower": resolving_power,
        "intensityLabel": y_axis_title.split(" (")[0],
        "unit": "%" if y_axis_title.endswith("(%)") else " 1/Th",
    }, allow_nan=False).replace("<", r"\u003c")

    fig.update_layout(hovermode="closest")
    fig.update_traces(hoverinfo="none")
    chart = fig.to_html(
        full_html=False,
        include_plotlyjs=True,
        config={
            "displayModeBar": True,
            "scrollZoom": False,
            "showAxisDragHandles": False,
        },
        post_script=_HOVER_SCRIPT.replace("__HOVER_METADATA__", metadata),
    )
    html(
        '<style>body { margin: 0; }</style>' + chart,
        height=fig.layout.height,
    )
