"""Streamlit entry point for the Cluster Mass Spectrum Simulator."""
import streamlit as st

from ui.cluster_generation import render_cluster_generation
from ui.isotope_selection import render_isotope_selection
from ui.overlap_selection import render_overlap_selection
from ui.visualization import render_visualization_section
from ui.labels import cluster_name

from pathlib import Path


st.set_page_config(
    page_title="ClusterSpec",
    page_icon="assets/clusterspec_icon.png",
    layout="wide"
)
st.title("ClusterSpec")

st.caption("Cluster & Molecular Mass Spectrum Simulator")


# 1. Choose elements and generate possible cluster compositions.
elements, max_counts, clusters = render_cluster_generation()

# 2. Choose which natural isotopes to include.
selected_isotopes = render_isotope_selection(elements)

# 3. Find and display pairs of potentially overlapping isotopologues.
render_overlap_selection(clusters, elements, max_counts, selected_isotopes)

# 4. Visualize individual isotope patterns or a hypothetical mixture.
render_visualization_section(
    clusters=clusters,
    elements=elements,
    selected_isotopes=selected_isotopes,
    cluster_name=cluster_name,
)
