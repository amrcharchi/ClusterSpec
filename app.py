"""Streamlit entry point for the Cluster Mass Spectrum Simulator."""
import streamlit as st

from ui.cluster_generation import render_cluster_generation
from ui.isotope_selection import render_isotope_selection
from ui.overlap_selection import render_overlap_selection
from ui.visualization import render_spectrum_visualization
from ui.labels import cluster_name

from pathlib import Path


st.set_page_config(
    page_title="ClusterSpec",
    page_icon="assets/clusterspec_icon.png",
    layout="wide"
)
st.image(
    "assets/clusterspec_logo.png" \
    "",
    width=550
)

st.caption("Cluster & Molecular Mass Spectrum Simulator")


# 1. Choose elements and generate possible cluster compositions.
elements, max_counts, clusters = render_cluster_generation()

# 2. Choose which natural isotopes to include.
selected_isotopes = render_isotope_selection(elements)

# 3. Find and display pairs of potentially overlapping isotopologues.
render_overlap_selection(clusters, elements, max_counts, selected_isotopes)

# 4. Visualize individual isotope patterns or a hypothetical mixture.
st.header("4. Mass spectrum visualization")
selected_clusters = st.multiselect(
    "Choose cluster spectra to display:",
    options=clusters,
    format_func=lambda cluster: cluster_name(cluster, elements),
    placeholder="Select one or more clusters...",
    key="plot_clusters"
)
R = st.slider(
    "Mass resolving power (R = m / FWHM)",
    min_value=500,
    max_value=10000,
    value=3000,
    step=100
)
render_spectrum_visualization(
    selected_clusters=selected_clusters,
    elements=elements,
    selected_isotopes=selected_isotopes,
    R=R,
    cluster_name=cluster_name
)
