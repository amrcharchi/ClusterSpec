from hashlib import sha256

import streamlit as st

from model.isotopes import isotope_pattern_general
from model.overlaps import (
    OverlapSearchLimitExceeded,
    count_distinct_isotopologues_for_clusters,
    find_cluster_overlaps_general,
    group_overlap_regions,
)
from ui.labels import cluster_name, isotopologue_label


def render_overlap_selection(clusters, elements, max_counts, selected_isotopes):
    # A selected count tuple means something different when the molecular
    # core or element order changes. Clear previous selections on any change
    # to the generated set so the plot never silently shows a different ion.
    context = (tuple(elements), tuple(clusters))
    if st.session_state.get("previous_cluster_context") != context:
        st.session_state["plot_clusters"] = []
        st.session_state["last_overlap_selection"] = None
        st.session_state["previous_cluster_context"] = context

    # Store the clusters selected for plotting
    if "plot_clusters" not in st.session_state:
        st.session_state.plot_clusters = []

    # Remove selections that are no longer valid
    st.session_state.plot_clusters = [
        cluster
        for cluster in st.session_state.plot_clusters
        if cluster in clusters
    ]

    st.header("3. Potential mass overlaps")

    tolerance = st.number_input(
        "Mass proximity threshold (Th)",
        min_value=0.001,
        max_value=0.500,
        value=0.050,
        step=0.010,
        format="%.3f"
    )

    # Optional reporting filter: probabilities are fractions of the selected
    # isotope distribution for each cluster. 
    probability_filter_enabled = st.checkbox(
        "Filter low-probability isotopologues",
        value=False,
        key="overlap_probability_filter_enabled",
    )
    minimum_probability_percent = 0.0
    if probability_filter_enabled:
        minimum_probability_percent = st.selectbox(
            "Minimum isotopologue probability (%)",
            options=(0.001, 0.01, 0.1, 1.0),
            index=1,
            format_func=lambda percentage: f"{percentage:g}%",
            key="overlap_minimum_probability_percent",
        )
        st.caption(
            "A peak-pair match is retained only if BOTH isotopologues have "
            "probability at least the selected percentage within their own "
            "cluster's isotope distribution. This is not an instrument "
            "detection limit."
        )
    else:
        st.caption("All theoretical isotopologues are included in overlap searching.")

    # The model takes a probability fraction, not a percentage.
    min_probability = minimum_probability_percent / 100.0

    # Count the actual generated compositions, including fixed-core series.
    distinct_count = count_distinct_isotopologues_for_clusters(
        clusters, elements, selected_isotopes
    )
    st.caption(
        f"Estimated distinct isotopologues across {len(clusters):,} clusters: "
        f"{distinct_count:,}."
    )

    # Protect the app against excessively large allocations without blocking the rest of the interface.
    MAX_ISOTOPOLOGUES = 100_000
    CONFIRMATION_THRESHOLD = 15_000
    MAX_REPORTED_MATCHES = 50_000

    if distinct_count > MAX_ISOTOPOLOGUES:
        st.warning(
            f"This selection has {distinct_count:,} distinct isotopologues. "
            f"Automatic overlap searching is limited to "
            f"{MAX_ISOTOPOLOGUES:,} to protect the app's memory. "
            "Reduce atom counts or selected isotopes to search overlaps. "
            "Spectrum visualization remains available below."
        )
        return

    if distinct_count > CONFIRMATION_THRESHOLD:
        selection_signature = (
            tuple(clusters),
            tuple(elements),
            tuple(max_counts),
            tuple(
                (element, tuple(sorted(selected_isotopes[element])))
                for element in elements
            ),
            float(tolerance),
            float(min_probability),
        )
        if st.session_state.get("approved_large_overlap_search") != selection_signature:
            st.warning(
                f"This search contains {distinct_count:,} distinct isotopologues "
                "and may take some time. Run it only if needed."
            )
            if st.button("Run this large overlap search"):
                st.session_state["approved_large_overlap_search"] = selection_signature
            else:
                return

    try:
        with st.spinner("Finding nearby isotope masses..."):
            overlaps, all_pair_counts = find_cluster_overlaps_general(
                clusters,
                elements,
                selected_isotopes,
                tolerance,
                max_matches=MAX_REPORTED_MATCHES,
                min_probability=min_probability,
                return_unfiltered_counts=True,
            )
    except OverlapSearchLimitExceeded:
        st.warning(
            f"More than {MAX_REPORTED_MATCHES:,} peak-pair overlaps were found. "
            "Narrow the mass proximity threshold, apply a probability "
            "filter, or reduce selected atom counts or isotopes to inspect "
            "a manageable set. No partial overlap results are shown."
        )
        return

    overlap_rows = []

    for (cluster1, cluster2), matches in overlaps.items():

        regions = group_overlap_regions(matches)

        row = {
            "Cluster 1": cluster_name(cluster1, elements),
            "Cluster 2": cluster_name(cluster2, elements),
            "Peak-pair matches": len(matches),
        }
        if probability_filter_enabled:
            row["All theoretical matches"] = all_pair_counts[(cluster1, cluster2)]
        row["Overlap regions"] = len(regions)
        row["Minimum Δm (Th)"] = round(
            min(match[4] for match in matches), 6
        )
        overlap_rows.append(row)


    total_unfiltered_matches = sum(all_pair_counts.values())
    total_filtered_matches = sum(len(matches) for matches in overlaps.values())

    if probability_filter_enabled:
        st.write(
            f"**Peak-pair matches (after filter):** {total_filtered_matches:,} "
            f"of {total_unfiltered_matches:,} theoretical matches"
        )
        st.write(
            f"**Overlapping cluster pairs (after filter):** {len(overlaps):,} "
            f"of {len(all_pair_counts):,}"
        )
    else:
        st.write("Number of overlapping cluster pairs:", len(overlaps))
        st.write(f"Total theoretical peak-pair matches: {total_unfiltered_matches:,}")

    st.caption(
        "A peak-pair match is a pair of theoretical isotopologues closer "
        "than the mass threshold; it is not necessarily a distinct or "
        "experimentally resolved spectral peak."
    )


    if overlap_rows:

        # Show up to six table rows before scrolling
        row_height = 35
        visible_rows = min(len(overlap_rows), 6)
        table_height = 40 + visible_rows * row_height

        overlap_pairs = list(overlaps.keys())

        # Use a short, stable key tied to the exact set of generated species.
        key_material = repr((
            tuple(clusters), tuple(elements), selected_isotopes,
            tolerance, min_probability,
        ))
        table_key = "overlap_" + sha256(key_material.encode()).hexdigest()[:20]

        table_event = st.dataframe(
            overlap_rows,
            width="stretch",
            hide_index=True,
            height=table_height,
            row_height=row_height,
            on_select="rerun",
            selection_mode="single-row",
            key=table_key
        )

        selected_rows = table_event.selection.rows

        if selected_rows:

            row_index = selected_rows[0]
            selected_pair = overlap_pairs[row_index]

            # Preserve the existing click-to-plot behavior
            selection_id = (table_key, row_index)

            if st.session_state.get("last_overlap_selection") != selection_id:

                st.session_state.plot_clusters = list(selected_pair)
                st.session_state.last_overlap_selection = selection_id

            cluster1, cluster2 = selected_pair

            # Group the matches belonging to this pair
            regions = group_overlap_regions(
                overlaps[selected_pair]
            )

            with st.expander("Isotopologue details"):

                region_index = st.selectbox(
                    "Choose a mass region:",
                    options=list(range(len(regions))),
                    format_func=lambda i: (
                        f"Region {i + 1}: "
                        f"{regions[i]['start']:.4f}–"
                        f"{regions[i]['end']:.4f} Th "
                        f"({len(regions[i]['matches'])} matches)"
                    ),
                    key=f"region_{elements}_{selected_pair}_{min_probability}"
                )

                region = regions[region_index]

                # Obtain the probabilities of each isotopologue
                pattern1 = isotope_pattern_general(
                    dict(zip(elements, cluster1)),
                    selected_isotopes
                )

                pattern2 = isotope_pattern_general(
                    dict(zip(elements, cluster2)),
                    selected_isotopes
                )

                detail_rows = []

                for iso1, iso2, mass1, mass2, difference in region["matches"]:

                    prob1 = pattern1[iso1][1]
                    prob2 = pattern2[iso2][1]

                    detail_rows.append({
                        "Isotopologue 1": isotopologue_label(iso1),
                        "Isotopologue 2": isotopologue_label(iso2),
                        "Mass 1 (Th)": round(mass1, 6),
                        "Mass 2 (Th)": round(mass2, 6),
                        "Δm (Th)": round(difference, 6),
                        "P1 (%)": f"{100 * prob1:.6g}",
                        "P2 (%)": f"{100 * prob2:.6g}"
                    })

                st.dataframe(
                    detail_rows,
                    width="stretch",
                    hide_index=True,
                    row_height=35,
                    height=40 + 35 * min(len(detail_rows), 6)
                )

                st.caption(
                    "Probabilities are relative to each cluster's "
                    "selected isotope distribution. Mass regions are "
                    "grouped using nearby theoretical peak positions, "
                    "not instrumental resolution."
                )

        else:
            st.session_state.last_overlap_selection = None
            st.caption(
                "Select a cluster pair to display its spectra "
                "and inspect the contributing isotopologues."
            )

    else:
        if probability_filter_enabled and total_unfiltered_matches > 0:
            st.info(
                "There are theoretical mass overlaps, but no peak pairs "
                "meet this probability threshold. Try a lower percentage "
                "or turn off the probability filter."
            )
        else:
            st.info("No nearby isotope peaks were found.")
