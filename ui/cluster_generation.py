import periodictable as pt
import streamlit as st

from model.clusters import generate_clusters_general
from model.molecular_cores import (
    FormulaError, generate_molecular_families, parse_molecular_formula,
)
from model.isotopes import get_element_isotopes


def render_cluster_generation():
    st.header("1. Cluster generation")
    generation_mode = st.radio(
        "Generation mode",
        ["Element combinations", "Molecular core + related cluster families"],
        horizontal=True,
        index=None,  
    )

    if generation_mode is None:
        st.info("Choose a generation mode to begin.")
        st.stop()
    if generation_mode == "Molecular core + related cluster families":
        return render_molecular_families()
    return render_element_combinations()


def render_molecular_families():
    """Build bare-core, bare-metal, and core+metal families as requested."""
    st.caption(
        "Choose the families to generate. Enter a molecular core (such as "
        "C60) and one or two variable elements (such as Nb). As elsewhere "
        "in this simulator, displayed species are assumed singly charged +1 ions."
    )

    core_formula = st.text_input(
        "Fixed molecular core formula",
        value="",
        key="fixed_core_formula_v3",
        placeholder="e.g. C60, CH4, Fe2O3 (optional for bare-metal-only mode)",
    )

    # The user may leave the core blank when requesting only bare metal ions.
    core_composition = {}
    if core_formula.strip():
        try:
            core_composition = parse_molecular_formula(core_formula)
        except FormulaError as exc:
            st.warning(str(exc))
            st.stop()

        if sum(core_composition.values()) > 300:
            st.warning("The molecular core currently supports at most 300 atoms.")
            st.stop()
        for element in core_composition:
            try:
                available_isotopes = get_element_isotopes(element)
            except (ValueError, KeyError, AttributeError):
                available_isotopes = {}
            if not available_isotopes:
                st.warning(
                    f"No natural isotope information is available for '{element}'."
                )
                st.stop()

    available_additions = [
        element.symbol for element in pt.elements
        if element.number > 0
        and element.symbol not in core_composition
        and get_element_isotopes(element.symbol)
    ]

    # Remove stale choices if the user changes the molecular core.
    chosen_key = "core_extra_elements_v3"
    if chosen_key in st.session_state:
        st.session_state[chosen_key] = [
            element for element in st.session_state[chosen_key]
            if element in available_additions
        ][:2]

    added_elements = st.multiselect(
        "Variable element(s) (up to two)",
        options=available_additions,
        default=[],
        max_selections=2,
        key=chosen_key,
        placeholder="Choose Nb, Fe, O, ...",
    )

    added_max_counts = {}
    missing_maxima = []
    if added_elements:
        columns = st.columns(len(added_elements))
        for column, element in zip(columns, added_elements):
            with column:
                maximum = st.number_input(
                    f"Maximum added {element}",
                    min_value=1,
                    max_value=100,
                    value=1,  
                    step=1,
                    key=f"core_max_v3_{element}",
                    placeholder="Enter maximum",
                )
                if maximum is None:
                    missing_maxima.append(element)
                else:
                    added_max_counts[element] = int(maximum)

    st.write("**Include species:**")
    include_bare_core = st.checkbox(
        "Bare molecular core (e.g. C60+)",
        value=False,
        key="include_bare_core_v3",
    )
    include_bare_metals = st.checkbox(
        "Bare variable-element clusters (e.g. Nb+, Nb2+, Nb3+)",
        value=False,
        key="include_bare_metals_v3",
    )
    include_core_metals = st.checkbox(
        "Molecular core + variable-element clusters "
        "(e.g. C60Nb+, C60Nb2+)",
        value=False,
        key="include_core_metals_v3",
    )

    if not any((include_bare_core, include_bare_metals, include_core_metals)):
        st.info("Select at least one family of species to generate.")
        st.stop()

    if (include_bare_core or include_core_metals) and not core_composition:
        st.info("Enter a fixed molecular core for the selected families.")
        st.stop()

    if (include_bare_metals or include_core_metals) and not added_elements:
        st.info("Choose at least one variable element for the selected families.")
        st.stop()

    if (include_bare_metals or include_core_metals) and missing_maxima:
        st.info("Enter a maximum atom count for every selected variable element.")
        st.stop()

    if len(added_elements) == 2 and (include_bare_metals or include_core_metals):
        st.caption(
            "With two variable elements, the selected families also include "
            "mixed compositions (e.g. NbFe+ and C60NbFe+)."
        )

    elements, max_counts, clusters = generate_molecular_families(
        core_composition,
        added_max_counts,
        include_bare_core=include_bare_core,
        include_bare_additions=include_bare_metals,
        include_core_additions=include_core_metals,
    )
    if not clusters:
        st.info("No species were generated. Check the selected families and maxima.")
        st.stop()

    from ui.labels import cluster_name
    st.write("Number of generated species:", len(clusters))
    with st.expander("Generated molecular ions", expanded=len(clusters) <= 12):
        names = [cluster_name(cluster, elements) for cluster in clusters[:100]]
        st.write(", ".join(names))
        if len(clusters) > 100:
            st.caption(f"Showing the first 100 of {len(clusters):,} species.")
    return elements, max_counts, clusters


def render_element_combinations():

    # Elements with available natural isotope abundance data
    available_elements = [
        element.symbol
        for element in pt.elements
        if element.number > 0
        and get_element_isotopes(element.symbol)
    ]


    # Shared selection for search and periodic table
    if "chosen_elements" not in st.session_state:
        st.session_state.chosen_elements = []


    def toggle_element(symbol):

        selected = st.session_state.chosen_elements.copy()

        if symbol in selected:
            selected.remove(symbol)

        elif len(selected) < 3:
            selected.append(symbol)

        st.session_state.chosen_elements = selected


    # Option 1: Search box
    st.multiselect(
        "Search for elements:",
        options=available_elements,
        max_selections=3,
        key="chosen_elements",
        placeholder="Type an element symbol..."
    )


    # Option 2: Periodic table
    with st.expander("Or choose elements from the periodic table"):

        st.caption(
            "Click to select or deselect an element. "
            "Selected elements are highlighted."
        )

        periodic_rows = [
            ["H"] + [""] * 16 + ["He"],

            ["Li", "Be"] + [""] * 10
            + "B C N O F Ne".split(),

            ["Na", "Mg"] + [""] * 10
            + "Al Si P S Cl Ar".split(),

            "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr".split(),

            "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe".split(),

            "Cs Ba Ln Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn".split(),

            "Fr Ra An Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split(),

            [""] * 2
            + "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split()
            + [""],

            [""] * 2
            + "Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr".split()
            + [""]
        ]

        for row_index, row in enumerate(periodic_rows):

            if row_index == 7:
                st.caption("Lanthanides")

            if row_index == 8:
                st.caption("Actinides")

            columns = st.columns(18, gap="small")

            for i, symbol in enumerate(row):

                if not symbol:
                    continue

                if symbol in ("Ln", "An"):
                    columns[i].caption(symbol)
                    continue

                selected = symbol in st.session_state.chosen_elements

                disabled = (
                    symbol not in available_elements
                    or (
                        len(st.session_state.chosen_elements) >= 3
                        and not selected
                    )
                )

                columns[i].button(
                    symbol,
                    key=f"pt_{symbol}",
                    on_click=toggle_element,
                    args=(symbol,),
                    type="primary" if selected else "secondary",
                    disabled=disabled,
                    width="stretch"
                )

        st.caption(
            "Greyed-out elements have no available natural isotope "
            "distribution in the current model, or the selection is full."
        )


    elements = st.session_state.chosen_elements

    if not elements:
        st.info("Select up to three elements to begin.")
        st.stop()


    st.subheader("Maximum atom counts")

    max_counts = []

    columns = st.columns(len(elements))

    for column, element in zip(columns, elements):
        with column:
            count = st.number_input(
            f"Maximum {element}",
            min_value=0,
            max_value=10,
            value=1,
            step=1,
            placeholder="Enter maximum",
            key=f"max_{element}"
        )
            max_counts.append(count)

    if any(count is None for count in max_counts):
        st.info("Enter the maximum number of atoms for each selected element.")
        st.stop()

    clusters = generate_clusters_general(elements, max_counts)

    st.write("Number of possible clusters:", len(clusters))

    return elements, max_counts, clusters
