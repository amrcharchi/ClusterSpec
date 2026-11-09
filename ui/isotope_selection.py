import streamlit as st

from model.isotopes import get_element_isotopes


def render_isotope_selection(elements):
    st.header("2. Isotope selection")

    if "isotope_table_version" not in st.session_state:
        st.session_state.isotope_table_version = 0


    # Called when a user changes an isotope dropdown
    def isotope_dropdown_changed():
        st.session_state.isotope_table_version += 1


    # Called when a user checks or unchecks an isotope in the table
    def isotope_table_changed(editor_key, row_lookup):

        edits = st.session_state[editor_key]["edited_rows"]

        for row_index, changes in edits.items():

            if "Included" not in changes:
                continue

            element, mass_number = row_lookup[int(row_index)]
            widget_key = f"isotopes_{element}"

            chosen = set(st.session_state[widget_key])

            if changes["Included"]:
                chosen.add(mass_number)
            else:
                chosen.discard(mass_number)

            st.session_state[widget_key] = sorted(chosen)

        # Refresh the table with the updated selections
        st.session_state.isotope_table_version += 1



    selected_isotopes = {}
    isotope_rows = []
    row_lookup = []

    for element in elements:

        isotope_data = get_element_isotopes(element)
        available_isotopes = sorted(isotope_data.keys())

        widget_key = f"isotopes_{element}"

        # Default: only 16O for oxygen, all isotopes otherwise
        if widget_key not in st.session_state:

            if element == "O" and 16 in available_isotopes:
                st.session_state[widget_key] = [16]
            else:
                st.session_state[widget_key] = available_isotopes

        chosen = st.multiselect(
            f"Select {element} isotopes:",
            options=available_isotopes,
            format_func=lambda n, e=element: f"{n}{e}",
            key=widget_key,
            on_change=isotope_dropdown_changed
        )

        selected_isotopes[element] = chosen

        for mass_number, (mass, abundance) in isotope_data.items():

            isotope_rows.append({
                "Element": element,
                "Isotope": f"{mass_number}{element}",
                "Exact mass (u)": mass,
                "Natural abundance (%)": abundance,
                "Included": mass_number in chosen
            })

            row_lookup.append((element, mass_number))


    st.subheader("Isotope information")

    editor_key = (
        f"isotope_editor_{'_'.join(elements)}_"
        f"{st.session_state.isotope_table_version}"
    )

    st.data_editor(
        isotope_rows,
        key=editor_key,
        width="stretch",
        hide_index=True,
        num_rows="fixed",
        disabled=[
            "Element",
            "Isotope",
            "Exact mass (u)",
            "Natural abundance (%)"
        ],
        column_config={
            "Exact mass (u)": st.column_config.NumberColumn(
                format="%.6f"
            ),
            "Natural abundance (%)": st.column_config.NumberColumn(
                format="%.4f"
            ),
            "Included": st.column_config.CheckboxColumn(
                "Included",
                help="Include this isotope in the simulation"
            )
        },
        on_change=isotope_table_changed,
        args=(editor_key, row_lookup)
    )

    # Require at least one isotope for each selected element
    for element in elements:
        if not selected_isotopes[element]:
            st.warning(f"Select at least one isotope of {element}.")
            st.stop()

    return selected_isotopes
