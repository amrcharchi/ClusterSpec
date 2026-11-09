"""Natural isotope data and probabilities of distinct cluster isotopologues."""

from collections import Counter
from itertools import combinations_with_replacement, product
from math import factorial, fsum, prod

import periodictable as pt


def get_element_isotopes(symbol):
    """Return {mass_number: (exact_mass, natural_abundance_percent)}."""
    element = pt.elements.symbol(symbol)
    isotope_data = {}
    for mass_number in element.isotopes:
        isotope = element[mass_number]
        mass, abundance = isotope.mass, isotope.abundance
        if mass is not None and abundance is not None and abundance > 0:
            isotope_data[mass_number] = (mass, abundance)
    return isotope_data


def isotope_pattern_general(composition, selected_isotopes):
    """Return {isotopologue: [mass, probability]} for a cluster.

    Instead of enumerating all ordered atom assignments (k**n for n atoms
    and k isotopes), visit each *distinct* isotopic composition once:
    C(n+k-1, k-1) possibilities. Account for equivalent permutations with
    the multinomial coefficient n! / (n1! * n2! * ...).

    The returned key and value shapes match the original model, including
    probabilities renormalized over the isotopes selected by the user.
    """
    per_element_options = []

    for element, n_atoms in composition.items():
        if n_atoms == 0:
            continue

        isotope_data = get_element_isotopes(element)
        chosen = sorted(selected_isotopes[element])
        if not chosen:
            raise ValueError(f"No isotopes selected for {element}")

        abundance_total = sum(isotope_data[n][1] for n in chosen)
        masses = {n: isotope_data[n][0] for n in chosen}
        probabilities = {
            n: isotope_data[n][1] / abundance_total for n in chosen
        }

        options = []
        for mass_numbers in combinations_with_replacement(chosen, n_atoms):
            counts = Counter(mass_numbers)
            multiplicity = factorial(n_atoms)
            for count in counts.values():
                multiplicity //= factorial(count)

            probability = float(multiplicity)
            for mass_number, count in counts.items():
                probability *= probabilities[mass_number] ** count

            key = tuple((element, n) for n in mass_numbers)
            mass = fsum(masses[n] for n in mass_numbers)
            options.append((key, mass, probability))

        per_element_options.append(options)

    combined = {}
    for parts in product(*per_element_options):
        key = tuple(sorted(isotope for part in parts for isotope in part[0]))
        mass = fsum(part[1] for part in parts)
        probability = prod(part[2] for part in parts)
        combined[key] = [mass, probability]

    return combined
