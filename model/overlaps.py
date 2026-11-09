"""Efficient search for nearby theoretical cluster masses."""

from bisect import bisect_left, bisect_right
from itertools import combinations
from math import comb

from .isotopes import isotope_pattern_general


class OverlapSearchLimitExceeded(RuntimeError):
    """Raised if the volume of overlap matches exceeds a safe output budget."""


def count_distinct_isotopologues(elements, max_counts, selected_isotopes):
    """Count unique isotopologues across all nonempty cluster compositions.

    For each element with k selected isotopes and maximum n atoms, the
    number of compositions including zero atoms is C(n + k, k).
    The all-zero cluster is excluded.
    """
    total = 1
    for element, maximum in zip(elements, max_counts):
        k = len(selected_isotopes[element])
        if k < 1:
            raise ValueError(f"No isotopes selected for {element}")
        total *= comb(maximum + k, k)
    return total - 1


def count_distinct_isotopologues_for_clusters(clusters, elements, selected_isotopes):
    """Count isotopologues of exactly the generated compositions.

    Unlike count_distinct_isotopologues (for a complete 0..maximum grid),
    this function works for a fixed-core series or any explicit cluster list.
    Each k-isotope element appearing n times has C(n+k-1, k-1)
    distinct assignments. The blank composition is excluded.
    """
    total = 0
    for cluster in clusters:
        if len(cluster) != len(elements):
            raise ValueError("The cluster tuple must match the element list.")
        if not any(cluster):
            continue
        count = 1
        for element, n_atoms in zip(elements, cluster):
            isotope_count = len(selected_isotopes[element])
            if isotope_count == 0:
                raise ValueError(f"No isotopes selected for {element}")
            count *= comb(n_atoms + isotope_count - 1, isotope_count - 1)
        total += count
    return total


def find_cluster_overlaps_general(
    clusters, elements, selected_isotopes, tolerance=0.05, *, max_matches=None,
    min_probability=0.0, return_unfiltered_counts=False,
):
    #Find theoretical mass-overlap pairs, optionally filtering weak isotopologues.

    if not 0.0 <= min_probability <= 1.0:
        raise ValueError("min_probability must be a fraction between 0 and 1.")
    if tolerance <= 0:
        raise ValueError("tolerance must be positive.")

    all_patterns = {}
    sorted_patterns = {}

    for cluster in clusters:
        composition = dict(zip(elements, cluster))
        pattern = isotope_pattern_general(composition, selected_isotopes)
        all_patterns[cluster] = pattern

        full_entries = sorted(
            ((mass, isotopes, probability)
             for isotopes, (mass, probability) in pattern.items()),
            key=lambda entry: entry[0],
        )
        kept_entries = [
            (mass, isotopes)
            for mass, isotopes, probability in full_entries
            if probability >= min_probability
        ]
        sorted_patterns[cluster] = (
            [mass for mass, _isotopes, _prob in full_entries],
            full_entries,
            [mass for mass, _isotopes in kept_entries],
            kept_entries,
        )

    overlapping_pairs = {}
    all_pair_counts = {}
    total_matches = 0

    for cluster1, cluster2 in combinations(clusters, 2):
        matches = []
        full_masses2, full_entries2, kept_masses2, kept_entries2 = (
            sorted_patterns[cluster2]
        )
        unfiltered_count = 0

        for isotopes1, (mass1, probability1) in all_patterns[cluster1].items():
            # When filtering, count the original theoretical matches separately,
            # so the UI can display "all" vs "retained" without constructing
            # potentially huge unfiltered match lists.
            if return_unfiltered_counts and min_probability > 0.0:
                raw_low = bisect_left(full_masses2, mass1 - tolerance)
                raw_high = bisect_right(full_masses2, mass1 + tolerance)
                unfiltered_count += sum(
                    abs(mass1 - mass2) < tolerance
                    for mass2, _isotopes2, _probability2
                    in full_entries2[raw_low:raw_high]
                )

            if probability1 < min_probability:
                continue

            low = bisect_left(kept_masses2, mass1 - tolerance)
            high = bisect_right(kept_masses2, mass1 + tolerance)

            for mass2, isotopes2 in kept_entries2[low:high]:
                difference = abs(mass1 - mass2)
                if difference < tolerance:
                    total_matches += 1
                    if max_matches is not None and total_matches > max_matches:
                        raise OverlapSearchLimitExceeded(
                            f"The search exceeds {max_matches:,} matching peak pairs."
                        )
                    matches.append(
                        (isotopes1, isotopes2, mass1, mass2, difference)
                    )

        if return_unfiltered_counts:
            # With no filter, the retained set *is* the full theoretical set.
            if min_probability == 0.0:
                unfiltered_count = len(matches)
            if unfiltered_count:
                all_pair_counts[(cluster1, cluster2)] = unfiltered_count

        if matches:
            overlapping_pairs[(cluster1, cluster2)] = matches

    if return_unfiltered_counts:
        return overlapping_pairs, all_pair_counts
    return overlapping_pairs


def group_overlap_regions(matches):
    sorted_matches = sorted(matches, key=lambda match: min(match[2], match[3]))
    regions = []
    for match in sorted_matches:
        mass1, mass2 = match[2], match[3]
        start, end = min(mass1, mass2), max(mass1, mass2)
        if not regions or start > regions[-1]["end"]:
            regions.append({"start": start, "end": end, "matches": [match]})
        else:
            regions[-1]["end"] = max(regions[-1]["end"], end)
            regions[-1]["matches"].append(match)
    return regions
