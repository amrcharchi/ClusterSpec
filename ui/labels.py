from collections import Counter

def cluster_name(cluster, elements):

    return "".join(
        f"{element}{n if n > 1 else ''}"
        for element, n in zip(elements, cluster)
        if n > 0
    ) + "+"



def isotopologue_label(isotopes):

    counts = Counter(isotopes)
    parts = []

    for (element, mass_number), number in sorted(counts.items()):

        label = f"{mass_number}{element}"

        if number > 1:
            label += f"×{number}"

        parts.append(label)

    return " · ".join(parts)
