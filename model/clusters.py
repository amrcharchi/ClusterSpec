from itertools import product

def generate_clusters_general(elements, max_counts):

    if not 1 <= len(elements) <= 3:
        raise ValueError("Select between 1 and 3 elements.")

    if len(elements) != len(max_counts):
        raise ValueError("Each element needs a maximum atom count.")

    if len(set(elements)) != len(elements):
        raise ValueError("Elements must be different.")

    clusters = []

    for counts in product(
        *(range(n + 1) for n in max_counts)
    ):

        if sum(counts) == 0:
            continue

        clusters.append(counts)

    return clusters
