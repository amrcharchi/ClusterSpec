import numpy as np

from .isotopes import isotope_pattern_general

def gaussian_spectrum_general(
    cluster, elements, selected_isotopes, mass_axis, resolving_power
):

    
    composition = dict(zip(elements, cluster))

    # Calculate isotope masses and probabilities
    pattern = isotope_pattern_general(
        composition, selected_isotopes
    )

    spectrum = np.zeros_like(mass_axis, dtype=float)

    for isotopes, (mass, probability) in pattern.items():

        # Gaussian width from mass resolving power
        fwhm = mass / resolving_power
        sigma = fwhm / 2.355

        # Gaussian peak profile
        gaussian = (
            1 / (sigma * np.sqrt(2 * np.pi))
            * np.exp(-0.5 * ((mass_axis - mass) / sigma) ** 2)
        )

        spectrum += probability * gaussian

    return spectrum
