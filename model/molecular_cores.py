from collections import Counter
from itertools import product
import re

FORMULA_TOKEN = re.compile(r"[A-Z][a-z]?|[0-9]+|[()]")


class FormulaError(ValueError):
    """Raised when the molecular core is not a valid ordinary chemical formula."""


def parse_molecular_formula(formula):
    """Parse C60, CH4, Fe2O3, or grouped formulas like Ca(OH)2.
    """
    if not isinstance(formula, str) or not formula.strip():
        raise FormulaError("Enter a molecular core, for example C60 or CH4.")

    formula = formula.strip()
    tokens = []
    position = 0
    for match in FORMULA_TOKEN.finditer(formula):
        if match.start() != position:
            raise FormulaError(f"Unsupported formula syntax near '{formula[position:]}'.")
        tokens.append(match.group())
        position = match.end()
    if position != len(formula):
        raise FormulaError(f"Unsupported formula syntax near '{formula[position:]}'.")

    index = 0

    def positive_multiplier():
        nonlocal index
        if index < len(tokens) and tokens[index].isdigit():
            number = int(tokens[index])
            index += 1
            if number < 1:
                raise FormulaError("Atom counts must be positive integers.")
            return number
        return 1

    def parse_group(inside_parentheses=False):
        nonlocal index
        composition = Counter()
        order = []
        while index < len(tokens):
            token = tokens[index]
            if token == ")":
                if not inside_parentheses:
                    raise FormulaError("Unexpected closing parenthesis.")
                if not composition:
                    raise FormulaError("Empty parentheses are not allowed.")
                index += 1
                return composition, order
            if token == "(":
                index += 1
                subgroup, suborder = parse_group(True)
                multiplier = positive_multiplier()
                for symbol in suborder:
                    if symbol not in composition:
                        order.append(symbol)
                    composition[symbol] += multiplier * subgroup[symbol]
            elif re.fullmatch(r"[A-Z][a-z]?", token):
                index += 1
                multiplier = positive_multiplier()
                if token not in composition:
                    order.append(token)
                composition[token] += multiplier
            else:
                raise FormulaError(f"Unexpected token '{token}' in formula.")

        if inside_parentheses:
            raise FormulaError("Missing closing parenthesis.")
        return composition, order

    counts, ordered_elements = parse_group()
    if not counts:
        raise FormulaError("Enter at least one element in the molecular core.")
    return {symbol: counts[symbol] for symbol in ordered_elements}


def generate_core_series(core_composition, added_max_counts, include_bare_core=True):
    #Return (elements, max_counts, clusters) for a fixed-core series.

    if not core_composition or any(count <= 0 for count in core_composition.values()):
        raise ValueError("The fixed core must contain positive atom counts.")
    if set(core_composition) & set(added_max_counts):
        raise ValueError("Additional elements must differ from the core elements.")
    if any(not isinstance(count, int) or count < 0 for count in added_max_counts.values()):
        raise ValueError("Additional maximum atom counts must be nonnegative integers.")

    elements = list(core_composition) + list(added_max_counts)
    max_counts = list(core_composition.values()) + list(added_max_counts.values())
    fixed_counts = tuple(core_composition.values())
    generated = []
    for extra_counts in product(*(range(n + 1) for n in added_max_counts.values())):
        if not include_bare_core and not any(extra_counts):
            continue
        generated.append(fixed_counts + extra_counts)
    return elements, max_counts, generated



def generate_molecular_families(
    core_composition,
    added_max_counts,
    *,
    include_bare_core=False,
    include_bare_additions=False,
    include_core_additions=False,
):
    #Generate independently selectable core, bare atom, and attached families.
    
    core_composition = dict(core_composition or {})
    added_max_counts = dict(added_max_counts or {})

    if any(not isinstance(n, int) or isinstance(n, bool) or n < 1
           for n in core_composition.values()):
        raise ValueError("Core atom counts must be positive integers.")
    if any(not isinstance(n, int) or isinstance(n, bool) or n < 0
           for n in added_max_counts.values()):
        raise ValueError("Added-atom maximum counts must be nonnegative integers.")
    if set(core_composition) & set(added_max_counts):
        raise ValueError("Additional elements must differ from core elements.")

    needs_core = include_bare_core or include_core_additions
    needs_added = include_bare_additions or include_core_additions

    if needs_core and not core_composition:
        raise ValueError("Enter a molecular core to generate core-containing species.")
    if needs_added and not added_max_counts:
        raise ValueError("Choose at least one variable element for this family.")

    # The table of isotopes should only include elements used by the selected
    # families; a bare-Nb-only search needn't ask about C isotopes.
    core_elements = list(core_composition) if needs_core else []
    extra_elements = list(added_max_counts) if needs_added else []
    elements = core_elements + extra_elements
    maxima = ([core_composition[s] for s in core_elements] +
              [added_max_counts[s] for s in extra_elements])

    core_counts = tuple(core_composition[s] for s in core_elements)
    zeros_core = (0,) * len(core_elements)
    zeros_added = (0,) * len(extra_elements)
    species = []

    if include_bare_core:
        species.append(core_counts + zeros_added)

    if needs_added:
        combinations_of_added_atoms = [
            numbers
            for numbers in product(*(range(added_max_counts[s] + 1)
                                     for s in extra_elements))
            if any(numbers)
        ]
        if include_bare_additions:
            species.extend(zeros_core + numbers
                           for numbers in combinations_of_added_atoms)
        if include_core_additions:
            species.extend(core_counts + numbers
                           for numbers in combinations_of_added_atoms)

    return elements, maxima, species
