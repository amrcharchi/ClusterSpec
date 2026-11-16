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



class MolecularElements(list):
    """Element order with display labels keyed by unique elemental composition.

    Counts remain ordinary tuples for isotope calculations and dictionary keys.
    On a composition collision, the first generated family supplies the label.
    """

    def __init__(self, elements):
        super().__init__(elements)
        self.species_labels = {}


def generate_molecular_families(
    core_composition,
    added_max_counts,
    *,
    include_bare_core=False,
    include_bare_additions=False,
    include_core_additions=False,
    include_core_additions_reactant=False,
    reactant_formula=None,
    max_reactant_molecules=None,
):
    """Return (elements, actual maxima, unique elemental count tuples).

    Reactant units are expanded before calculation. Display labels preserve
    the molecular grouping; coincident compositions are emitted only once.
    Reactant inputs are ignored when their family is disabled.
    """
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

    needs_core = (include_bare_core or include_core_additions
                  or include_core_additions_reactant)
    needs_added = (include_bare_additions or include_core_additions
                   or include_core_additions_reactant)
    if needs_core and not core_composition:
        raise ValueError("Enter a molecular core to generate core-containing species.")
    if needs_added and not added_max_counts:
        raise ValueError("Choose at least one variable element for this family.")

    reactant = {}
    if include_core_additions_reactant:
        reactant = parse_molecular_formula(reactant_formula)
        if (not isinstance(max_reactant_molecules, int)
                or isinstance(max_reactant_molecules, bool)
                or max_reactant_molecules < 1):
            raise ValueError("Maximum reactant molecules must be a positive integer.")

    core = core_composition if needs_core else {}
    additions = added_max_counts if needs_added else {}
    elements = MolecularElements(dict.fromkeys([*core, *additions, *reactant]))
    species = []

    def formula(composition):
        return "".join(f"{s}{n if n > 1 else ''}"
                       for s, n in composition.items() if n)

    def append(composition, label):
        counts = tuple(composition.get(s, 0) for s in elements)
        if counts not in elements.species_labels:
            species.append(counts)
            elements.species_labels[counts] = label + "+"

    if include_bare_core:
        append(core, formula(core))

    bare_species, attached_species, reactant_species = [], [], []
    for numbers in product(*(range(n + 1) for n in additions.values())):
        if not any(numbers):
            continue
        extra = dict(zip(additions, numbers))
        attached = Counter(core) + Counter(extra)
        label = formula(core) + formula(extra)
        if include_bare_additions:
            bare_species.append((extra, formula(extra)))
        if include_core_additions:
            attached_species.append((attached, label))
        if include_core_additions_reactant:
            for units in range(1, max_reactant_molecules + 1):
                total = attached.copy()
                total.update({s: n * units for s, n in reactant.items()})
                suffix = str(units) if units > 1 else ""
                reactant_species.append((total, f"{label}({reactant_formula.strip()}){suffix}"))

    for composition, label in bare_species + attached_species + reactant_species:
        append(composition, label)
    maxima = [max((counts[i] for counts in species), default=0)
              for i in range(len(elements))]
    return elements, maxima, species
