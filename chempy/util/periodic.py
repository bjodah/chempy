# -*- coding: utf-8 -*-


# The data in '_relative_atomic_masses' is licensed under the CC-SA license
# https://en.wikipedia.org/wiki/Standard_atomic_weight#List_of_atomic_weights
_elements = (
    # [Symbol, Name, Relative atomic mass, uncertainty]
    ["H", "Hydrogen", 1.008, 0.0],
    ["He", "Helium", 4.002602, 2e-06],
    ["Li", "Lithium", 6.94, 0.0],
    ["Be", "Beryllium", 9.0121831, 5e-07],
    ["B", "Boron", 10.81, 0.0],
    ["C", "Carbon", 12.011, 0.0],
    ["N", "Nitrogen", 14.007, 0.0],
    ["O", "Oxygen", 15.999, 0.0],
    ["F", "Fluorine", 18.998403163, 6e-09],
    ["Ne", "Neon", 20.1797, 0.0006],
    ["Na", "Sodium", 22.98976928, 2e-08],
    ["Mg", "Magnesium", 24.305, 0.0],
    ["Al", "Aluminium", 26.9815384, 3e-07],
    ["Si", "Silicon", 28.085, 0.0],
    ["P", "Phosphorus", 30.973761998, 5e-09],
    ["S", "Sulfur", 32.06, 0.0],
    ["Cl", "Chlorine", 35.45, 0.0],
    ["Ar", "Argon", 39.95, 0.0],
    ["K", "Potassium", 39.0983, 0.0001],
    ["Ca", "Calcium", 40.078, 0.004],
    ["Sc", "Scandium", 44.955908, 5e-06],
    ["Ti", "Titanium", 47.867, 0.001],
    ["V", "Vanadium", 50.9415, 0.0001],
    ["Cr", "Chromium", 51.9961, 0.0006],
    ["Mn", "Manganese", 54.938043, 2e-06],
    ["Fe", "Iron", 55.845, 0.002],
    ["Co", "Cobalt", 58.933194, 3e-06],
    ["Ni", "Nickel", 58.6934, 0.0004],
    ["Cu", "Copper", 63.546, 0.003],
    ["Zn", "Zinc", 65.38, 0.02],
    ["Ga", "Gallium", 69.723, 0.001],
    ["Ge", "Germanium", 72.63, 0.008],
    ["As", "Arsenic", 74.921595, 6e-06],
    ["Se", "Selenium", 78.971, 0.008],
    ["Br", "Bromine", 79.904, 0.0],
    ["Kr", "Krypton", 83.798, 0.002],
    ["Rb", "Rubidium", 85.4678, 0.0003],
    ["Sr", "Strontium", 87.62, 0.01],
    ["Y", "Yttrium", 88.90584, 1e-05],
    ["Zr", "Zirconium", 91.224, 0.002],
    ["Nb", "Niobium", 92.90637, 1e-05],
    ["Mo", "Molybdenum", 95.95, 0.01],
    ["Tc", "Technetium", "[98]", 0.0],
    ["Ru", "Ruthenium", 101.07, 0.02],
    ["Rh", "Rhodium", 102.90549, 2e-06],
    ["Pd", "Palladium", 106.42, 0.01],
    ["Ag", "Silver", 107.8682, 0.0002],
    ["Cd", "Cadmium", 112.414, 0.004],
    ["In", "Indium", 114.818, 0.001],
    ["Sn", "Tin", 118.71, 0.007],
    ["Sb", "Antimony", 121.76, 0.001],
    ["Te", "Tellurium", 127.6, 0.03],
    ["I", "Iodine", 126.90447, 3e-05],
    ["Xe", "Xenon", 131.293, 0.006],
    ["Cs", "Caesium", 132.90545196, 6e-08],
    ["Ba", "Barium", 137.327, 0.007],
    ["La", "Lanthanum", 138.90547, 7e-05],
    ["Ce", "Cerium", 140.116, 0.001],
    ["Pr", "Praseodymium", 140.90766, 1e-05],
    ["Nd", "Neodymium", 144.242, 0.003],
    ["Pm", "Promethium", "[145]", 0.0],
    ["Sm", "Samarium", 150.36, 0.02],
    ["Eu", "Europium", 151.964, 0.001],
    ["Gd", "Gadolinium", 157.25, 0.03],
    ["Tb", "Terbium", 158.925354, 8e-06],
    ["Dy", "Dysprosium", 162.5, 0.001],
    ["Ho", "Holmium", 164.930328, 7e-06],
    ["Er", "Erbium", 167.259, 0.003],
    ["Tm", "Thulium", 168.934218, 6e-06],
    ["Yb", "Ytterbium", 173.045, 0.010],
    ["Lu", "Lutetium", 174.9668, 0.0001],
    ["Hf", "Hafnium", 178.486, 0.006],
    ["Ta", "Tantalum", 180.94788, 2e-05],
    ["W", "Tungsten", 183.84, 0.01],
    ["Re", "Rhenium", 186.207, 0.001],
    ["Os", "Osmium", 190.23, 0.03],
    ["Ir", "Iridium", 192.217, 0.002],
    ["Pt", "Platinum", 195.084, 0.009],
    ["Au", "Gold", 196.966570, 4e-06],
    ["Hg", "Mercury", 200.592, 0.003],
    ["Tl", "Thallium", 204.38, 0.0],
    ["Pb", "Lead", 207.2, 1.1],
    ["Bi", "Bismuth", 208.9804, 1e-05],
    ["Po", "Polonium", "[209]", 0.0],
    ["At", "Astatine", "[210]", 0.0],
    ["Rn", "Radon", "[222]", 0.0],
    ["Fr", "Francium", "[223]", 0.0],
    ["Ra", "Radium", "[226]", 0.0],
    ["Ac", "Actinium", "[227]", 0.0],
    ["Th", "Thorium", 232.0377, 0.0004],
    ["Pa", "Protactinium", 231.03588, 1e-05],
    ["U", "Uranium", 238.02891, 3e-05],
    ["Np", "Neptunium", "[237]", 0.0],
    ["Pu", "Plutonium", "[244]", 0.0],
    ["Am", "Americium", "[243]", 0.0],
    ["Cm", "Curium", "[247]", 0.0],
    ["Bk", "Berkelium", "[247]", 0.0],
    ["Cf", "Californium", "[251]", 0.0],
    ["Es", "Einsteinium", "[252]", 0.0],
    ["Fm", "Fermium", "[257]", 0.0],
    ["Md", "Mendelevium", "[258]", 0.0],
    ["No", "Nobelium", "[259]", 0.0],
    ["Lr", "Lawrencium", "[266]", 0.0],
    ["Rf", "Rutherfordium", "[267]", 0.0],
    ["Db", "Dubnium", "[268]", 0.0],
    ["Sg", "Seaborgium", "[269]", 0.0],
    ["Bh", "Bohrium", "[270]", 0.0],
    ["Hs", "Hassium", "[271]", 0.0],
    ["Mt", "Meitnerium", "[278]", 0.0],
    ["Ds", "Darmstadtium", "[281]", 0.0],
    ["Rg", "Roentgenium", "[282]", 0.0],
    ["Cn", "Copernicium", "[285]", 0.0],
    ["Nh", "Nihonium", "[286]", 0.0],
    ["Fl", "Flerovium", "[289]", 0.0],
    ["Mc", "Moscovium", "[290]", 0.0],
    ["Lv", "Livermorium", "[293]", 0.0],
    ["Ts", "Tennessine", "[294]", 0.0],
    ["Og", "Oganesson", "[294]", 0.0],
)

symbols = tuple(n[0] for n in _elements)
names = tuple(n[1] for n in _elements)
lower_names = tuple(n[1].lower() for n in _elements)

period_lengths = (2, 8, 8, 18, 18, 32, 32)
accum_period_lengths = (2, 10, 18, 36, 54, 86, 118)

# icosagens, crystallogens, pnictogens, chalcogens, halogens
groups = {g: tuple(x - 18 + g for x in accum_period_lengths[1:]) for g in range(13, 18)}
groups[1] = (1,) + tuple(x + 1 for x in accum_period_lengths[:-1])  # alkali metals
groups[2] = tuple(x + 2 for x in accum_period_lengths[:-1])  # alkaline earth metals
groups[18] = accum_period_lengths  # noble gases


def atomic_number(name):
    """Provide atomic number for a given element

    Parameters
    ----------
    name: str
        Full name or chemical symbol of an element

    Returns
    -------
    int
        Atomic number
    """
    try:
        return symbols.index(name.capitalize()) + 1
    except ValueError:
        return lower_names.index(name.lower()) + 1


def _get_relative_atomic_masses():
    for mass in tuple(element[2] for element in _elements):
        yield float(mass[1:-1]) if str(mass).startswith("[") else float(mass)


relative_atomic_masses = tuple(_get_relative_atomic_masses())


def mass_from_composition(composition):
    """Calculates molecular mass from atomic weights

    Parameters
    ----------
    composition: dict
        Dictionary mapping int (atomic number) to int (coefficient)

    Returns
    -------
    float
        molecular weight in atomic mass units


    Notes
    -----
    Atomic number 0 denotes charge or "net electron defficiency"

    Examples
    --------
    >>> '%.2f' % mass_from_composition({0: -1, 1: 1, 8: 1})
    '17.01'
    """
    mass = 0.0
    for k, v in composition.items():
        if k == 0:  # electron
            mass -= v * 5.489e-4
        else:
            mass += v * relative_atomic_masses[k - 1]
    return mass
