# -*- coding: utf-8 -*-


# The data in '_relative_atomic_masses' is licensed under the CC-SA license
# https://en.wikipedia.org/w/index.php?title=List_of_elements&oldid=700476748
_elements = (
    # [Symbol, Name, Relative atomic mass]
    ["H", "Hydrogen", "1.008"],
    ["He", "Helium", "4.002602(2)"],
    ["Li", "Lithium", "6.94"],
    ["Be", "Beryllium", "9.0121831(5)"],
    ["B", "Boron", "10.81"],
    ["C", "Carbon", "12.011"],
    ["N", "Nitrogen", "14.007"],
    ["O", "Oxygen", "15.999"],
    ["F", "Fluorine", "18.998403163(6)"],
    ["Ne", "Neon", "20.1797(6)"],
    ["Na", "Sodium", "22.98976928(2)"],
    ["Mg", "Magnesium", "24.305"],
    ["Al", "Aluminium", "26.9815385(7)"],
    ["Si", "Silicon", "28.085"],
    ["P", "Phosphorus", "30.973761998(5)"],
    ["S", "Sulfur", "32.06"],
    ["Cl", "Chlorine", "35.45"],
    ["Ar", "Argon", "39.948(1)"],
    ["K", "Potassium", "39.0983(1)"],
    ["Ca", "Calcium", "40.078(4)"],
    ["Sc", "Scandium", "44.955908(5)"],
    ["Ti", "Titanium", "47.867(1)"],
    ["V", "Vanadium", "50.9415(1)"],
    ["Cr", "Chromium", "51.9961(6)"],
    ["Mn", "Manganese", "54.938044(3)"],
    ["Fe", "Iron", "55.845(2)"],
    ["Co", "Cobalt", "58.933194(4)"],
    ["Ni", "Nickel", "58.6934(4)"],
    ["Cu", "Copper", "63.546(3)"],
    ["Zn", "Zinc", "65.38(2)"],
    ["Ga", "Gallium", "69.723(1)"],
    ["Ge", "Germanium", "72.630(8)"],
    ["As", "Arsenic", "74.921595(6)"],
    ["Se", "Selenium", "78.971(8)"],
    ["Br", "Bromine", "79.904"],
    ["Kr", "Krypton", "83.798(2)"],
    ["Rb", "Rubidium", "85.4678(3)"],
    ["Sr", "Strontium", "87.62(1)"],
    ["Y", "Yttrium", "88.90584(2)"],
    ["Zr", "Zirconium", "91.224(2)"],
    ["Nb", "Niobium", "92.90637(2)"],
    ["Mo", "Molybdenum", "95.95(1)"],
    ["Tc", "Technetium", "[98]"],
    ["Ru", "Ruthenium", "101.07(2)"],
    ["Rh", "Rhodium", "102.90550(2)"],
    ["Pd", "Palladium", "106.42(1)"],
    ["Ag", "Silver", "107.8682(2)"],
    ["Cd", "Cadmium", "112.414(4)"],
    ["In", "Indium", "114.818(1)"],
    ["Sn", "Tin", "118.710(7)"],
    ["Sb", "Antimony", "121.760(1)"],
    ["Te", "Tellurium", "127.60(3)"],
    ["I", "Iodine", "126.90447(3)"],
    ["Xe", "Xenon", "131.293(6)"],
    ["Cs", "Caesium", "132.90545196(6)"],
    ["Ba", "Barium", "137.327(7)"],
    ["La", "Lanthanum", "138.90547(7)"],
    ["Ce", "Cerium", "140.116(1)"],
    ["Pr", "Praseodymium", "140.90766(2)"],
    ["Nd", "Neodymium", "144.242(3)"],
    ["Pm", "Promethium", "[145]"],
    ["Sm", "Samarium", "150.36(2)"],
    ["Eu", "Europium", "151.964(1)"],
    ["Gd", "Gadolinium", "157.25(3)"],
    ["Tb", "Terbium", "158.92535(2)"],
    ["Dy", "Dysprosium", "162.500(1)"],
    ["Ho", "Holmium", "164.93033(2)"],
    ["Er", "Erbium", "167.259(3)"],
    ["Tm", "Thulium", "168.93422(2)"],
    ["Yb", "Ytterbium", "173.045(10)"],
    ["Lu", "Lutetium", "174.9668(1)"],
    ["Hf", "Hafnium", "178.49(2)"],
    ["Ta", "Tantalum", "180.94788(2)"],
    ["W", "Tungsten", "183.84(1)"],
    ["Re", "Rhenium", "186.207(1)"],
    ["Os", "Osmium", "190.23(3)"],
    ["Ir", "Iridium", "192.217(3)"],
    ["Pt", "Platinum", "195.084(9)"],
    ["Au", "Gold", "196.966569(5)"],
    ["Hg", "Mercury", "200.592(3)"],
    ["Tl", "Thallium", "204.38"],
    ["Pb", "Lead", "207.2(1)"],
    ["Bi", "Bismuth", "208.98040(1)"],
    ["Po", "Polonium", "[209]"],
    ["At", "Astatine", "[210]"],
    ["Rn", "Radon", "[222]"],
    ["Fr", "Francium", "[223]"],
    ["Ra", "Radium", "[226]"],
    ["Ac", "Actinium", "[227]"],
    ["Th", "Thorium", "232.0377(4)"],
    ["Pa", "Protactinium", "231.03588(2)"],
    ["U", "Uranium", "238.02891(3)"],
    ["Np", "Neptunium", "[237]"],
    ["Pu", "Plutonium", "[244]"],
    ["Am", "Americium", "[243]"],
    ["Cm", "Curium", "[247]"],
    ["Bk", "Berkelium", "[247]"],
    ["Cf", "Californium", "[251]"],
    ["Es", "Einsteinium", "[252]"],
    ["Fm", "Fermium", "[257]"],
    ["Md", "Mendelevium", "[258]"],
    ["No", "Nobelium", "[259]"],
    ["Lr", "Lawrencium", "[266]"],
    ["Rf", "Rutherfordium", "[267]"],
    ["Db", "Dubnium", "[268]"],
    ["Sg", "Seaborgium", "[269]"],
    ["Bh", "Bohrium", "[270]"],
    ["Hs", "Hassium", "[271]"],
    ["Mt", "Meitnerium", "[278]"],
    ["Ds", "Darmstadtium", "[281]"],
    ["Rg", "Roentgenium", "[282]"],
    ["Cn", "Copernicium", "[285]"],
    ["Nh", "Nihonium", "[286]"],
    ["Fl", "Flerovium", "[289]"],
    ["Mc", "Moscovium", "[290]"],
    ["Lv", "Livermorium", "[293]"],
    ["Ts", "Tennessine", "[294]"],
    ["Og", "Oganesson", "[294]"],
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
    for mass in tuple(n[2] for n in _elements):
        if mass.startswith("[") and mass.endswith("]"):
            yield float(mass[1:-1])
        elif "(" in mass:
            yield float(mass.split("(")[0])
        else:
            yield (float(mass))


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
