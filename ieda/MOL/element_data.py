#!/usr/bin/env python3

# Description
###############################################################################

""" Module to define the data of the elements.

This module contains a dictionary with the data of the elements. The data is
stored in a dictionary, where the keys are the element symbols and the values
are dictionaries with the data of each element. The data includes the atomic
number, the name, the number of electrons, the number of neutrons, the mass,
and the atomic radii.

"""

# Imports
###############################################################################

# License
###############################################################################
'''
IEDA
Authors: José de Anchieta de Oliveira Filho
[The Federal University of Rio de Janeiro]
Carlos Chagas Filho Institute of Biophysics
Laboratory for Molecular Modeling and Dynamics
Av. Carlos Chagas Filho 373 - CCS - bloco G1-19,
Cidade Universitária - Rio de Janeiro, RJ, CEP: 21941-902
E-mail address: anchieta.oliveira@biof.ufrj.br
This project is licensed under the MIT License.
'''

# Classes
###############################################################################

# Functions
###############################################################################
## Private ##

## Public ##

# Aliases
###############################################################################

'''
data = {
    "H": {"atomic_number": 1, "name": "Hydrogen", "electrons": 1, "neutrons": 0, "mass": 1.008, "atomic_radii": 0.25},
    "C": {"atomic_number": 6, "name": "Carbon", "electrons": 6, "neutrons": 6, "mass": 12.011, "atomic_radii": 0.7},
    "N": {"atomic_number": 7, "name": "Nitrogen", "electrons": 7, "neutrons": 7, "mass": 14.007, "atomic_radii": 0.65},
    "O": {"atomic_number": 8, "name": "Oxygen", "electrons": 8, "neutrons": 8, "mass": 15.999, "atomic_radii": 0.6},
    "F": {"atomic_number": 9, "name": "Fluorine", "electrons": 9, "neutrons": 10, "mass": 18.998, "atomic_radii": .5},
    "Na": {"atomic_number": 11, "name": "Sodium", "electrons": 11, "neutrons": 12, "mass": 22.99, "atomic_radii": 1.8},
    "P": {"atomic_number": 15, "name": "Phosphorus", "electrons": 15, "neutrons": 16, "mass": 30.973, "atomic_radii": 1.0},
    "S": {"atomic_number": 16, "name": "Sulfur", "electrons": 16, "neutrons": 16, "mass": 32.06, "atomic_radii": 1.0},
    "Cl": {"atomic_number": 17, "name": "Chlorine", "electrons": 17, "neutrons": 18, "mass": 35.45, "atomic_radii": 1.0},
    "Zn": {"atomic_number": 30, "name": "Zinc", "electrons": 11, "neutrons": 12, "mass": 22.99, "atomic_radii": 1.35},
}
'''

data = {
    "H": {
        "atomic_number": 1,
        "name": "Hydrogen",
        "electrons": 1,
        "neutrons": 0,
        "mass": 1.008,
        "atomic_radii": 0.53,
        "energy_levels": (1,),
        "electronegativity": 2.20,
        "electron_affinity": 0.754
    },
    "C": {
        "atomic_number": 6,
        "name": "Carbon",
        "electrons": 6,
        "neutrons": 6,
        "mass": 12.011,
        "atomic_radii": 0.70,
        "energy_levels": (2, 4),
        "electronegativity": 2.55,
        "electron_affinity": 1.262
    },
    "N": {
        "atomic_number": 7,
        "name": "Nitrogen",
        "electrons": 7,
        "neutrons": 7,
        "mass": 14.007,
        "atomic_radii": 0.65,
        "energy_levels": (2, 5),
        "electronegativity": 3.04,
        "electron_affinity": -0.07
    },
    "O": {
        "atomic_number": 8,
        "name": "Oxygen",
        "electrons": 8,
        "neutrons": 8,
        "mass": 15.999,
        "atomic_radii": 0.60,
        "energy_levels": (2, 6),
        "electronegativity": 3.44,
        "electron_affinity": 1.461
    },
    "F": {
        "atomic_number": 9,
        "name": "Fluorine",
        "electrons": 9,
        "neutrons": 10,
        "mass": 18.998,
        "atomic_radii": 0.50,
        "energy_levels": (2, 7),
        "electronegativity": 3.98,
        "electron_affinity": 3.401
    },
    "Na": {
        "atomic_number": 11,
        "name": "Sodium",
        "electrons": 11,
        "neutrons": 12,
        "mass": 22.990,
        "atomic_radii": 1.80,
        "energy_levels": (2, 8, 1),
        "electronegativity": 0.93,
        "electron_affinity": 0.548
    },
    "P": {
        "atomic_number": 15,
        "name": "Phosphorus",
        "electrons": 15,
        "neutrons": 16,
        "mass": 30.974,
        "atomic_radii": 1.00,
        "energy_levels": (2, 8, 5),
        "electronegativity": 2.19,
        "electron_affinity": 0.746
    },
    "S": {
        "atomic_number": 16,
        "name": "Sulfur",
        "electrons": 16,
        "neutrons": 16,
        "mass": 32.06,
        "atomic_radii": 1.00,
        "energy_levels": (2, 8, 6),
        "electronegativity": 2.58,
        "electron_affinity": 2.077
    },
    "Cl": {
        "atomic_number": 17,
        "name": "Chlorine",
        "electrons": 17,
        "neutrons": 18,
        "mass": 35.45,
        "atomic_radii": 1.00,
        "energy_levels": (2, 8, 7),
        "electronegativity": 3.16,
        "electron_affinity": 3.612
    },
    "Zn": {
        "atomic_number": 30,
        "name": "Zinc",
        "electrons": 30,
        "neutrons": 35,
        "mass": 65.38,
        "atomic_radii": 1.35,
        "energy_levels": (2, 8, 18, 2),
        "electronegativity": 1.65,
        "electron_affinity": 0
    }
}
