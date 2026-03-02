Theory Overview
===============

IEDA is based on quantum chemical descriptors derived from:

- Mulliken overlap populations
- Orbital-based approaches

Purpose
-------

The method identifies regions of intermolecular electronic density (IED),
especially in non-covalent interactions such as:

- Hydrogen bonds
- π–π stacking

Applications
------------

- Binding affinity analysis
- Drug design
- Biomolecular interaction studies

Theory
-------

Considering the normalized Molecular Orbital (MO) of a diatomic molecule, written approximately as a Linear Combination of Atomic Orbitals (LCAO), :math:`\chi_r` and :math:`\chi_s` refer to the r and s electrons in the k and l atoms, respectively:

.. math::

    \Phi = c_{r,k} \chi_{r,k} + c_{s,l} \chi_{s,l}

where :math:`c` is the contribution coefficient of the atomic orbital to the molecular orbital and :math:`\chi` is the wave function that describes the respective atomic orbital for the r or s electrons.

Mulliken :sup:`1,2` developed a quantitative analysis of the division of the electronic population, called (a) gross atomic populations, (b) net atomic populations and (c) overlap populations.

These LCAO populations were extensively studied by Mulliken :sup:`1,3,4`, resulting in important contributions to the study of covalent bonds, bond order, bond energy :sup:`2`, liquid atomic charges, Mulliken's electronegativity scale, dipole moment :sup:`3`, among others.

Here, these functions have been rearranged to investigate intermolecular interactions by electronic density between non-covalently bonded molecules and the correlation between them.

For a general case, two molecules can be represented as a set of electrons associated with the molecules A and B, where the Intermolecular Electronic Density (IED) can be quantified with the Mulliken overlap populations of the electrons set A in B, as:

.. math::

    IED(A,B) = \sum_{i}^{MO_{occ}} \sum_{r \in A} \sum_{s \in B}
    \left| c_{ir,k} \, c_{is,l} \right| S_{r,s}

being the sum of each MO, considering only the MOs occupied (:math:`MO_{occ}`) by electrons.

For each MO, the electronic density between the atoms of molecule A and B is calculated using Mulliken's population analysis :sup:`1`, expressed as the sum of the product modulus between the coefficient of the r orbital belonging to A and the coefficient of the s orbital belonging to B, multiplied by the overlap integral :math:`S` of the r and s atomic orbitals.

Another strategy for probing intermolecular electronic interactions is formulated at the orbital level. In this approach, an orbital-based method constructed from orbital coefficients is introduced to quantify the electronic co-occupation between non-covalently bonded molecular fragments, providing an alternative population-level perspective to the Mulliken-based analysis.

Thus, by representing the electronic population at the orbital level through :math:`|c_i|^2`, where :math:`c_i` is the coefficient in Equation (2), one can quantify the electronic density between two non-covalently bonded molecules, represented by two sets of atoms A and B, corresponding to each molecule, with:

.. math::

    IED(A,B) = \sum_{i}^{MO_{occ}} \sum_{r \in A} \sum_{s \in B}
    c_{ir,k}^{2} \, c_{is,l}^{2}


References
----------

.. [1] Mulliken, R. S. Electronic Population Analysis on LCAO–MO Molecular Wave Functions. I. *J. Chem. Phys.* 1955.
   DOI: `10.1063/1.1740588 <https://doi.org/10.1063/1.1740588>`_

.. [2] Mulliken, R. S. Electronic Structures of Molecules XI. Electroaffinity, Molecular Orbitals and Dipole Moments. *J. Chem. Phys.* 1935.
   DOI: `10.1063/1.1749731 <https://doi.org/10.1063/1.1749731>`_

.. [3] Mulliken, R. S. Electronic Population Analysis II. Overlap Populations, Bond Orders, and Covalent Bond Energies. *J. Chem. Phys.* 1955.
   DOI: `10.1063/1.1740589 <https://doi.org/10.1063/1.1740589>`_

.. [4] Mulliken, R. S. Electronic Structures of Molecules XI. Electroaffinity, Molecular Orbitals and Dipole Moments. *J. Chem. Phys.* 1935.
   DOI: `10.1063/1.1749731 <https://doi.org/10.1063/1.1749731>`_