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

    \mathrm{IED}_{AB}^{\mathrm{M}} = 2 \sum_{i}^{\mathrm{MO}_{\mathrm{occ}}} \sum_{r \in A} \sum_{s \in B} c_{ri}\,c_{si}\,S_{rs}.

where the summation runs over all occupied molecular orbitals (:math:`\mathrm{MO}_{\mathrm{occ}}`). For each occupied molecular orbital, the intermolecular electron density is obtained from the product 
of the coefficients of atomic orbital :math:`r` belonging to molecule :math:`A` and atomic orbital :math:`s` belonging to molecule :math:`B`, weighted by the corresponding overlap integral :math:`S_{rs}`. 
The prefactor of 2 accounts for doubly occupied molecular orbitals.

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