Selection Keywords
==================

This section describes the available keywords used to create atom selections.
Selections allow users to filter atoms based on structural, chemical, or spatial
properties of a molecular system.

Residue-based Selections
------------------------

These keywords select atoms based on residue information.

- **resid**
  
  Select atoms by residue sequence number.
    E.g., ``resid 10`` selects atoms in residue number 10.

- **resname**
  
  Select atoms by residue name.
    E.g., ``resname ALA`` selects atoms in alanine residues.

- **residue**
  
  Select atoms by internal residue ID.
    E.g., ``residue 5`` selects atoms in the residue with internal ID 5 (sequential numbering).


Atom-based Selections
---------------------

These keywords select atoms using atom-level identifiers.

- **atom**
  
  Select atoms by atom serial number.
    E.g., ``atom 15`` selects the atom with serial number 15.

- **atomid**
  
  Select atoms by internal atom ID.
    E.g., ``atomid 20`` selects the atom with internal ID 20 (sequential numbering).

- **name**
  
  Select atoms by atom name.
    E.g., ``name CA`` selects atoms named CA (alpha carbons in proteins).


Chain and Segment Selections
----------------------------

These keywords select atoms based on structural grouping.

- **chain**
  
  Select atoms by chain ID.
    E.g., ``chain A`` selects atoms in chain A.

- **segname**
  
  Select atoms by segment name.
    E.g., ``segname PROT`` selects atoms in the segment named PROT.

Numeric Filters
---------------

These keywords filter atoms based on numerical properties.

- **beta**
  
  Select atoms using B-factor (temperature factor) values.
  Supports comparison operators such as ``>``, ``<``, ``>=``, ``<=``, and ``==``.
    E.g., ``beta > 0.5`` selects atoms with B-factor greater than 0.5.

- **occupancy**
  
  Select atoms based on occupancy values using comparison operators.
    E.g., ``occupancy <= 0.8`` selects atoms with occupancy less than or equal to 0.8.


Spatial Selection
-----------------

Selections based on geometric relationships.

- **within**
  
  Select atoms within a specified distance from another selection.
    E.g., ``within 5 of resid 10`` selects atoms within 5 Å of residue number 10.


Logical and Relational Operators
--------------------------------

These operators allow combining or modifying selections.

- **and**
  
  Combine selections using intersection.
    E.g., ``protein and resname ALA`` selects atoms that are both in the protein and in alanine residues.

- **or**
  
  Combine selections using union.
    E.g., ``resname ALA or segname ALIG`` selects atoms that are either in alanine residues or in the segment named ALIG.

- **not**
  
  Exclude atoms from a selection.
    E.g., ``not water`` selects all atoms that are not part of water molecules.

- **same**
  
  Select atoms sharing a specific property with another selection.
    E.g., ``same residue as within 5 of resid 10`` selects atoms in the same residue as those within 5 Å of residue number 10.


Biomolecular Categories
-----------------------

Predefined selections for common biomolecular groups.

- **protein**
  
  Select atoms belonging to protein residues.
    E.g., ``protein`` selects all atoms in protein residues.

- **backbone**
  
  Select protein backbone atoms (e.g., N, CA, C, O).
    E.g., ``backbone`` selects all backbone atoms in proteins.

- **water**
  
  Select water molecules.
    E.g., ``water`` selects all atoms in water molecules.

- **lipids**
  
  Select lipid molecules.
    E.g., ``lipids`` selects all atoms in lipid molecules.

- **nucleic**
  
  Select nucleic acid residues.
    E.g., ``nucleic`` selects all atoms in nucleic acid residues.


General Selection
-----------------

- **all**
  
  Select all atoms in the system.

