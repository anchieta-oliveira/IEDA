import numpy as np
import pytest
from ieda.MOL.selection import Selection


class TestSelection:
    def test_selection_all(self, mini_pdb):
        sel = Selection(selection="all", mol=mini_pdb)
        assert sel.result is not None
        assert len(sel.result.data) == 13

    def test_selection_chain_a(self, mini_pdb):
        sel = Selection(selection="chain A", mol=mini_pdb)
        assert len(sel.result.data) == 9

    def test_selection_chain_x(self, mini_pdb):
        sel = Selection(selection="chain X", mol=mini_pdb)
        assert len(sel.result.data) == 4

    def test_selection_resname_ala(self, mini_pdb):
        sel = Selection(selection="resname ALA", mol=mini_pdb)
        assert len(sel.result.data) == 5

    def test_selection_resname_gly(self, mini_pdb):
        sel = Selection(selection="resname GLY", mol=mini_pdb)
        assert len(sel.result.data) == 4

    def test_selection_resname_lig(self, mini_pdb):
        sel = Selection(selection="resname LIG", mol=mini_pdb)
        assert len(sel.result.data) == 4

    def test_selection_protein(self, mini_pdb):
        sel = Selection(selection="protein", mol=mini_pdb)
        assert len(sel.result.data) == 9

    def test_selection_within(self, mini_pdb):
        sel = Selection(selection="within 3 of resid 1", mol=mini_pdb)
        assert len(sel.result.data) >= 1

    def test_selection_not(self, mini_pdb):
        sel_all = Selection(selection="all", mol=mini_pdb)
        sel_not = Selection(selection="not chain X", mol=mini_pdb)
        assert len(sel_not.result.data) == len(sel_all.result.data) - 4

    def test_selection_and(self, mini_pdb):
        sel = Selection(selection="protein and chain A", mol=mini_pdb)
        assert len(sel.result.data) == 9

    def test_selection_or(self, mini_pdb):
        sel = Selection(selection="resname ALA or resname GLY", mol=mini_pdb)
        assert len(sel.result.data) == 9

    def test_selection_atom_name(self, mini_pdb):
        sel = Selection(selection="name CA", mol=mini_pdb)
        assert len(sel.result.data) == 3

    def test_selection_resid(self, mini_pdb):
        sel = Selection(selection="resid 1", mol=mini_pdb)
        assert len(sel.result.data) == 9

    def test_selection_beta_equal(self, mini_pdb):
        sel = Selection(selection="beta = 0", mol=mini_pdb)
        assert len(sel.result.data) == 13

    def test_selection_beta_greater(self, mini_pdb):
        sel = Selection(selection="beta > 0.5", mol=mini_pdb)
        assert len(sel.result.data) == 0

    def test_selection_beta_less(self, mini_pdb):
        sel = Selection(selection="beta <= 1", mol=mini_pdb)
        assert len(sel.result.data) == 13

    def test_selection_invalid_operator(self, mini_pdb):
        with pytest.raises(ValueError):
            Selection(selection="beta >> 1", mol=mini_pdb)

    def test_selection_none_mol(self):
        with pytest.raises((ValueError, AttributeError)):
            Selection(selection="all", mol=None)

    def test_selection_same(self, mini_pdb):
        sel = Selection(selection="same resid as chain X", mol=mini_pdb)
        assert len(sel.result.data) >= 1

    def test_selection_backbone(self, mini_pdb):
        sel = Selection(selection="backbone", mol=mini_pdb)
        for name in sel.result.data.names:
            assert name in ("N", "CA", "C", "O")
