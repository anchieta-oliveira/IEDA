import os
import numpy as np
import pytest
from ieda.MOL.PDB import PDB


class TestPDB:
    def test_read_pdb(self, mini_pdb):
        assert mini_pdb.data is not None
        assert len(mini_pdb.data) == 13

    def test_atom_count(self, mini_pdb):
        assert len(mini_pdb.atoms) == 13

    def test_chain_ids(self, mini_pdb):
        chains = set(mini_pdb.data.chainids)
        assert chains == {"A", "X"}

    def test_residue_names(self, mini_pdb):
        resnames = set(mini_pdb.data.resnames)
        assert resnames == {"ALA", "GLY", "LIG"}

    def test_coordinates(self, mini_pdb):
        first = mini_pdb.data[0]
        assert first.xs == 1.0
        assert first.ys == 0.0
        assert first.zs == 0.0

    def test_write_pdb(self, mini_pdb, tmp_path):
        out = tmp_path / "test.pdb"
        mini_pdb.write(str(out))
        assert out.exists()
        reread = PDB(str(out))
        assert len(reread.data) == 13

    def test_get_center(self, mini_pdb):
        center = mini_pdb.get_center()
        assert len(center) == 3
        assert isinstance(center, np.ndarray)

    def test_get_distance_matrix(self, mini_pdb):
        dist = mini_pdb.get_distance_matrix()
        assert dist.shape == (13, 13)
        assert dist[0, 0] == 0.0
        assert dist[0, 1] > 0.0

    def test_from_pdb_string(self, mini_pdb):
        text = mini_pdb.get_text()
        assert "ATOM" in text

    def test_invalid_path(self):
        with pytest.raises((FileNotFoundError, ValueError)):
            PDB("nonexistent.pdb")

    def test_add_atoms(self, mini_pdb):
        from ieda.MOL.atom import Atom
        atom = Atom(id=14, name="N", resname="ALA", chain="A", resid=2,
                    coordinates=(5.0, 5.0, 5.0))
        mini_pdb.add_atoms([atom])
        assert len(mini_pdb.data) == 14

    def test_show_software_validation(self, mini_pdb):
        with pytest.raises(ValueError):
            mini_pdb.show(software="invalid")

    def test_move_center_to(self, mini_pdb):
        mini_pdb.move_center_to([0.0, 0.0, 0.0])

    def test_rotate_x(self, mini_pdb):
        mini_pdb.rotate(90, "x")
