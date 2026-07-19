import numpy as np
import pytest
from ieda.MOL.atom import Atom


class TestAtom:
    def test_default_creation(self):
        atom = Atom()
        assert atom.id == 0
        assert atom.name == ""
        assert atom.element == ""
        assert isinstance(atom.coordinates, object)

    def test_creation_with_values(self):
        atom = Atom(id=1, name="CA", resname="ALA", chain="A",
                     resid=10, coordinates=(1.0, 2.0, 3.0))
        assert atom.id == 1
        assert atom.name == "CA"
        assert atom.resname == "ALA"
        assert atom.chain == "A"
        assert atom.resid == 10
        assert atom.coordinates.x == 1.0
        assert atom.coordinates.y == 2.0
        assert atom.coordinates.z == 3.0

    def test_element_explicit(self):
        atom = Atom(name="ZN", element="ZN")
        assert atom.element == "ZN"

    def test_atomic_number_guessed_from_name(self):
        atom = Atom(name="N")
        assert atom.atomic_number == 7

    def test_atomic_number_guessed_from_symbol(self):
        atom = Atom(name="CA", symbol="C")
        assert atom.atomic_number == 6

    def test_set_coordinates_method(self):
        atom = Atom()
        atom.set_coordinates(4.0, 5.0, 6.0)
        assert atom.coordinates.x == 4.0
        assert atom.coordinates.y == 5.0
        assert atom.coordinates.z == 6.0

    def test_coordinates_default(self):
        atom = Atom()
        assert atom.coordinates.x == 0.0
        assert atom.coordinates.y == 0.0
        assert atom.coordinates.z == 0.0

    def test_bfactor_default(self):
        atom = Atom()
        assert atom.bfactor == 0.0

    def test_chain_default(self):
        atom = Atom()
        assert atom.chain == "X"
