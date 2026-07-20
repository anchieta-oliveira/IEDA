import pytest
import os
from ieda.MOL.PDB import PDB
from ieda.QM.MO import MO


@pytest.fixture
def mini_pdb():
    path = os.path.join(os.path.dirname(__file__), "data", "mini.pdb")
    return PDB(path=path)


@pytest.fixture
def mo_objects():
    return [
        MO(id=0, occupation=2.0),
        MO(id=2, occupation=0.0),
        MO(id=1, occupation=0.0),
    ]
