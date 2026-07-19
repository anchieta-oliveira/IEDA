import numpy as np
import pytest
from ieda.ieda import IEDA


class TestIEDA:
    def test_init(self):
        ieda = IEDA()
        assert ieda is not None

    def test_get_index_homo_lumo(self, mo_objects):
        ieda = IEDA()
        homo, lumo = ieda._IEDA__get_index_homo_lumo(mo_objects)
        assert homo == 0
        assert lumo == 1

    def test_get_index_homo_lumo_all_occupied(self, mo_objects):
        for mo in mo_objects:
            if mo.id == 2:
                mo.occupation = 1.0
        ieda = IEDA()
        homo, lumo = ieda._IEDA__get_index_homo_lumo(mo_objects)
        assert homo == 1
        assert lumo == 2

    def test_two_sel_with_synthetic_data(self):
        n_aos = 4
        rng = np.random.RandomState(42)
        mo_coeffs = rng.randn(3, n_aos).astype(np.float64)
        s_matrix = np.eye(n_aos, dtype=np.float64)
        ao_atomindex = np.array([1, 1, 2, 2], dtype=np.int32)
        ats_a = np.array([1], dtype=np.int32)
        ats_b = np.array([2], dtype=np.int32)

        ao_ids_a = np.array([i for i, ao_id in enumerate(ao_atomindex) if ao_id in ats_a])
        ao_ids_b = np.array([i for i, ao_id in enumerate(ao_atomindex) if ao_id in ats_b])

        from ieda.core import intermolecular_eletron_density_two_selection
        result = intermolecular_eletron_density_two_selection(
            ao_ids_a, ao_ids_b, mo_coeffs, s_matrix
        )
        assert isinstance(result, float)

    def test_arrays_to_nested_dict(self):
        ieda = IEDA()
        matrix = np.array([[1.0, 2.0], [3.0, 4.0]])
        ats = np.array([1, 2])
        result = ieda._IEDA__arrays_to_nested_dict(matrix, ats)
        assert result[1][1] == 1.0
        assert result[1][2] == 2.0
        assert result[2][2] == 4.0
