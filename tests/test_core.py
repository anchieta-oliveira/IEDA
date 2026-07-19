import numpy as np
import pytest
from ieda.core import (
    intermolecular_eletron_density_two_selection,
    matrix_intermolecular_eletron_density_numba,
    two_selection_type_ao,
    read_file,
)


class TestCoreFunctions:
    def test_two_selection_identity(self):
        ao_ids_a = np.array([0, 1], dtype=np.int32)
        ao_ids_b = np.array([2, 3], dtype=np.int32)
        mo_coeffs = np.array([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ], dtype=np.float64)
        s_matrix = np.eye(4, dtype=np.float64)
        result = intermolecular_eletron_density_two_selection(
            ao_ids_a, ao_ids_b, mo_coeffs, s_matrix
        )
        assert result == 0.0

    def test_two_selection_simple(self):
        ao_ids_a = np.array([0], dtype=np.int32)
        ao_ids_b = np.array([1], dtype=np.int32)
        mo_coeffs = np.array([
            [1.0, 1.0, 0.0, 0.0],
        ], dtype=np.float64)
        s_matrix = np.array([
            [1.0, 0.5, 0.0, 0.0],
            [0.5, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ], dtype=np.float64)
        result = intermolecular_eletron_density_two_selection(
            ao_ids_a, ao_ids_b, mo_coeffs, s_matrix
        )
        expected = 2.0 * (1.0 * 1.0 * 0.5)
        assert abs(result - expected) < 1e-10

    def test_two_selection_zero_overlap(self):
        ao_ids_a = np.array([0], dtype=np.int32)
        ao_ids_b = np.array([2], dtype=np.int32)
        mo_coeffs = np.array([
            [1.0, 0.0, 1.0, 0.0],
        ], dtype=np.float64)
        s_matrix = np.eye(4, dtype=np.float64)
        result = intermolecular_eletron_density_two_selection(
            ao_ids_a, ao_ids_b, mo_coeffs, s_matrix
        )
        expected = 2.0 * (1.0 * 1.0 * 0.0)
        assert abs(result - expected) < 1e-10

    def test_two_selection_symmetric(self):
        ao_ids_a = np.array([0, 1], dtype=np.int32)
        ao_ids_b = np.array([2, 3], dtype=np.int32)
        rng = np.random.RandomState(123)
        mo_coeffs = rng.randn(2, 4).astype(np.float64)
        s_matrix = np.eye(4, dtype=np.float64)
        s_matrix[0, 2] = s_matrix[2, 0] = 0.3
        result_ab = intermolecular_eletron_density_two_selection(
            ao_ids_a, ao_ids_b, mo_coeffs, s_matrix
        )
        result_ba = intermolecular_eletron_density_two_selection(
            ao_ids_b, ao_ids_a, mo_coeffs, s_matrix
        )
        assert abs(result_ab - result_ba) < 1e-10

    def test_matrix_numba_small(self):
        ats = np.array([1, 2], dtype=np.int32)
        ao_atomindex = np.array([1, 1, 2, 2], dtype=np.int32)
        mo_coeffs = np.array([
            [1.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.5],
        ], dtype=np.float64)
        s_matrix = np.eye(4, dtype=np.float64)
        result = matrix_intermolecular_eletron_density_numba(
            ats, ao_atomindex, mo_coeffs, s_matrix
        )
        assert result.shape == (2, 2)
        assert result[0, 1] == result[1, 0]

    def test_matrix_numba_diagonal_positive(self):
        ats = np.array([1, 2], dtype=np.int32)
        ao_atomindex = np.array([1, 1, 2, 2], dtype=np.int32)
        mo_coeffs = np.array([
            [0.5, 0.3, 0.2, 0.1],
            [0.1, 0.2, 0.3, 0.4],
        ], dtype=np.float64)
        s_matrix = np.eye(4, dtype=np.float64)
        result = matrix_intermolecular_eletron_density_numba(
            ats, ao_atomindex, mo_coeffs, s_matrix
        )
        assert np.all(result >= 0)

    def test_two_selection_type_ao_filtered(self):
        ao_ids_a = np.array([0], dtype=np.int32)
        ao_ids_b = np.array([1], dtype=np.int32)
        mo_coeffs = np.array([
            [1.0, 1.0, 0.0, 0.0],
        ], dtype=np.float64)
        s_matrix = np.array([
            [1.0, 0.5, 0.0, 0.0],
            [0.5, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ], dtype=np.float64)
        result = two_selection_type_ao(
            ao_ids_a, ao_ids_b, mo_coeffs, s_matrix
        )
        assert result.shape[1] == 4

    def test_read_file_txt(self, tmp_path):
        test_file = tmp_path / "test.txt"
        test_file.write_text("hello world")
        content = read_file(str(test_file))
        assert content == "hello world"
