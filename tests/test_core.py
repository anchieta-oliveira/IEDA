import numpy as np
import pytest
from ieda.QM.aux import AUX
from ieda.core import (
    intermolecular_eletron_density_two_selection,
    matrix_intermolecular_eletron_density_numba,
    two_selection_type_ao,
    read_file,
)


def reference_ied_matrix(ats, ao_atomindex, mo_coefficients, s_matrix):
    """Original AO lookup strategy, retained only to verify result parity."""
    result = np.zeros((ats.size, ats.size))
    for iat_a, at_a in enumerate(ats):
        ao_ids_a = np.where(ao_atomindex == at_a)[0]
        for iat_b in range(iat_a, ats.size):
            ao_ids_b = np.where(ao_atomindex == ats[iat_b])[0]
            value = 0.0
            for mo in mo_coefficients:
                contribution = 0.0
                for ci in ao_ids_a:
                    for cj in ao_ids_b:
                        contribution += mo[ci] * mo[cj] * s_matrix[ci, cj]
                value += 2 * contribution
            result[iat_a, iat_b] = value
            result[iat_b, iat_a] = value
    return result


class TestCoreFunctions:
    def test_aux_reads_only_occupied_mos_for_ieda(self, tmp_path):
        aux_file = tmp_path / "minimal.aux"
        aux_file.write_text(
            "ATOM_SYMTYPE[002]=\n"
            " S S\n"
            "AO_ATOMINDEX[002]=\n"
            " 1 2\n"
            "OVERLAP_MATRIX[00003]=\n"
            "# Lower half triangle only\n"
            " 1.0 0.5 1.0\n"
            "LMO_VECTORS[00004]=\n"
            " 0.1 0.2 0.3 0.4\n"
            "LMO_ENERGY_LEVELS[002]=\n"
            " -1.0 0.5\n"
            "MOLECULAR_ORBITAL_OCCUPANCIES[00002]=\n"
            " 2.0 0.0\n"
            "END OF MOPAC FILE\n"
        )
        aux = AUX()
        aux.read_ieda_data(str(aux_file))

        assert aux.lmo_vectors == []
        assert aux.ieda_mo_coefficients.dtype == np.float32
        assert aux.ieda_mo_coefficients.flags.c_contiguous
        assert aux.ieda_s_matrix.flags.c_contiguous
        np.testing.assert_array_equal(aux.ao_atomindex, [1, 2])
        np.testing.assert_allclose(aux.ieda_mo_coefficients, [[0.1, 0.2]])
        np.testing.assert_allclose(aux.ieda_s_matrix, [[1.0, 0.5], [0.5, 1.0]])

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

    def test_two_selection_float32_matches_float64(self):
        rng = np.random.RandomState(9)
        ao_ids_a = np.array([0, 2], dtype=np.int32)
        ao_ids_b = np.array([1, 3], dtype=np.int32)
        mo_coeffs = rng.randn(5, 4)
        s_matrix = rng.randn(4, 4)
        s_matrix = (s_matrix + s_matrix.T) / 2

        expected = intermolecular_eletron_density_two_selection(
            ao_ids_a, ao_ids_b, mo_coeffs, s_matrix
        )
        result = intermolecular_eletron_density_two_selection(
            ao_ids_a, ao_ids_b, mo_coeffs.astype(np.float32), s_matrix.astype(np.float32)
        )

        np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-6)

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

    def test_matrix_numba_matches_previous_ao_lookup(self):
        rng = np.random.RandomState(42)
        ats = np.array([1, 2, 3], dtype=np.int32)
        ao_atomindex = np.array([1, 2, 1, 3, 2, 3], dtype=np.int32)
        mo_coeffs = rng.randn(4, 6).astype(np.float64)
        s_matrix = rng.randn(6, 6).astype(np.float64)
        s_matrix = (s_matrix + s_matrix.T) / 2

        expected = reference_ied_matrix(ats, ao_atomindex, mo_coeffs, s_matrix)
        result = matrix_intermolecular_eletron_density_numba(
            ats, ao_atomindex, mo_coeffs, s_matrix
        )

        np.testing.assert_array_equal(result, expected)

    def test_matrix_numba_float32_uses_float32_output(self):
        rng = np.random.RandomState(7)
        ats = np.array([1, 2, 3], dtype=np.int32)
        ao_atomindex = np.array([1, 2, 1, 3, 2, 3], dtype=np.int32)
        mo_coeffs = rng.randn(4, 6).astype(np.float32)
        s_matrix = rng.randn(6, 6).astype(np.float32)
        s_matrix = (s_matrix + s_matrix.T) / np.float32(2.0)

        result = matrix_intermolecular_eletron_density_numba(
            ats, ao_atomindex, mo_coeffs, s_matrix
        )
        expected = reference_ied_matrix(
            ats, ao_atomindex, mo_coeffs, s_matrix
        ).astype(np.float32)

        assert result.dtype == np.float32
        np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-6)

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
