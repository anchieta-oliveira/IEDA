import os
import math
import numpy as np
from ieda.QM.overlap_matrix import OverlapMatrix


class TestOverlapMatrix:
    def test_default_creation(self):
        sm = OverlapMatrix()
        assert sm.title == ""
        assert sm.lower_half_triangle == []
        assert sm.matrix == []

    def test_empty_matrix(self):
        sm = OverlapMatrix(matrix=[])
        assert sm.matrix == []

    def test_normalization_factor(self):
        from ieda.QM.overlap_matrix import normalization_factor
        n = normalization_factor(1.0, 0, 0, 0)
        assert n > 0
        expected = (2.0 / math.pi) ** 0.75
        assert abs(n - expected) < 1e-10

    def test_boys_function_small_x(self):
        from ieda.QM.overlap_matrix import boys_function
        b = boys_function(0, 1e-10)
        assert abs(b - 1.0) < 1e-6

    def test_boys_function_large_x(self):
        from ieda.QM.overlap_matrix import boys_function
        b = boys_function(0, 2.0)
        assert b > 0

    def test_generate_angular_momentum(self):
        from ieda.QM.overlap_matrix import generate_angular_momentum
        result = generate_angular_momentum(2)
        assert len(result) >= 1

    def test_read_from_multiwfn(self):
        content = "     1\n 1  0.5\n     2\n 1  0.1\n 2  0.2\n"
        path = "/tmp/_test_intmat.txt"
        with open(path, "w") as f:
            f.write(content)
        sm = OverlapMatrix()
        try:
            result = sm.read_from_multiwfn(path)
            expected = [[0.5, 0.1],
                        [0.1, 0.2]]
            assert np.allclose(result, expected)
        finally:
            os.unlink(path)
