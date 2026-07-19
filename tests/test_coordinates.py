import numpy as np
import pytest
from ieda.MOL.coordinates import Coordinates


class TestCoordinates:
    def test_default_creation(self):
        coord = Coordinates()
        assert coord.x == 0.0
        assert coord.y == 0.0
        assert coord.z == 0.0

    def test_creation_with_values(self):
        coord = Coordinates(1.0, 2.0, 3.0)
        assert coord.x == 1.0
        assert coord.y == 2.0
        assert coord.z == 3.0

    def test_measure_distance(self):
        a = Coordinates(0.0, 0.0, 0.0)
        b = Coordinates(3.0, 4.0, 0.0)
        assert abs(a.measure_distance(b.x, b.y, b.z) - 5.0) < 1e-10

    def test_measure_distance_symmetry(self):
        a = Coordinates(1.0, 2.0, 3.0)
        b = Coordinates(4.0, 5.0, 6.0)
        d1 = a.measure_distance(b.x, b.y, b.z)
        d2 = b.measure_distance(a.x, a.y, a.z)
        assert abs(d1 - d2) < 1e-10

    def test_rotate_x(self):
        coord = Coordinates(0.0, 1.0, 0.0)
        coord.rotate("x", 90)
        assert abs(coord.z - 1.0) < 1e-10
        assert abs(coord.y) < 1e-10

    def test_rotate_y(self):
        coord = Coordinates(1.0, 0.0, 0.0)
        coord.rotate("y", 90)
        assert abs(coord.z + 1.0) < 1e-10
        assert abs(coord.x) < 1e-10

    def test_rotate_z(self):
        coord = Coordinates(0.0, 1.0, 0.0)
        coord.rotate("z", 90)
        assert abs(coord.x + 1.0) < 1e-10
        assert abs(coord.y) < 1e-10

    def test_move_to_distance_from(self):
        coord = Coordinates(0.0, 0.0, 0.0)
        coord.move_to_distance_from(5.0, 10.0, 0.0, 0.0)
        assert abs(coord.measure_distance(10.0, 0.0, 0.0) - 5.0) < 1e-6

    def test_get_set(self):
        coord = Coordinates(7.0, 8.0, 9.0)
        assert coord.get_x() == 7.0
        assert coord.get_y() == 8.0
        assert coord.get_z() == 9.0

    def test_get_tuple(self):
        coord = Coordinates(1.0, 2.0, 3.0)
        assert coord.get_tuple() == (1.0, 2.0, 3.0)

    def test_get_array(self):
        coord = Coordinates(1.0, 2.0, 3.0)
        arr = coord.get_array()
        assert np.array_equal(arr, [1.0, 2.0, 3.0])

    def test_set_coordinates(self):
        coord = Coordinates()
        coord.set_coordinates(4.0, 5.0, 6.0)
        assert coord.x == 4.0
        assert coord.y == 5.0
        assert coord.z == 6.0

    def test_move_to(self):
        coord = Coordinates(1.0, 1.0, 1.0)
        coord.move_to(7.0, 8.0, 9.0)
        assert coord.x == 7.0
        assert coord.y == 8.0
        assert coord.z == 9.0

    def test_invalid_axis_raises(self):
        coord = Coordinates()
        with pytest.raises(ValueError):
            coord.rotate("w", 90)
