from ieda.QM.MO import MO


class TestMO:
    def test_default_creation(self):
        mo = MO()
        assert mo.id == 0
        assert mo.energy == 0.0
        assert mo.spin == ""
        assert mo.occupation == 0.0
        assert mo.coefficients == []
        assert mo.ao_number == []

    def test_creation_with_values(self):
        mo = MO(title="1a", id=1, energy=-0.5, spin="Alpha", occupation=2.0)
        assert mo.title == "1a"
        assert mo.id == 1
        assert mo.energy == -0.5
        assert mo.spin == "Alpha"
        assert mo.occupation == 2.0

    def test_coefficients_float_type(self):
        mo = MO()
        mo.coefficients = [0.5, 0.3, 0.1]
        assert isinstance(mo.coefficients[0], float)

    def test_occupation_float_type(self):
        mo = MO(occupation=2.0)
        assert isinstance(mo.occupation, float)
        assert mo.occupation == 2.0

    def test_occupation_zero(self):
        mo = MO(occupation=0.0)
        assert mo.occupation == 0.0

    def test_ao_number(self):
        mo = MO()
        mo.ao_number = [1, 2, 3, 4]
        assert len(mo.ao_number) == 4

    def test_symtype(self):
        mo = MO()
        mo.symtype = ["s", "p", "d"]
        assert mo.symtype == ["s", "p", "d"]

    def test_ao_atomindex(self):
        mo = MO()
        mo.ao_atomindex = [1, 1, 2, 2]
        assert len(mo.ao_atomindex) == 4
