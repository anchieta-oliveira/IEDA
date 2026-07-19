import pytest
from ieda.inputParse import InputParser


class TestInputParser:
    def test_make_input_simple(self):
        parser = InputParser()
        cmd = "--matrix - pdb:test.pdb - qm:test.out - qm_sof:orca"
        result = parser.make_input(cmd)
        assert "matrix" in result
        assert result["matrix"]["pdb"] == "test.pdb"
        assert result["matrix"]["qm"] == "test.out"
        assert result["matrix"]["qm_sof"] == "orca"

    def test_make_input_nested(self):
        parser = InputParser()
        cmd = ("--two_sel - sel_a:chain A - sel_b:resname LIG "
               "- pdb:test.pdb - qm:test.out - qm_sof:orca")
        result = parser.make_input(cmd)
        assert "two_sel" in result
        assert result["two_sel"]["sel_a"] == "chain A"

    def test_make_input_with_subargs(self):
        parser = InputParser()
        cmd = ("--plot_heatmap - ied:data.npy - pdb:test.pdb "
               "+++intramol:False +++pep_bond:True")
        result = parser.make_input(cmd)
        assert "plot_heatmap" in result
        assert result["plot_heatmap"]["ied"] == "data.npy"

    def test_make_input_multiple_pipelines(self):
        parser = InputParser()
        cmd = ("--matrix - pdb:a.pdb - qm:a.out - qm_sof:orca "
               ": "
               "--two_sel - sel_a:chain A - sel_b:chain B - pdb:b.pdb "
               "- qm:b.out - qm_sof:molden")
        result = parser.make_input(cmd)
        assert "matrix" in result
        assert "two_sel" in result

    def test_file_args(self, tmp_path):
        config = tmp_path / "config.txt"
        config.write_text(
            "--matrix - pdb:test.pdb - qm:test.out - qm_sof:orca\n"
            "# this is a comment\n"
            "--two_sel - sel_a:all - sel_b:resname LIG - pdb:t.pdb - qm:t.out - qm_sof:xtb\n"
        )
        parser = InputParser()
        result = parser.file_args(str(config))
        assert "matrix" in result
        assert "two_sel" in result

    def test_empty_input(self):
        parser = InputParser()
        result = parser.make_input("")
        assert result == {}

    def test_comment_ignored(self):
        parser = InputParser()
        cmd = "--matrix - pdb:test.pdb - qm:test.out # this is a comment"
        result = parser.make_input(cmd)
        assert "matrix" in result
