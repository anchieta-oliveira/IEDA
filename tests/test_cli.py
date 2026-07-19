import subprocess
import sys
import pytest


class TestCLI:
    def test_help_exit_code(self):
        result = subprocess.run(
            [sys.executable, "-m", "main", "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0

    def test_matrix_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "main", "matrix", "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0

    def test_two_sel_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "main", "two_sel", "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0

    def test_map_3d_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "main", "map_3D", "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0

    def test_plot_heatmap_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "main", "plot_heatmap", "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0

    def test_radial_distribution_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "main", "radial_distribution", "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0

    def test_config_file_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "main", "config_file", "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0
