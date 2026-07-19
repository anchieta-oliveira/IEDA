import os
import numpy as np
import pytest
from ieda.pipeline import Pipeline


class TestPipeline:
    def test_init(self):
        pipe = Pipeline(verbose=False)
        assert pipe is not None

    def test_str_to_bool_true(self):
        pipe = Pipeline()
        assert pipe._Pipeline__str_to_bool("True")
        assert pipe._Pipeline__str_to_bool("true")
        assert pipe._Pipeline__str_to_bool("yes")
        assert pipe._Pipeline__str_to_bool("1")

    def test_str_to_bool_false(self):
        pipe = Pipeline()
        assert not pipe._Pipeline__str_to_bool("False")
        assert not pipe._Pipeline__str_to_bool("false")
        assert not pipe._Pipeline__str_to_bool("no")
        assert not pipe._Pipeline__str_to_bool("0")

    def test_read_df_npy(self, tmp_path):
        pipe = Pipeline()
        data = {"1": {"2": 0.5}}
        path = tmp_path / "test.json"
        import json
        with open(path, "w") as f:
            json.dump(data, f)
        df = pipe._Pipeline__read_df(str(path))
        assert df is not None

    def test_config_file_not_found(self):
        pipe = Pipeline()
        with pytest.raises((FileNotFoundError, ValueError)):
            pipe.config_file("nonexistent.txt")

    def test_matrix_with_empty_args(self):
        pipe = Pipeline(verbose=False)
        with pytest.raises(Exception):
            pipe.matrix()

    def test_map_3d_with_empty_args(self):
        pipe = Pipeline(verbose=False)
        with pytest.raises(Exception):
            pipe.map_3D()

    def test_two_sel_with_empty_args(self):
        pipe = Pipeline(verbose=False)
        with pytest.raises(Exception):
            pipe.two_sel()

    def test_plot_heatmap_with_empty_args(self):
        pipe = Pipeline(verbose=False)
        with pytest.raises(Exception):
            pipe.plot_heatmap()
