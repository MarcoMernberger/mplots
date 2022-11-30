# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for mplots.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import pathlib
import sys

root = pathlib.Path(".").parent.parent
sys.path.append(str(root / "src"))
import pytest
import pandas as pd


@pytest.fixture
def test_sample_frame():
    df = pd.DataFrame(
        {
            "sample1": [1, 1000, 15],
            "sample2": [2, 1000, 20],
            "sample3": [-2, 5000, -30],
            "sample4": [-3, 6000, -20],
        }
    )
    df.index = ["gen1", "gen2", "gen3"]
    return df


@pytest.fixture
def test_fc_frame():
    df = pd.DataFrame({"log2FC": [-1, 1, 15, 0.5], "fdr": [0.003, 0.05, 0.02, 0.5]})
    df.index = ["gen1", "gen2", "gen3", "gen4"]
    return df
