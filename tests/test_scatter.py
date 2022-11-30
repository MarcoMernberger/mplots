#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
import unittest
from matplotlib.figure import Figure
from matplotlib.colors import Colormap
from mplots import scatter
from pandas import DataFrame

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def test_volcano_calc_columns(test_fc_frame):
    print(test_fc_frame.columns)
    df = scatter.volcano_calc(test_fc_frame, logFC_column="log2FC", fdr_column="fdr")
    assert isinstance(df, DataFrame)
    assert "group" in df.columns
    assert "logFC" in df.columns
    assert "-log10(p-value)" in df.columns


def test_volcano_calc_group(test_fc_frame):
    df = scatter.volcano_calc(test_fc_frame, logFC_column="log2FC", fdr_column="fdr")
    unittest.TestCase().assertListEqual(list(df["group"].values), ["blue", "red", "red", "grey"])


def test_volcano_plot(test_fc_frame):
    df = scatter.volcano_calc(test_fc_frame, logFC_column="log2FC", fdr_column="fdr")
    assert isinstance(df, DataFrame)
    fig = scatter.volcano_plot(df)
    print(type(fig))
    assert isinstance(fig, Figure)
