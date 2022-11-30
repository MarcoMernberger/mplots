#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
from matplotlib.figure import Figure
from matplotlib.colors import Colormap
from mplots import heatmap


__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def test_heatmap(test_sample_frame):
    fig = heatmap(test_sample_frame)
    assert isinstance(fig, Figure)


def test_heatmap_colormap(test_sample_frame):
    colormap = "hot"
    fig = heatmap(test_sample_frame, colormap)
    assert isinstance(fig, Figure)
    cdict = {
        "red": [(0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
        "green": [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.75, 1.0, 1.0), (1.0, 1.0, 1.0)],
        "blue": [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 1.0, 1.0)],
    }
    colormap = matplotlib.colors.LinearSegmentedColormap("test", cdict)
    fig = heatmap(test_sample_frame, colormap)
    assert isinstance(fig, Figure)
