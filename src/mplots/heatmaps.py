#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""heatmaps.py: Contains some basic reusable heatmap functions."""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pandas import DataFrame
from matplotlib.figure import Figure
from matplotlib.colors import Colormap
from typing import Union


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def heatmap(
    df_plot: DataFrame,
    cmap: Union[str, Colormap] = "seismic",
    center: bool = True,
    shrink: float = 0.5,
    figsize=(10, 10),
    **kwargs
) -> Figure:
    f = plt.figure(figsize=figsize)
    if center:
        vmax = df_plot.abs().max().max()
        vmin = -vmax
        norm = matplotlib.colors.Normalize(vmin, vmax)
    else:
        norm = matplotlib.colors.Normalize()
    plt.imshow(df_plot, cmap=cmap, norm=norm, **kwargs)
    plt.yticks(np.arange(df_plot.shape[0]), df_plot.index.values)
    plt.xticks(np.arange(df_plot.shape[1]), df_plot.columns.values, rotation=75)
    plt.colorbar(shrink=shrink)
    return f
