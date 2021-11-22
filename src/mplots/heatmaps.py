#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""heatmaps.py: Contains some basic reusable heatmap functions."""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pandas import DataFrame
from matplotlib.figure import Figure
from matplotlib.colors import Colormap
from typing import Union, Tuple


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def heatmap(
    df_plot: DataFrame,
    cmap: Union[str, Colormap] = "seismic",
    center: bool = True,
    shrink: float = 0.5,
    figsize: Tuple[int, int] = (10, 10),
    **kwargs
) -> Figure:
    """
    A simple heatmap.

    Returns a simple heatmap that plots values of a given DataFrame. Figure size,
    colormap, legend shrinkage parameter can be supplied. additional keyword
    arguments are supplied to matplotlib's imshow function. X and Y labels
    correspond to DataFrame column and row index.

    Parameters
    ----------
    df_plot : DataFrame
        DataFrame with values to be plotted.
    cmap : Union[str, Colormap], optional
        The ColorMap to use, by default "seismic".
    center : bool, optional
        Center the colorbar on zero, by default True.
    shrink : float, optional
        Shrinkage for colorbar, by default 0.5.
    figsize : Tuple[int, int], optional
        Figure size in inches, by default (10, 10).

    Returns
    -------
    Figure
        The plotted figure.
    """
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
