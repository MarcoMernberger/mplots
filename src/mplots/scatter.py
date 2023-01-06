#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""plots.py: Contains some routine plots."""
from mbf.genomics.genes import Genes
from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Union
from pypipegraph import Job
from pandas import DataFrame
from .jobs import MPPlotJob
import pandas as pd
import pypipegraph as ppg
import scipy.stats as st


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


import pandas as pd
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pypipegraph as ppg
from pathlib import Path
from matplotlib.figure import Figure


def volcano_calc(
    df: DataFrame,
    fc_threshold: float = 1,
    alpha: float = 0.05,
    logFC_column: str = "logFC",
    fdr_column: str = "p-value",
) -> DataFrame:
    """
    Prepares a given DataFrame for volcano plot.

    It adds a column 'group' to the dataframe, stratifying data points into
    significant and non-significant groups coded by color, renames the fold change
    column and fdr column.

    Parameters
    ----------
    df : DataFrame
        DataFrame with data points.
    fc_threshold : float, optional
        logFC threshold for meaningful regulation, by default 1.
    alpha : float, optional
        FDR threshold, by default 0.05.
    logFC_column : str, optional
        Column name of logFC column, by default "logFC".
    fdr_column : str, optional
        Column name of FDR column, by default "p-value".

    Returns
    -------
    DataFrame
        DataFrame with group variable containing the colors for the plot.
    """
    df = df.rename(columns={logFC_column: "logFC", fdr_column: "-log10(p-value)"})
    df["group"] = ["grey"] * len(df)
    df["group"][(df["logFC"].values >= fc_threshold) & (df["-log10(p-value)"] <= alpha)] = "red"
    df["group"][(df["logFC"].values <= -fc_threshold) & (df["-log10(p-value)"] <= alpha)] = "blue"
    df["logFC"] = -np.log10(df["logFC"])
    return df


def volcano_plot(
    df,
    logFC_column: str = "logFC",
    fdr_column: str = "-log10(p-value)",
    alpha: float = 0.05,
    fc_threshold: float = 1.0,
    **kwargs,
) -> Figure:
    """
    Plots a volcano plot.

    Plots a volcano plot and returns a matplotlib figure. It expects a DataFrame
    with a given log FC column, an fdr column and a group column.

    Parameters
    ----------
    df : _type_
        DataFrame with data points.
    logFC_column : str, optional
        Column name of logFC column, by default "logFC"
    fdr_column : str, optional
        Column name of FDR column, by default "-log10(p-value)"
    alpha : float, optional
        Threshold for FDR, by default 0.05.

    Returns
    -------
    Figure
        Matplotlib figure with volcano plot.
    """
    labels = kwargs.get("labels", {"grey": "non-sign.", "red": "up", "blue": "down"})
    figsize = kwargs.get("figsize", (10, 10))
    title = kwargs.get("title", "Volcano")
    fig = plt.figure(figsize=figsize)
    xlabel = kwargs.get("xlabel", logFC_column)
    ylabel = kwargs.get("ylabel", r"-log10($p_{corrected}$)")
    for color, df_sub in df.groupby("group"):
        plt.plot(
            df_sub[logFC_column].values,
            df_sub[fdr_column].values,
            ls="",
            marker="o",
            color=color,
            label=labels[color],
        )
    plt.axhline(-np.log10(alpha), color="lightgrey")
    plt.axvline(-fc_threshold, color="lightgrey")
    plt.axvline(fc_threshold, color="lightgrey")
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.legend()
    plt.title(title)
    return fig
