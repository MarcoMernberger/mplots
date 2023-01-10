#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""plots.py: Contains some routine plots."""


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


import numpy as np
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as grid
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt

from pandas import DataFrame
from matplotlib.figure import Figure
from matplotlib.colors import Colormap
from pypipegraph2 import Job
from .functions import default_cycler
from typing import Callable
from .functions import default_cycler
from cycler import Cycler

# import pypipegraph2 as ppg
# from mbf.genomics.genes import Genes
# from pathlib import Path
# from typing import Optional, Callable, List, Dict, Tuple, Union
# from .jobs import MPPlotJob
# from typing import Union, Tuple
# from pathlib import Path


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
    df["-log10(p-value)"] = -np.log10(df["-log10(p-value)"])
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
    fontsize = kwargs.get("fontsize", 10)
    fontsize_title = kwargs.get("fontsize_title", fontsize)
    fontsize_ticks = kwargs.get("fontsize_ticks", fontsize)
    fontsize_legend = kwargs.get("fontsize_legend", fontsize)
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
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.xticks(fontsize=fontsize_ticks)
    plt.yticks(fontsize=fontsize_ticks)
    plt.legend(fontsize=fontsize_legend)
    plt.title(title, fontsize=fontsize_title)
    return fig


def generate_dr_plot(
    df: DataFrame,
    title: str = None,
    class_label_column: str = None,
    label_function: Callable = lambda x: x,
    custom_cycler=None,
    **params,
) -> Figure:
    """
    This assumes that self.transformed_matrix is an array-like object with shape (n_samples, n_components)
    """
    fontsize_title = params.get("fontsize_title", 12)
    custom_order = params.get("custom_order", None)
    show_names = params.get("show_names", False)
    x_label = params.get("xlabel", None)
    y_label = params.get("ylabel", None)
    dpi = params.get("dpi", 100)
    fig_x = params.get("fig_x", 8)
    fig_y = params.get("fig_y", 8)
    marker_size = params.get("marker_size", None)
    mfc = params.get("mfc", None)
    if not isinstance(custom_cycler, Cycler):
        custom_cycler = default_cycler()
    f = plt.figure(figsize=(fig_x, fig_y), dpi=dpi, frameon=True, edgecolor="k", linewidth=2)
    axe = plt.gca()
    axe.set_prop_cycle(custom_cycler)
    columns_to_use = list(df.columns.values)
    if x_label is None:
        x_label = f"{columns_to_use[0]}"
    if y_label is None:
        y_label = f"{columns_to_use[1]}"
    class_labels = class_label_column in df.columns
    labels = "labels" in df.columns
    if class_labels:
        columns_to_use.remove(class_label_column)
        if custom_order is not None:
            df["custom_order"] = [
                custom_order.find(label) for label in df[class_label_column].values
            ]
            df = df.sort_values("custom_order")
    else:
        df[class_label_column] = [""] * len(df)
    if not labels:
        df["labels"] = df.index
    else:
        columns_to_use.remove("labels")
    dimensions = params.get("dimension", len(df.columns))
    if dimensions < 2:
        raise ValueError(
            f"No 2D projection possible with only {dimensions} components, set k >= 2."
        )
    if len(columns_to_use) > 2:
        columns_to_use = columns_to_use[:2]

    for label, df_sub in df.groupby(class_label_column):
        plt.plot(
            df_sub[columns_to_use[0]].values,
            df_sub[columns_to_use[1]].values,
            markersize=marker_size,
            mfc=mfc,
            alpha=0.8,
            label=label,
            linestyle="None",
        )
    if title is not None:
        plt.title(title, fontsize=fontsize_title)
    elif class_labels:
        plt.title("Transformed samples with classes", {"fontsize": fontsize_title})
    else:
        plt.title("Transformed samples without classes", {"fontsize": fontsize_title})
    xmin = df[columns_to_use[0]].values.min()
    ymin = df[columns_to_use[1]].values.min()
    xmax = df[columns_to_use[0]].values.max()
    ymax = df[columns_to_use[1]].values.max()
    plt.gca().set_xlim([1.3 * xmin, 1.3 * xmax])
    plt.gca().set_ylim([1.3 * ymin, 1.3 * ymax])
    plt.gca().set_xlabel(x_label)
    plt.gca().set_ylabel(y_label)
    if class_labels:
        plt.gca().legend(loc="best")
    if show_names:
        for i, row in df.iterrows():
            plt.annotate(
                label_function(row["labels"]),
                xy=(row[columns_to_use[0]], row[columns_to_use[1]]),
                xytext=(-1, 1),
                textcoords="offset points",
                ha="right",
                va="bottom",
                size=8,
            )
    return f


def plot_correlation(df, column_x, column_y, pearson=True, **kwargs):
    fontsize = kwargs.get("fontsize", 10)
    fontsize_title = kwargs.get("fontsize_title", fontsize)
    fontsize_ticks = kwargs.get("fontsize_ticks", fontsize)
    fontsize_legend = kwargs.get("fontsize_legend", fontsize)
    fontsize_label = kwargs.get("fontsize_label", fontsize)
    fontsize_text = kwargs.get("fontsize_text", fontsize)
    figsize = kwargs.get("figsize", (6, 6))
    xlabel = kwargs.get("xlabel", column_x)
    ylabel = kwargs.get("ylabel", column_y)
    title = kwargs.get(
        "title",
        f"Correlation of {column_x} and {column_y}",
    )
    fig = plt.figure(figsize=figsize)
    plt.plot(
        df[column_x],
        df[column_y],
        ls="",
        marker=".",
    )
    if pearson:
        rho, p = st.pearsonr(df[column_x], df[column_y])
    else:
        rho, p = st.spearmanr(df[column_x], df[column_y])
    plt.xlabel(xlabel, fontsize=fontsize_label)
    plt.ylabel(ylabel, fontsize=fontsize_label)
    plt.title(title, fontsize=fontsize_title)
    plt.tight_layout()
    plt.text(0.02, 0.94, f"Pearson R = {rho:.3f}\np = {p:.3f}", transform=plt.gca().transAxes)
    return fig
