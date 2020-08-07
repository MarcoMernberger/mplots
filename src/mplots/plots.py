#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""plots.py: Contains some routine plots."""
from mbf_genomics.genes import Genes
from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Union
from pypipegraph import Job
from pandas import DataFrame
from .pipegraph import MPPlotJob
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


def volcanoplot(
    genes_or_df: Union[Genes, DataFrame],
    fdr_column: str,
    logFC_column: str,
    significance_threshold: float = 0.05,
    fc_threhold: float = 1,
    outfile: Path = None,
    dependencies: List[Job] = [],
    **kwargs,
):
    """
    Draws a volcane plot for a given Genes object or dataframe.

    This returns a job that plots a default volcano plot that can take a
    Dataframe or Genes object as input. Column names for FDR and logFC columns
    must be supplied.
    Optional plot arguments can be supplied via kwargs, accepted so far are
    name, which is used as plot title and and figsize to adjust the figure size.

    Parameters
    ----------
    genes_or_df : Union[Genes, DataFrame]
        The data object that holds the numbers to plot.
    fdr_column : str
        The column name of the FDR values.
    logFC_column : str
        The column name of the logFC values.
    significance_threshold : float, optional
        The significance threshold, by default 0.05. A line is displayed at this
        value.
    fc_threhold : float, optional
        The fold change threshold, by default 1.
    outfile : Path, optional
        The path to the result file to be generated, by default None.
    dependencies : List[Job], optional
        List of job dependencies, by default [].

    Returns
    -------
    Job
        The FilegeneratingJob that creates the file.
    """
    deps = dependencies
    if isinstance(genes_or_df, DataFrame):
        default_name = "default"
        deps.append(genes_or_df.load())
        results_dir = Path("results/plots")
    else:
        default_name = f"{genes_or_df.name}"
        results_dir = genes_or_df.result_dir
    name = kwargs.get("name", default_name)
    figsize = kwargs.get("figsize", (10, 10))
    title = kwargs.get("title", name)
    if outfile is None:
        outfile = results_dir / f"{name}_volcano.svg"
    if isinstance(outfile, str):
        outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    calc_args = [fdr_column, logFC_column, significance_threshold, fc_threhold]

    def calc():
        if isinstance(genes_or_df, DataFrame):
            df = genes_or_df
        else:
            df = genes_or_df.df.copy()
        df = df.rename(
            {logFC_column: "log2 fold change", fdr_column: "-log10(p-value)"}
        )
        df["group"] = ["grey"] * len(df)
        df["group"][(df[logFC_column].values >= 1) & (df[fdr_column] <= 0.05)] = "red"
        df["group"][(df[logFC_column].values <= -1) & (df[fdr_column] <= 0.05)] = "blue"
        return df

    def plot(df):
        colors = ["grey", "red", "blue"]
        labels = {"grey": "non-sign.", "red": "up", "blue": "down"}
        fig = plt.figure(figsize=figsize)
        for group in colors:
            df_sub = df[df["group"] == group]
            plt.plot(
                df_sub[logFC_column].values,
                -np.log10(df_sub[fdr_column].values),
                ls="",
                marker="o",
                color=group,
                label=labels[group],
            )
        plt.axhline(-np.log10(0.05), color="lightgrey")
        plt.axvline(-1, color="lightgrey")
        plt.axvline(1, color="lightgrey")
        plt.ylabel(r"-log($p_{corrected}$)")
        plt.xlabel("log2FC")
        plt.legend()
        plt.title(title)
        return fig

    return MPPlotJob(outfile, calc, plot, calc_args=calc_args, plot_args=sorted(list(kwargs.items()))).depends_on(deps)


def logFC_correlation(
    name: str,
    genes_or_df: Union[Genes, DataFrame],
    cond1: str,
    cond2: str,
    fc_col1: str,
    fc_col2: str,
    mangler: Callable,
    logfc: float = 1,
    fdr: float = .05,
    outfile: Path = None,
    dependencies: List[Job] = [],
    **kwargs,
):
    """
    Draws a correlation plot for all genes in the DataFrame.

    This returns a job that plots a simple correlation plot that can take a
    Dataframe or Genes object as input.
    LogFC thresholds are indicated with dashed lines and FDR is color coded.
    Optional plot arguments can be supplied via kwargs, accepted so far are
    name, which is used as plot title and and figsize to adjust the figure size.

    Parameters
    ----------
    genes_or_df : Union[Genes, DataFrame]
        The data object that holds the numbers to plot.
    cond1 : str
        Name of the first condition for axis labels.
    cond2 : str
        Name of the second condition for axis labels.
    fc_col1 : str
        logFC column name for first condition.
    fc_col2 : str
        logFC column name for second condition.
    mangler : Callable
        Mangler function for comparison column name for axis labels.
    logfc : float, optional
        logFC threshold, by default 1
    fdr : float, optional
        FDR threshold, by default .05
    outfile : Path, optional
        The path to the result file to be generated, by default None.
    dependencies : List[Job], optional
        List of job dependencies, by default [].

    Returns
    -------
    Job
        The FilegeneratingJob that creates the file.
    """
    deps = dependencies
    if isinstance(genes_or_df, DataFrame):
        default_name = "default"
        deps.append(genes_or_df.load())
        results_dir = Path("results/plots")
    else:
        default_name = f"{genes_or_df.name}"
        results_dir = genes_or_df.result_dir
    name = kwargs.get("name", default_name)
    figsize = kwargs.get("figsize", (10, 10))
    title = kwargs.get("title", f"Correlation of Log2FC on {name}\n(logFC threshold = {logfc}, FDR threshold = {fdr})")
    if outfile is None:
        outfile = results_dir / f"{name}_corr.svg"
    if isinstance(outfile, str):
        outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    calc_args = [fc_col1, fc_col2, logfc, fdr]

    def calc():
        if isinstance(genes_or_df, DataFrame):
            df = genes_or_df
        else:
            df = genes_or_df.df.copy()
        return df

    def plot(df):
        result_dir = Path("results/plots")
        result_dir.mkdir(parents=True, exist_ok=True)
        fig = plt.figure(figsize=figsize)
        rho, p = st.pearsonr(df[fc_col1], df[fc_col2])
        fdr_col1 = fc_col1.replace("log2FC", "FDR")
        fdr_col2 = fc_col2.replace("log2FC", "FDR")
        sig1 = (df[fc_col1].abs() >= logfc) & (df[fdr_col1].abs() < fdr) 
        sig2 = (df[fc_col2].abs() >= logfc) & (df[fdr_col2].abs() < fdr)
        sig_both = sig1 & sig2
        no_sig = ~ sig_both
        for sig, label, color in [
            (no_sig, "non-sig.", "grey"),
            (sig1, cond1, "b"),
            (sig2, cond2, "r"),
            (sig_both, "both", "purple")
        ]:
            genes_selected = df[sig]
            plt.plot(genes_selected[fc_col1], genes_selected[fc_col2], ls="", marker=".", color=color, label=label)
        plt.gca().axhline(1, ls=(0, (5, 10)), lw=.5, color="k")
        plt.gca().axhline(-1, ls=(0, (5, 10)), lw=.5, color="k")
        plt.gca().axvline(1, ls=(0, (5, 10)), lw=.5, color="k")
        plt.gca().axvline(-1, ls=(0, (5, 10)), lw=.5, color="k")
        plt.xlabel(mangler(fc_col1))
        plt.ylabel(mangler(fc_col2))
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.text(0.02, .96, f"Pearson R = {rho:.3f}\np = {p:.3f}", transform=plt.gca().transAxes)
        return fig

    return MPPlotJob(outfile, calc, plot, calc_args=calc_args, plot_args=sorted(list(kwargs.items()))).depends_on(deps)
