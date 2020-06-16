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
