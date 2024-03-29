#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""plots.py: Contains some routine plots."""
from mbf.genomics.genes import Genes
from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Union
from pypipegraph import Job
from pandas import DataFrame
import pandas as pd
import pypipegraph as ppg
import scipy.stats as st
import inspect
from abc import ABC, abstractmethod
from .jobs import MPPlotJob
from .scatter import volcano_calc, volcano_plot
from .newplot import style_wrapper, add_function_wrapper

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


def volcano(
    genes_or_df: Union[Genes, DataFrame],
    logFC_column: str,
    p_column: str,
    fdr_column: str,
    significance_threshold: float = 0.05,
    fc_threhold: float = 1.0,
    outfile: Optional[Path] = None,
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
    p_column : str
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
    calc_args = [p_column, logFC_column, significance_threshold, fc_threhold]
    if isinstance(genes_or_df, DataFrame):
        default_name = "default"
        deps.append(genes_or_df.load())
        results_dir = Path(kwargs.get("result_dir", "results/plots"))
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

    def calc():
        if isinstance(genes_or_df, DataFrame):
            df = genes_or_df
        else:
            df = genes_or_df.df.copy()
        df = volcano_calc(
            df,
            fc_threshold=fc_threhold,
            alpha=significance_threshold,
            logFC_column=logFC_column,
            p_column=p_column,
            fdr_column=fdr_column,
        )
        return df

    def plot(df):
        return volcano_plot(df, figsize=figsize, title=title)

    job = MPPlotJob(
        outfile,
        calc,
        plot,
        plot_args=sorted(list(kwargs.items())),
        calc_args=calc_args,
        dependencies=deps,
    )
    return job


def volcanoplot(
    genes_or_df: Union[Genes, DataFrame],
    fdr_column: str,
    logFC_column: str,
    significance_threshold: float = 0.05,
    fc_threhold: float = 1.0,
    outfile: Optional[Path] = None,
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
        results_dir = Path("results/plots")
    else:
        default_name = f"{genes_or_df.name}"
        results_dir = genes_or_df.result_dir
        deps.append(genes_or_df.load())
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
        df = df.rename({logFC_column: "log2 fold change", fdr_column: "-log10(p-value)"})
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

    return MPPlotJob(
        outfile, calc, plot, calc_args=calc_args, plot_args=sorted(list(kwargs.items()))
    ).depends_on(deps)


def logFC_correlation(
    name: str,
    genes_or_df: Union[Genes, DataFrame],
    cond1: str,
    cond2: str,
    fc_col1: str,
    fc_col2: str,
    mangler: Callable,
    logfc: float = 1,
    fdr: float = 0.05,
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
    title = kwargs.get(
        "title",
        f"Correlation of Log2FC on {name}\n(logFC threshold = {logfc}, FDR threshold = {fdr})",
    )
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
        no_sig = ~sig_both
        for sig, label, color in [
            (no_sig, "non-sig.", "grey"),
            (sig1, cond1, "b"),
            (sig2, cond2, "r"),
            (sig_both, "both", "purple"),
        ]:
            genes_selected = df[sig]
            plt.plot(
                genes_selected[fc_col1],
                genes_selected[fc_col2],
                ls="",
                marker=".",
                color=color,
                label=label,
            )
        plt.gca().axhline(1, ls=(0, (5, 10)), lw=0.5, color="k")
        plt.gca().axhline(-1, ls=(0, (5, 10)), lw=0.5, color="k")
        plt.gca().axvline(1, ls=(0, (5, 10)), lw=0.5, color="k")
        plt.gca().axvline(-1, ls=(0, (5, 10)), lw=0.5, color="k")
        plt.xlabel(mangler(fc_col1))
        plt.ylabel(mangler(fc_col2))
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.text(0.02, 0.96, f"Pearson R = {rho:.3f}\np = {p:.3f}", transform=plt.gca().transAxes)
        return fig

    return MPPlotJob(
        outfile, calc, plot, calc_args=calc_args, plot_args=sorted(list(kwargs.items()))
    ).depends_on(deps)


# a function to plot the differential genes in comparison
def plot_fc_scatter_vioin(
    df: DataFrame,
    fc_columns_dict: Dict[str, str],
    fdr_columns_dict: Dict[str, str],
    reference_column: str,
    alpha: float = 0.05,
):
    """Plots a scatter/violin plot with labels indicating values that are >0 or <0 in a given reference column"""
    genes_up_in_ref = df[df[reference_column] > 0].index
    genes_down_in_ref = df[df[reference_column] < 0].index

    def _make_join_index(df):
        df = df.reset_index()
        df["join"] = df[["gene_stable_id", "variable"]].agg(",".join, axis=1)
        df = df.set_index("join")
        return df

    rng = np.random.default_rng(seed=12345)
    sample_names = list(fc_columns_dict.values())
    df_plot = df.rename(columns=fc_columns_dict)[sample_names].melt(
        value_name="logFC", ignore_index=False
    )
    df_plot = _make_join_index(df_plot)
    other_df = df.rename(columns=fdr_columns_dict)[sample_names].melt(
        value_name="FDR", ignore_index=False
    )
    other_df = _make_join_index(other_df)
    df_plot = df_plot.join(other_df[["FDR"]])

    group_to_numbers = dict(zip(df_plot.variable.unique(), np.arange(len(fc_columns_dict)) + 1))
    df_plot["num"] = df_plot.variable.map(group_to_numbers)
    f = plt.figure(figsize=(10, 6))
    for label, index in [
        (f"up in {fc_columns_dict[reference_column]}", genes_up_in_ref),
        (f"down in {fc_columns_dict[reference_column]}", genes_down_in_ref),
    ]:
        df_sub = df_plot[df_plot["gene_stable_id"].isin(index)]
        color = next(plt.gca()._get_lines.prop_cycler)["color"]
        for suffix, marker, df_sub_fdr in [
            ("(*)", ".", df_sub[df_sub.FDR < alpha]),
            ("(n.s.)", "x", df_sub[df_sub.FDR >= alpha]),
        ]:
            xs = rng.uniform(low=-0.3, high=0.3, size=df_sub_fdr.shape[0]) + df_sub_fdr.num
            plt.scatter(xs, df_sub_fdr.logFC, marker=marker, label=f"{label} {suffix}", color=color)
    labels = list(fc_columns_dict.values())
    plt.violinplot(df[list(fc_columns_dict.keys())])
    plt.axhline(y=1, color="r", linestyle="-", alpha=0.3)
    plt.axhline(y=0, color="r", linestyle="-", alpha=0.3)
    plt.axhline(y=-1, color="r", linestyle="-", alpha=0.3)
    plt.xticks(np.arange(len(labels)) + 1, labels)
    plt.ylabel("log2FC")
    plt.legend()
    plt.close()
    return f


def plot_expression(df: DataFrame, tpms: List[str], logfc_column: str):
    """
    Plots expression values (TPM) and means per gene.
    Left are genes upregulated in logfc_column, right are down-regulated genes.
    """
    ll = len(df.index.unique()) / 5
    f, axes = plt.subplots(1, 2, figsize=(14, ll))
    ii = 0
    for tick_position, direction, df_difference in [
        ("left", "up", df[df[logfc_column] > 0]),
        ("right", "down", df[df[logfc_column] < 0]),
    ]:
        plt.sca(axes[ii])
        df_plot = df_difference[tpms + ["mean(SKI-)", "mean(SKI+)"]].melt(
            var_name="sample", value_name="TPM", ignore_index=False
        )
        mapping = dict(zip(df_difference.index.values, np.arange(df_difference.shape[0])))
        df_plot["y"] = df_plot.index.map(mapping)
        df_plot = df_plot.fillna(0)
        for off, gr, mean_col in [(-0.05, "SKI-", "mean(SKI-)"), (0.05, "SKI+", "mean(SKI+)")]:
            color = next(plt.gca()._get_lines.prop_cycler)["color"]
            df_sub = df_plot[df_plot["sample"].str.startswith(gr)]
            plt.scatter(df_sub.TPM, df_sub.y + off, color=color, label=gr, alpha=0.3, s=10)
            df_sub = df_plot[df_plot["sample"] == mean_col]
            plt.scatter(
                df_sub.TPM, df_sub.y + off, color=color, label=mean_col, marker="|", alpha=1, s=100
            )
            plt.gca().yaxis.set_ticks_position(tick_position)
            for yy in df_sub.y:
                plt.axhline(y=yy + 0.5, color="grey", linestyle="-", alpha=0.3, linewidth=0.5)
        plt.yticks(ticks=df_sub.y, labels=df_sub.index)
        plt.ylim([-1, df_sub.shape[0] + 1])
        plt.title(direction)
        plt.xlabel("TPM")
        plt.legend(loc="lower right")
        plt.ylabel("Gene")
        ii += 1
    return f


def plot_means(df, **kwargs):
    figsize = kwargs.get("figsize", (10, 10))
    means = (
        df_counts.apply(["mean", "min", "max"])
        .sort_values(axis=1, by="mean", ascending=False)
        .transpose()
    )
    fig = plt.figure(figsize=figsize)
    plt.plot(means, label=means.columns.values)
    plt.yscale("log")
    plt.legend()
    plt.xticks([])
    plt.tight_layout()
    return fig


def plot_boxplots(df, **kwargs):
    figsize = kwargs.pop("figsize", (10, 5))
    sym = kwargs.pop("sym", "b.")
    rotation = kwargs.pop("rotation", 75)
    ylabel = kwargs.pop("ylabel", None)
    title = kwargs.pop("title", None)
    flierprops = kwargs.pop("flierprops", {"markersize": 2})
    notch = kwargs.pop("notch", True)
    figure = plt.figure(figsize=figsize)
    plt.boxplot(df, patch_artist=True, labels=df.columns, sym=sym, flierprops=flierprops, **kwargs)
    plt.xticks(rotation=rotation)
    plt.ylabel(ylabel)
    plt.title(title)
    return figure


def plot_hist(df, bins=50, topx=6, **kwargs):
    topx = min(topx, len(df.columns))
    title = kwargs.pop("title", None)
    figsize = kwargs.pop("figsize", None)
    f = plt.figure(figsize=figsize)
    df[df.columns[:topx]].hist(bins=bins)
    plt.tight_layout()
    plt.title(title)
    return f


"""
Plots - wrapper for matplotlib plots

we need a plot function.
sometimes we need a calc function as well.
we need a way to extract keyword arguments from kwargs that are accepted by each method.

should be callable

__call__ should accept *kwargs
we need a way to handle default params
we can have same kwargs for different plotlib functions --> solution: set all the same and introduce special args

should now acceptable arguments
"""


class Plot(ABC):
    def __init__(self, name):
        self.name = name
        self.parameters = parameters

    @abstractmethod
    def __call__(self):
        pass

    @abstractmethod
    def print_parameter(self):
        ...

    def __handle_kwargs(self, kwargs, function) -> dict:
        args = inspect.getfullargspec(function).args
        new_kwargs = {key: kwargs[key] for key in args if key in kwargs}
        return new_kwargs

    @abstractmethod
    def calc(self, **kwargs):
        pass

    @abstractmethod
    def plot(self, **kwargs):
        pass

    def plot_func_arguments(self):
        pass


class Box(Plot):
    def __init__(self):
        super().__init__("Boxplot")

    def __call__(df, **kwargs):
        return figure

    def plot_boxplots(df, **kwargs):
        figsize = kwargs.get("figsize", (10, 5))
        sym = kwargs.get("sym", "b.")
        rotation = kwargs.get("rotation", 75)
        title = kwargs.get("title", None)
        flierprops = kwargs.get("flierprops", {"markersize": 2})
        notch = kwargs.get("notch", True)
        figure = plt.figure(figsize=figsize)
        plot_kwargs = __handle_kwargs(kwargs, plt.boxplot)
        plt.boxplot(
            df, patch_artist=True, labels=df.columns, sym=sym, flierprops=flierprops, **plot_kwargs
        )
        plt.xticks(rotation=rotation)
        plt.title(title)
        return figure

    def get_args(self, **kwargs):
        figsize = kwargs.get("figsize", (10, 5))
        fontsize = kwargs.get("fontsize", 10)
        fontsize_title = kwargs.get("fontsize_title", fontsize)
        fontsize_ticks = kwargs.get("fontsize_ticks", fontsize)
        fontsize_legend = kwargs.get("fontsize_legend", fontsize)
        title = kwargs.get("title", "Volcano")


# f = plot_boxplots((df_counts +1).apply("log2"), title="log2(raw count + 1)")


# kwargs = {"figsize": (10, 12), "bla": 3, "notch": True, "label": 12}
# print(__handle_kwargs(kwargs, plt.figure))
# print(__handle_kwargs(kwargs, plt.boxplot))
# print(__handle_kwargs(kwargs, plt.xlabel))
# inspect.getfullargspec(plt.plot)


def save_figure(f, folder, name, bbox_inches="tight"):
    folder.mkdir(exist_ok=True, parents=True)
    for suffix in [".png", ".svg", ".pdf"]:
        f.savefig(folder / (name + suffix), bbox_inches=bbox_inches)


@style_wrapper
@add_function_wrapper
def plot_empty(*args, **kwargs):
    if len(args) > 0:
        t = args[0]
    else:
        t = "No data to plot"
    plt.text(0.4, 0.5, t)


if __name__ == "__main__":
    print(inspect.getfullargspec(plt.plot))
