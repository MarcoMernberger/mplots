#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""heatmaps.py: Contains some basic reusable heatmap functions."""
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import Colormap
from cycler import cycler
import matplotlib.gridspec as grid
from pandas import DataFrame
from typing import Union, Tuple
from .newplot import style_wrapper, add_function_wrapper

matplotlib.use("Agg")


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


@style_wrapper
@add_function_wrapper
def plot_correlation_map(df, **kwargs):
    df_corr = df.corr(method="pearson")
    plt.imshow(df_corr, **kwargs)
    plt.colorbar(shrink=0.2)


@style_wrapper
@add_function_wrapper
def hmap(df):
    plt.imshow(df, cmap="seismic", aspect="equal")
    plt.xticks(np.arange(df.shape[1]), labels=df.columns, rotation=90)
    plt.colorbar(shrink=0.5)


@style_wrapper
def multihmap(*dfs, title="MultiHeatmap"):
    f = plt.figure(constrained_layout=True)
    xgrid = sum([d.shape[1] for d in dfs])
    gs = f.add_gridspec(11, xgrid)
    dfs = list(dfs)
    ncol = len(dfs)
    if ncol == 1:
        axes = [axes]
    x = 0
    for ii, df in enumerate(dfs):
        newx = x + df.shape[1]
        ax = f.add_subplot(gs[:, x:newx])
        if ii > 0:
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.tick_params(left=False, labelleft=False)
        else:
            plt.yticks(np.arange(df.shape[0]), labels=df.index)
        plt.sca(ax)
        plt.imshow(df, aspect="auto", cmap="seismic")
        plt.xticks(np.arange(df.shape[1]), labels=df.columns, rotation=90)
        x = newx
    plt.suptitle(title)
    plt.tight_layout()
    return f


def calc_correlation_group_heatmap(df, method="pearson"):
    return df.corr(method="pearson")


# def plot_correlation_group_heatmap(df_corr, **kwargs):
#    title, legend = get_args(title, legend)
#    print(title, legend)
#    plt.title(title=tt)
#    plt.imshow(df_corr, kwargs)
#    plt.xlabel(df_corr.columns)
#    plt.ylabel(df_corr.columns)


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


# from mdataframe
def generate_heatmap_figure(
    df,
    title=None,
    label_function=None,
    display_linkage_column=None,
    display_linkage_row=None,
    legend_location=None,
    show_column_label=True,
    show_row_label=True,
    **params
):
    """
    Plots a heatmap.
    @param df dataframe with values to plot. Plots the whole df, labels are df.columns and df.index.
    @display_linkage_column name of column clustering that has and linkage attribute to be plotted
    @display_linkage_row name of row clustering that has and linkage attribute to be plotted
    @legend_location where the color scale is shown
    @show_row_label plot row labels
    @show_column_label plot column labels
    """
    # set some parameters
    rows_per_inch = params.get("rows_per_inch", 6)
    columns_per_inch = params.get("columns_per_inch", 2)
    dpi = params.get("dpi", 300)
    fontsizex = params.get("fontsizex", 8)
    fontsizey = params.get("fontsizey", math.floor(71 / rows_per_inch))
    inch_for_dendro_top = params.get("inch_for_dendro_top", 1)
    inch_for_dendro_left = params.get("inch_for_dendro_left", 1)
    colormap = params.get("colormap", "seismic")
    vmax = params.get("vmax", df.max().max())
    vmin = params.get("vmin", df.min().min())
    aspect = params.get("aspect", "auto")
    dodge = params.get("dodge", 0)
    if "row_label_inches" in params:
        row_label_inches = params["row_label_inches"]
    else:
        max_row_label = np.max([len(str(label)) for label in df.index.values])
        row_label_inches = math.ceil(fontsizey * (1 / 72) * max_row_label / dpi)
    if "column_label_inches" in params:
        column_label_inches = params["column_label_inches"]
    else:
        max_column_label = np.max([len(str(label)) for label in df.columns.values])
        column_label_inches = math.ceil(fontsizex * 2.6 * max_column_label / dpi)  # 1
    fontsize_title = params.get("fontsize_title", 12)
    # set the plot dimensions
    len_x = len(df.columns)
    len_y = len(df.index)
    inches_top = dodge
    inches_bottom = dodge + column_label_inches
    inches_left = dodge
    inches_right = dodge + row_label_inches  # 1
    start_dendro_left = dodge
    start_dendro_top = dodge

    # check and set color scale position and adjust plot dimensions accordingly
    allowed_legend_locations = [
        "r",
        "b",
        "t",
        "l",
        "br",
        "bl",
        "rb",
        "lb",
        "tr",
        "tl",
        "rt",
        "lt",
    ]
    if legend_location is not None and legend_location not in allowed_legend_locations:
        raise ValueError("Legend location must be one of {}.".format(allowed_legend_locations))
    if not show_row_label:
        inches_right = dodge
        if legend_location in ["br", "rb", "tr", "rt"]:
            inches_right = dodge + 1
    if not show_column_label:
        inches_bottom = dodge
        if legend_location in ["br", "rb", "bl", "lb"]:
            inches_bottom += 1
    if legend_location == "r":
        inches_right += 1
    elif legend_location == "l":
        inches_left += 1
        start_dendro_left = 1
    elif legend_location == "b":
        inches_bottom += 1
    elif legend_location == "t":
        start_dendro_top = 1
        inches_top += 1
    if display_linkage_row is not None:
        inches_left += inch_for_dendro_left
    elif legend_location in ["bl", "lb", "tl", "lt"]:
        inches_left += 1
    if display_linkage_column is not None:
        inches_top += inch_for_dendro_top
    elif legend_location in ["tl", "tr", "rt", "lt"]:
        inches_top += 1

    # create the figure and grid
    data_inches_x = len_x / columns_per_inch
    data_inches_y = len_y / rows_per_inch
    grid_x = len_x + columns_per_inch * inches_left + columns_per_inch * inches_right
    grid_y = len_y + rows_per_inch * inches_top + rows_per_inch * inches_bottom
    fig_x = data_inches_x + inches_right + inches_left
    fig_y = data_inches_y + inches_bottom + inches_top
    check_y = fig_y * dpi
    check_x = fig_x * dpi
    if check_x > 60000 or check_y > 60000:
        raise ValueError(
            "The number of pixels is too large: {}x{}, both dimensions should be below 60k. Try decreasing the dpi or using multipage plot.".format(
                check_x, check_y
            )
        )
    f = plt.figure(figsize=(fig_x, fig_y), dpi=dpi, frameon=True, edgecolor="k", linewidth=2)
    gridspec = grid.GridSpec(grid_y, grid_x)
    gridspec.update(wspace=0, hspace=0)

    # create the axes
    axes = []
    if display_linkage_column is not None:
        ax_dendro_column = plt.subplot(
            gridspec[
                start_dendro_top * rows_per_inch : start_dendro_top * rows_per_inch
                + rows_per_inch * inch_for_dendro_top,
                columns_per_inch * inches_left : len_x + columns_per_inch * inches_left,
            ]
        )
        axes.append(ax_dendro_column)
    if display_linkage_row is not None:
        ax_dendro_row = plt.subplot(
            gridspec[
                rows_per_inch * inches_top : len_y + rows_per_inch * inches_top,
                start_dendro_left * columns_per_inch : start_dendro_left * columns_per_inch
                + columns_per_inch * inch_for_dendro_left,
            ]
        )
        axes.append(ax_dendro_row)
    if show_column_label:
        ax_labels_column = plt.subplot(
            gridspec[
                rows_per_inch * (inches_top)
                + len_y : rows_per_inch * (inches_top + inches_bottom)
                + len_y,
                columns_per_inch * inches_left : len_x + columns_per_inch * inches_left,
            ]
        )
        axes.append(ax_labels_column)
    if show_row_label:
        ax_labels_row = plt.subplot(
            gridspec[
                rows_per_inch * inches_top : len_y + rows_per_inch * inches_top,
                columns_per_inch * inches_left
                + len_x : columns_per_inch * (inches_left + 1)
                + len_x,
            ]
        )
        axes.append(ax_labels_row)
    ax_matrix = plt.subplot(
        gridspec[
            rows_per_inch * inches_top : len_y + rows_per_inch * inches_top,
            columns_per_inch * inches_left : len_x + columns_per_inch * inches_left,
        ]
    )

    # format main axis
    ax_matrix.set_xticks(np.arange(len_x))
    ax_matrix.set_yticks(np.arange(len_y))
    ax_matrix.tick_params(
        bottom=True,
        right=True,
        left=False,
        top=False,
        labelbottom=show_column_label,
        labeltop=False,
        labelright=show_row_label,
        labelleft=False,
    )
    norm = matplotlib.colors.Normalize(vmin, vmax)

    # add the labels
    rows = list(df.index.values)
    columns = list(df.columns.values)
    if label_function is not None:
        rows = [label_function(label) for label in rows]
        columns = [label_function(label) for label in columns]
    ax_matrix.set_yticklabels(rows, fontsize=fontsizey)
    ax_matrix.set_xticklabels(columns, rotation=75, ha="right", fontsize=fontsizex)

    # actually plot the heatmap
    im = ax_matrix.imshow(df, cmap=colormap, norm=norm, aspect=aspect)

    # add the color scale axis
    if legend_location is not None:
        cbar_tickparams = {
            "left": False,
            "right": False,
            "labelright": False,
            "labelleft": False,
            "top": False,
            "bottom": False,
            "labeltop": False,
            "labelbottom": False,
            "labelsize": "small",
        }
        orientation = "vertical"
        if legend_location == "br" or legend_location == "rb":
            ax_color = plt.subplot(
                gridspec[
                    grid_y - rows_per_inch * inches_bottom : grid_y,
                    grid_x - columns_per_inch * inches_right : grid_x,
                ]
            )
            bounds = list(ax_color.get_position().bounds)
            w = bounds[2] / 5
            h = bounds[3] * 0.95
            bounds[0] = bounds[0] + 0.7 * w
            bounds[1] = bounds[1] + 0.025 * bounds[3]
            bounds[2] = w
            bounds[3] = h
            cbar_tickparams["right"] = True
            cbar_tickparams["labelright"] = True
        elif legend_location == "bl" or legend_location == "lb":
            ax_color = plt.subplot(
                gridspec[
                    grid_y - rows_per_inch * inches_bottom : grid_y,
                    0:columns_per_inch,
                ]
            )
            bounds = list(ax_color.get_position().bounds)
            w = bounds[2] / 5
            h = bounds[3] * 0.95
            bounds[0] = bounds[0] + bounds[2] - 1.5 * w
            bounds[1] = bounds[1] + 0.025 * bounds[3]
            bounds[2] = w
            bounds[3] = h
            cbar_tickparams["left"] = True
            cbar_tickparams["labelleft"] = True

        elif legend_location == "tr" or legend_location == "rt":
            ax_color = plt.subplot(
                gridspec[
                    0:rows_per_inch,
                    grid_x - columns_per_inch * inches_right : grid_x,
                ]
            )
            bounds = list(ax_color.get_position().bounds)
            h = 0.9 * bounds[3]
            w = bounds[2] / 5
            bounds[0] = bounds[0] + 0.7 * w
            bounds[1] = bounds[1] + 0.06 * bounds[3]
            bounds[2] = w
            bounds[3] = h
            cbar_tickparams["right"] = True
            cbar_tickparams["labelright"] = True
        elif legend_location == "tl" or legend_location == "lt":
            ax_color = plt.subplot(gridspec[0:rows_per_inch, 0:columns_per_inch])
            bounds = list(ax_color.get_position().bounds)
            h = 0.9 * bounds[3]
            w = bounds[2] / 5
            bounds[0] = bounds[0] + bounds[2] - 1.5 * w
            bounds[1] = bounds[1] + 0.06 * bounds[3]
            bounds[2] = w
            bounds[3] = h
            cbar_tickparams["left"] = True
            cbar_tickparams["labelleft"] = True
        elif legend_location == "r":
            ax_color = plt.subplot(
                gridspec[
                    rows_per_inch * inches_top : len_y + rows_per_inch * inches_top,
                    grid_x - columns_per_inch : grid_x,
                ]
            )
            bounds = list(ax_color.get_position().bounds)
            w = bounds[2] / 5
            bounds[2] = w
            bounds[0] += 0.5 * w
            cbar_tickparams["right"] = True
            cbar_tickparams["labelright"] = True
        elif legend_location == "l":
            ax_color = plt.subplot(
                gridspec[
                    rows_per_inch * inches_top : len_y + rows_per_inch * inches_top,
                    0:columns_per_inch,
                ]
            )
            bounds = list(ax_color.get_position().bounds)
            w = bounds[2] / 5
            bounds[0] = bounds[0] + bounds[2] - 1.5 * w
            bounds[2] = w
            cbar_tickparams["left"] = True
            cbar_tickparams["labelleft"] = True
        elif legend_location == "t":
            ax_color = plt.subplot(
                gridspec[
                    0:rows_per_inch,
                    columns_per_inch * inches_left : len_x + columns_per_inch * inches_left,
                ]
            )
            bounds = list(ax_color.get_position().bounds)
            orientation = "horizontal"
            h = bounds[3] / rows_per_inch
            bounds[1] = bounds[1] + 0.5 * h
            bounds[3] = h
            cbar_tickparams["top"] = True
            cbar_tickparams["labeltop"] = True

        elif legend_location == "b":
            ax_color = plt.subplot(
                gridspec[
                    grid_y - rows_per_inch : grid_y,
                    columns_per_inch * inches_left : len_x + columns_per_inch * inches_left,
                ]
            )
            bounds = list(ax_color.get_position().bounds)
            orientation = "horizontal"
            h = bounds[3] / rows_per_inch
            bounds[1] = bounds[1] + bounds[3] - 1.5 * h
            bounds[3] = h
            cbar_tickparams["bottom"] = True
            cbar_tickparams["labelbottom"] = True

        ax_color.tick_params(
            bottom=False,
            right=False,
            left=False,
            top=False,
            labelbottom=False,
            labeltop=False,
            labelright=False,
            labelleft=False,
        )
        axes.append(ax_color)
        # since the color scale should be slightly smaller than the reserved grid axis,
        # we need to add another smaller axis, which in total makes:
        # ... seven evil axes
        ax_color2 = f.add_axes(bounds)
        cbar = plt.colorbar(im, cax=ax_color2, orientation=orientation)
        cbar.ax.tick_params(**cbar_tickparams)

    # remove spines and ticks
    for ax in axes:
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)
        ax.tick_params(bottom=False, right=False, left=False, top=False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

    # set the title
    if title is not None:
        title_len = len(title)
        title_size = math.ceil(fontsize_title * 2.6 * title_len / dpi)
        if title_size > fig_x:
            fontsize_title = max(1, math.floor((fig_x * dpi) / (2.6 * len(title))))
        plt.suptitle(title, fontsize=fontsize_title)
    plt.close()
    return f


def generate_heatmap_simple_figure(
    df, title=None, show_column_label=True, show_row_label=False, label_function=None, **params
):
    """
    Plots a heatmap. No fancy shit.
    @df dataframe with values to plot. Plots the whole df, labels are df.columns no index.
    @show_column_label plot column labels
    """
    fontsize_title = params.get("fontsize_title", 12)
    fontsize_columns = params.get("fontsize_column", 12)
    fontsize_rows = params.get("fontsize_rows", 12)
    dpi = params.get("dpi", 100)
    fig_x = params.get("fig_x", df.shape[1])
    fig_y = params.get("fig_y", max(10, len(df) * 4 / dpi))
    guess_shrink = 0.5
    if fig_y > 30:
        guess_shrink = 0.1
    shrink = params.get("shrink", guess_shrink)
    f = plt.figure(figsize=(fig_x, fig_y), dpi=dpi, frameon=True, edgecolor="k", linewidth=2)
    colormap = params.get("colormap", "seismic")
    vmax = params.get("vmax", df.abs().max().max())
    vmin = -vmax  # to ensure zero midpoint
    aspect = params.get("aspect", "auto")
    norm = matplotlib.colors.Normalize(vmin, vmax)
    im = plt.gca().imshow(df, cmap=colormap, norm=norm, aspect=aspect)
    plt.colorbar(im, shrink=shrink)
    #  add the labels
    if show_column_label:
        columns = df.columns.values
        if label_function is not None:
            columns = [label_function(label) for label in df.columns.values]
        plt.gca().set_xticks(range(len(df.columns)))
        plt.gca().set_xticklabels(columns, rotation=75, ha="right", fontsize=fontsize_columns)
    else:
        plt.gca().set_xticklabels([])
    if show_row_label:
        rows = df.index
        if label_function is not None:
            rows = [label_function(label) for label in df.index.values]
        plt.gca().set_yticks(range(len(df.index)))
        plt.gca().set_yticklabels(rows, ha="right", fontsize=fontsize_rows)
    else:
        plt.gca().set_yticklabels([])
    #  set the title
    if title is not None:
        plt.title(title, fontsize=fontsize_title)
    plt.tight_layout(pad=1, w_pad=1, h_pad=1)
    return f


def plot_empty(outfile, msg="Empty DataFrame"):
    fig = plt.figure()
    plt.text(0.5, 0.5, msg, ha="left", va="center")
    fig.savefig(outfile)


@style_wrapper
@add_function_wrapper
def plot_correlation(df, correlation="pearson"):
    df = df.corr(correlation)
    plt.imshow(df)
    t, l = np.arange(df.shape[1]), df.columns
    plt.yticks(t, l)
    plt.xticks(t, l, rotation=90)
    plt.colorbar(shrink=0.5)
