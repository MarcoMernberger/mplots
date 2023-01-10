#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""bars.py: Contains ...."""

import matplotlib
matplotlib.use("agg")
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def calc_ora(df, topx: int):
    df = df[df["Benjamini"] <= 0.05].sort_values("Benjamini", ascending=True)
    df = df[: min(len(df), topx)]
    return df


def plot_ora(df, **kwargs):
    figsize = kwargs.get("figsize", (10, 10))
    title = kwargs.get("title", "Enriched gene sets")
    fig = plt.figure(figsize=figsize)
    xlabel = kwargs.get("xlabel", r"$-log_{10}$ (corrected p-value)")
    ylabel = kwargs.get("ylabel", "enriched gene sets")
    fontsize = kwargs.get("fontsize", 10)
    fontsize_title = kwargs.get("fontsize_title", fontsize)
    if len(df):
        df["corr p_value"] = -1 * np.log10(df["Benjamini"].values)
        df = df.sort_values("corr p_value", ascending=True)
        df["Set"] = pd.Categorical(df["Set"].values, list(set(df["Set"].values)))
        y = range(len(df))
        plt.barh(
            y,
            df["corr p_value"].values,
            height=0.8,
            tick_label=df["Set"].values,
        )
        lims = plt.ylim()
        plt.gca().axvline([-1 * math.log(0.05, 10)], ymin=lims[0], ymax=lims[1], color="k")
        plt.xlabel(xlabel, fontsize=fontsize)
        plt.ylabel(ylabel, fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.title(title, fontsize=fontsize_title)
    else:
        plt.title("No significant entries")
    plt.tight_layout()
    return fig
