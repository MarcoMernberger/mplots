#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import matplotlib as mpl
import matplotlib.pyplot as plt
import inspect
import warnings
import pandas as pd
from mplots.newplot import compact_wrapper, style_wrapper, add_function_wrapper
from pandas import DataFrame


@pytest.fixture
def df_test():
    return pd.DataFrame({"A": [1, 2, 3], "B": [4, 3, 4]})


@compact_wrapper
def myplot(df, **kwargs):
    plt.plot(df.A, df.B, **kwargs)


@style_wrapper
@add_function_wrapper
def myplot2(df, **kwargs):
    plt.bar(df.A, df.B, **kwargs)


@style_wrapper
@add_function_wrapper
def myplot3(df, **kwargs):
    plt.imshow(df, **kwargs)


def test_compact_wrapper(df_test):
    f = myplot(
        df_test,
        xlabel=4,
        label="myline",
        styleparams={"lines.marker": "o"},
        title={"label": "test", "fontsize": 20},
        ylabel="y",
        linewidth=10,
        ls="--",
    )
    assert isinstance(f, plt.Figure)


def test_style_and_function_wrapper(df_test):
    f = myplot2(
        df_test,
        xlabel=4,
        label="myline",
        styleparams={"lines.marker": "o"},
        title={"label": "test", "fontsize": 20},
        ylabel="y",
        linewidth=10,
        ls="--",
        styles=["Solarize_Light2"],
    )
    assert isinstance(f, plt.Figure)


def test_heatmap(df_test):
    f = myplot3(df_test, figsize=(10, 10))
    assert isinstance(f, plt.Figure)
