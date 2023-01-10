#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
functions.py: Contains default functions and helper for plot generation and
customization.
"""
import matplotlib
from cycler import cycler, Cycler
from typing import Dict, List, Any

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def default_marker(label: str) -> str:
    return "o"


def default_cycler() -> cycler:
    custom_cycler = cycler(
        color=[
            "b",
            "g",
            "r",
            "c",
            "k",
            "m",
            "y",
            "grey",
            "darkblue",
            "darkgreen",
            "darkred",
            "darkcyan",
            "darkviolet",
            "gold",
            "slategrey",
        ]
    ) + cycler(marker=["o", "v", "^", "*", "s", "<", ">", "+", "o", "v", "^", "*", "s", "<", ">"])
    return custom_cycler


def custom_cycler(color_marker: Dict[str, List[Any]]) -> Cycler:
    custom_cycler = cycler(**color_marker)
    return matplotlib.rcsetup.cycler(custom_cycler)
