#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""pipegraph.py: Contains methods and classes for pypipoegraph compatibility."""

from pypipegraph import Job
from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any
import pypipegraph as ppg

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class MPlotJob(ppg.PlotJob):
    def __init__(self, calculation_function, plot_function):
        pass

