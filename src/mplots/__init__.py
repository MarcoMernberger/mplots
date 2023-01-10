# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound


from .plots import *
from .jobs import MPPlotJob
from .heatmaps import heatmap
from .scatter import volcano_plot, volcano_calc, generate_dr_plot, plot_correlation
from .functions import *
from .bars import *
