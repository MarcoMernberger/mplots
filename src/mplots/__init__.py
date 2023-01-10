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
<<<<<<< HEAD
from .scatter import volcano_plot, volcano_calc, plot_correlation
=======
from .scatter import volcano_plot, volcano_calc, generate_dr_plot, plot_correlation
from .functions import *
from .bars import *
>>>>>>> 2d1fedbfbe7a84d951a0318e1f06909feb5b14ae
