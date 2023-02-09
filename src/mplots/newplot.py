import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import inspect
import warnings


class CustomizationManager:
    """A class to set the default customization values"""

    def __init__(self, **kwargs):
        raise NotImplementedError


class Plot:
    def __init__(self, plot_function):
        raise NotImplementedError
        self.plot_function = plot_function
        self._customization_manager = CustomizationManager()


STYLE_SHORTCUTS = {
    "ls": "lines.linestyle",
    "figsize": "figure.figsize",
    "linestyle": "lines.linestyle",
    "linewidth": "lines.linewidth",
    "marker": "lines.marker",
    "markeredgecolor": "lines.markeredgecolor",
    "markeredgewidth": "lines.markeredgewidth",
    "markerfacecolor": "lines.markerfacecolor",
    "markersize": "lines.markersize",
}


def call_additional_func(func, params):
    if isinstance(params, dict):
        func(**params)
    elif isinstance(params, list):
        func(*params)
    else:
        func(params)


def add_calls(**kwargs):
    remaining = {}
    for key in kwargs:
        if hasattr(plt, key):
            func = getattr(plt, key)
            if callable(getattr(plt, key)):
                call_additional_func(func, kwargs[key])
        else:
            remaining[key] = kwargs[key]
    return remaining


def create_style_dict(**kwargs):
    local_style = {}
    remaining = {}
    for key in kwargs:
        if key == "styleparams":
            local_style.update(kwargs[key])
        elif key in mpl.rcParams.keys():
            local_style[key] = kwargs[key]
        elif key in STYLE_SHORTCUTS:
            local_style[STYLE_SHORTCUTS[key]] = kwargs[key]
        else:
            remaining[key] = kwargs[key]
    return local_style, remaining


def signature_name_is_pyplot_function(key):
    return hasattr(plt, key)


def name_warning(key):
    warnings.warn(
        f"The parameter {key} clashes with pyplot functions. Adding this python function via keyword argument is disabled.",
        SyntaxWarning,
    )
    raise NotImplementedError


def check_parameter_name_clash(func):
    for key in inspect.signature(func).parameters:
        if signature_name_is_pyplot_function(key):
            name_warning(key)


def compact_wrapper(func):
    def wrapper(*args, **kwargs):
        check_parameter_name_clash(func)
        local_style, remaining_kwargs = create_style_dict(**kwargs)
        with plt.style.context(["Solarize_Light2", local_style]):
            f = plt.figure()
            remaining_kwargs = add_calls(**remaining_kwargs)
            func(*args, **remaining_kwargs)
            return f

    return wrapper


def add_function_wrapper(func):
    def wrapper(*args, **kwargs):
        check_parameter_name_clash(func)
        f = plt.figure()
        remaining_kwargs = add_calls(**kwargs)
        func(*args, **remaining_kwargs)
        return f

    return wrapper


def style_wrapper(func):
    def wrapper(*args, **kwargs):
        styles_to_use = kwargs.pop("styles", [])
        local_style, remaining_kwargs = create_style_dict(**kwargs)
        styles_to_use.append(local_style)
        with plt.style.context(styles_to_use):
            return func(*args, **remaining_kwargs)

    return wrapper


"""
wanted behaviour:

def box(df, parameter, **kwargs):
    plt.boxplot()

box = BoxPlot


"""


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


def print_mpl_rc_params():
    for key in mpl.rcParams.keys():
        print(key.split(".")[-1], key)
