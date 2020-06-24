#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""pipegraph.py: Contains methods and classes for pypipoegraph compatibility."""

from pypipegraph import Job, PlotJob, ppg_exceptions, FileGeneratingJob, ParameterInvariant, FunctionInvariant
from pypipegraph.job import was_inited_before
from pathlib import Path
from typing import Callable, List, Union, Any
import pypipegraph.util as util
import pandas as pd
from matplotlib.figure import Figure
import pypipegraph as ppg
import os
import pickle


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class MPPlotJob(PlotJob):
    """Calculate some data for plotting, cache it in cache/output_filename, and plot from there.
    creates two jobs, a plot_job (this one) and a cache_job (FileGeneratingJob, in self.cache_job),

    Just like a regulat plot job, but it uses matplotlib.
    """

    def __init__(
        self,
        output_filename: Path,
        calc_function: Callable,
        plot_function: Callable,
        skip_table: bool = False,
        skip_caching: bool = False,
        calc_args: List[Any] = None,
        plot_args: List[Any] = None,
    ) -> None:
        """
        A plot job that separates data calculation and actual plotting in order
        to reduce runtime. Just like a regular pypipegraph.PlotJob, but using
        matplotlib instead of ggplot.

        As this uses matplotlib as a plotting library, the DataFrame that is
        calculated, follows a different structure. Instead of having one data
        point per row, this assumes a table structure as in Genes for example.
        Plot parameters that alter the appearance of the figure can be supplied
        as keyword arguments and are used to check, wether the figure needs
        recalculation.

        Parameters
        ----------
        output_filename : Path
            Path to output file.
        calc_function : Callable
            Callable that returns the DataFrame used for plotting.
        plot_function : Callable
            Callable that returns a matplotlib.pyplot.figure and gets the
            DataFrame from calc_function as sole parameter.
        skip_table : bool, optional
            If True, no DataFrame will be written to file, by default False.
        skip_caching : bool, optional
            If True, caching will be skipped, by default False.
        calc_args : List[Any], optional
            List of arguments on which the DataFrame calculation depends, by
            default None.
        plot_args : List[Any], optional
            List of arguments on which the DataFrame calculation depends, by
            default None.

        Raises
        ------
        ValueError
            If a MPPlotjob with the same id has been called before.
        ValueError
            If an unknown file extension is specified.
        ppg_exceptions.JobContractError
            If the calc_function did not return a DataFrame or dict.
        ppg_exceptions.JobContractError
            If the plot_function did not return a matplotlib.pyplot.figure.
        """
        if was_inited_before(self, MPPlotJob):
            if (
                skip_caching != self.skip_caching
                or skip_table != self.skip_table
                or calc_args != self.calc_args
                or plot_args != self.plot_args
            ):
                raise ValueError(
                    "MPPlotJob(%s) called twice with different parameters" % self.job_id
                )
            return
        if not (
            output_filename.suffix == ".png"
            or output_filename.suffix == ".pdf"
            or output_filename.suffix == ".svg"
        ):
            raise ValueError(
                f"Don't know how to create this file {output_filename}, must end on .png or .pdf or .svg")

        self.output_filename = output_filename
        self.filenames = [output_filename]
        self._check_for_filename_collisions()
        self.table_filename = self.output_filename.parent / (self.output_filename.name + ".tsv")
        self.calc_function = calc_function
        self.plot_function = plot_function
        self.skip_caching = skip_caching
        self.skip_table = skip_table
        self._fiddle = None
        self.calc_args = calc_args if calc_args is not None else []
        self.plot_args = plot_args if plot_args is not None else []
        self.plot_args += self.calc_args
        if not self.skip_caching:
            if output_filename.is_absolute():
                self.cache_filename = Path(util.global_pipegraph.cache_folder) / output_filename.relative_to("/")
            else:
                self.cache_filename = Path(util.global_pipegraph.cache_folder) / output_filename
            self.cache_filename = self.cache_filename.with_suffix(".calc")
            self.cache_filename.parent.mkdir(exist_ok=True, parents=True)

            def run_calc():
                df = calc_function()
                if not isinstance(df, pd.DataFrame):
                    do_raise = True
                    if isinstance(df, dict):  # might be a list dfs...
                        do_raise = False
                        for x in df.values():
                            if not isinstance(x, pd.DataFrame):
                                do_raise = True
                                break
                    if do_raise:
                        raise ppg_exceptions.JobContractError(
                            f"{output_filename}.calc_function did not return a DataFrame (or dict of such), was {type(df)}."
                        )
                try:
                    os.makedirs(os.path.dirname(self.cache_filename))
                except OSError:
                    pass
                of = open(self.cache_filename, "wb")
                pickle.dump(df, of, pickle.HIGHEST_PROTOCOL)
                of.close()

        def run_plot():
            df = self.get_data()
            fig = plot_function(df)
            if not isinstance(fig, Figure):
                raise ppg_exceptions.JobContractError(
                    f"%{output_filename}.plot_function did not return a matplotlib.pyplot.figure object, was {type(fig)}."
                )
            if self._fiddle is not None and callable(self._fiddle):
                self._fiddle(fig)
            fig.savefig(output_filename)

        FileGeneratingJob.__init__(self, output_filename, run_plot)
        if plot_args is not None:
            Job.depends_on(
                self, ParameterInvariant(str(self.output_filename) + "_plot_params", self.plot_args)
            )
        if not self.skip_caching:
            cache_job = FileGeneratingJob(self.cache_filename, run_calc)
            if calc_args is not None:
                cache_job.depends_on(
                    ParameterInvariant(str(self.output_filename) + "_calc_params", self.calc_args)
                )
            Job.depends_on(self, cache_job)
            self.cache_job = cache_job
        if not skip_table:

            def dump_table():
                import pandas as pd

                df = self.get_data()
                if isinstance(df, pd.DataFrame):
                    df.to_csv(self.table_filename, sep="\t")
                else:
                    with open(self.table_filename, "w") as op:
                        for key, dframe in df.items():
                            op.write("#%s\n" % key)
                            dframe.to_csv(op, sep="\t")

            table_gen_job = FileGeneratingJob(self.table_filename, dump_table)
            if not self.skip_caching:
                table_gen_job.depends_on(cache_job)
            self.table_job = table_gen_job
        else:
            self.table_job = None

    def add_another_plot(self, output_filename: Path, plot_function: Callable, plot_args: List[Any] = None) -> Job:
        """
        Add another plot job that runs on the same data as the original one.

        No recalculation is done.

        Parameters
        ----------
        output_filename : Path
            Path to generated file.
        plot_function : Callable
            Callable that returns the new figure.
        plot_args : List[Any], optional
            List of arguments on which the DataFrame calculation depends, by
            default None.

        Returns
        -------
        Job
            The FileGeneratingJob that produces the file.
        """

        def run_plot():
            df = self.get_data()
            fig = plot_function(df)
            if not isinstance(fig, Figure):
                raise ppg_exceptions.JobContractError(
                    f"%{output_filename}.plot_function did not return a matplotlib.pyplot.figure object, was {type(fig)}."
                )
            fig.savefig(output_filename)

        job = FileGeneratingJob(output_filename, run_plot)
        if plot_args is not None:
            job.depends_on(
                ParameterInvariant(str(output_filename) + "_params", plot_args)
            )
        job.depends_on(FunctionInvariant(str(output_filename) + "_plotfunc", plot_function))
        job.depends_on(self.cache_job)
        return job

    def plot_depends_on(self, job_or_list):
        if isinstance(job_or_list, Job):
            self.prerequisites.add(job_or_list)
        elif hasattr(job_or_list, "__iter__"):
            for job in job_or_list:
                self.prerequisites.add(job)
        else:
            raise ValueError(f"Unexpected joblist: {type(job_or_list)}.")