import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pypipegraph as ppg
from pathlib import Path


def volcanoplot(
    genes,
    fdr_column,
    logFC_column,
    significance_threshold=.01,
    fc_threhold=1,
    outfile=None,
    dependencies=[],
    name=None,
    overlapping_genes=None,
    with_label=False,
    **kwargs,
):
    """
    overlapping_genes: genes object that should get another different color. First color: regulated genes; second color: regulated genes that are also in overlapping genes; third color: overlapping genes NOT regulated
    """
    name = kwargs.get("name", f"{genes.name}")
    figsize = kwargs.get("figsize", (10, 10))
    title = name if name is not None else genes.name
    if outfile is None:
        outfile = genes.result_dir / f"{title}_volcano.svg"
    if isinstance(outfile, str):
        outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    deps = [genes.load()] + dependencies
    if overlapping_genes is not None:
        deps.append(overlapping_genes.load())

    def calc():
        df = genes.df.copy()
        df = df.rename({foldchange_column_name: "log2 fold change", benjamini_column: "-log10(p-value)"})
        return df

    def plot(df):
        figure = plt.figure(figsize=figsize)
        df["group"] = ["grey"]*len(df)
        df["group"][(df[logFC_column].values >= 1) & (df[fdr_column] <= .05)] = "red"
        df["group"][(df[logFC_column].values <= -1) & (df[fdr_column] <= .05)] = "blue"
        colors = ["grey", "red", "blue"]
        labels = {"grey": "non-sign.", "red": "up", "blue": "down"}
        fig = plt.figure(figsize=(10, 10))
        for group in colors:
            df_sub = df[df["group"] == group]
            plt.plot(
                df_sub[logFC_column].values, 
                -np.log10(df_sub[fdr_column].values), 
                ls="", 
                marker="o", 
                color=group, 
                label=labels[group]
                )
        plt.axhline(-np.log10(.05), color="lightgrey")
        plt.axvline(-1, color="lightgrey")
        plt.axvline(1, color="lightgrey")
        plt.ylabel(r"-log($p_{corrected}$)")
        plt.xlabel("log2FC")
        plt.legend()
        plt.title(title)
        return figure
    if ppg.inside_ppg():
        return ppg.MatPlotJob(
            outfile, calc, plot, render_args={"width": 12, "height": 9}
        ).depends_on(deps)
    else:
        