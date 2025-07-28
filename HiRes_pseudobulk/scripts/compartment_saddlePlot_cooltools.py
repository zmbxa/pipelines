# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess

# Import python package for working with cooler files and tools for analysis
import cooler
import cooltools.lib.plotting
from packaging import version
import cooltools
import bioframe

import argparse

### get args
parser = argparse.ArgumentParser(description="This script is for calculate PC1 using cooltools")
parser.add_argument("-f","--file",help="REQUIRED, cooler files that you want to calculate on, like ./HiRes_pseudobulk.mcool",dest="file")
parser.add_argument("-r","--resolution",help="REQUIRED, cooler file resolution you want to use, default is 50000",dest="res",default=50000)
parser.add_argument("-t","--threads",help="optional, number of processor you want to use, default is 10",dest="threads",default=10)
args = parser.parse_args()

resolution=args.res
cooler_file=args.file
nProcessor=args.threads


from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import Normalize, LogNorm
from cytoolz import merge
def saddleplot_label(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    """
    Generate a saddle plot.
    Parameters

    ----------
    track : pd.DataFrame

        See cooltools.digitize() for details.
    saddledata : 2D array-like

        Saddle matrix produced by `make_saddle`. It will include 2 flanking

        rows/columns for outlier signal values, thus the shape should be

        `(n+2, n+2)`.
    cmap : str or matplotlib colormap

        Colormap to use for plotting the saddle heatmap

    scale : str

        Color scaling to use for plotting the saddle heatmap: log or linear

    vmin, vmax : float

        Value limits for coloring the saddle heatmap

    color : matplotlib color value

        Face color for margin bar plots

    fig : matplotlib Figure, optional

        Specified figure to plot on. A new figure is created if none is

        provided.
    fig_kws : dict, optional

        Passed on to `plt.Figure()`
    heatmap_kws : dict, optional

        Passed on to `ax.imshow()`
    margin_kws : dict, optional

        Passed on to `ax.bar()` and `ax.barh()`
    cbar_kws : dict, optional

        Passed on to `plt.colorbar()`
    subplot_spec : GridSpec object

        Specify a subregion of a figure to using a GridSpec.
    Returns

    -------
    Dictionary of axes objects.
    """

    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    from matplotlib.colors import Normalize, LogNorm

    from matplotlib import ticker

    import matplotlib.pyplot as plt

    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])

        def __call__(self, x, pos=None):
            if x not in [vmin, 1, vmax]:
                return ""
            else:
                return "{x:g}".format(x=x)

    track_value_col = track.columns[3]
    track_values = track[track_value_col].values

    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange

    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    if qrange is not None:
        lo, hi = qrange

        binedges = np.linspace(lo, hi, n_bins + 1)

    # Barplot of mean values and saddledata are flanked by outlier bins

    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata

    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    # Layout

    if subplot_spec is not None:
        GridSpec = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    grid = {}
    gs = GridSpec(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    # Figure

    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    # Heatmap

    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only linear and log color scaling is supported")

    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    # Margins

    margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})
    # left margin hist

    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])

    plt.barh(
        binedges, height=1/len(binedges), width=groupmean, align="edge", **margin_kws

    )

    plt.xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr

    plt.ylim(hi, lo)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    # top margin hist

    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])

    plt.bar(
        binedges, width=1/len(binedges), height=groupmean, align="edge", **margin_kws

    )

    plt.xlim(lo, hi)
    # plt.ylim(plt.ylim())  # correct

    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

    # Colorbar

    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
    if scale == "linear" and vmin is not None and vmax is not None:
        grid["ax_cbar"] = cb = plt.colorbar(img, **cbar_kws)
        decimal = 10

        nsegments = 5

        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal

        cb.set_ticks(cd_ticks)
    else:
        cb = plt.colorbar(img, format=MinOneMaxFormatter(), cax=grid["ax_cbar"], **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    # extra settings

    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid['ax_heatmap'].grid(False)
    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    # Calculate and annotate the mean values in the top-left and bottom-right 10% corners

    n = C.shape[0]
    top_left_mean = np.mean(C[:n//10, :n//10])
    bottom_right_mean = np.mean(C[-n//10:, -n//10:])

    ax.text(0.05, 0.9, f'{top_left_mean:.2f}', transform=ax.transAxes, color='w')
    ax.text(0.85, 0.05, f'{bottom_right_mean:.2f}', transform=ax.transAxes, color='w')

    return grid


mm10_genome = bioframe.load_fasta('/storage/zhangyanxiaoLab/share/fasta/mm10.fa');

base_name = os.path.basename(cooler_file).split('.mc')[0]
print(base_name)
# Load the cooler file
clr = cooler.Cooler(f'{cooler_file}::resolutions/{resolution}')
# Get the bins
bins = clr.bins()[:]
# Compute GC content
gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], mm10_genome)
# Define the view DataFrame
view_df = pd.DataFrame({'chrom': clr.chromnames,'start': 0,
'end': clr.chromsizes.values,'name': clr.chromnames})
# Compute eigenvectors
cis_eigs = cooltools.eigs_cis(clr,gc_cov,view_df=view_df,n_eigs=3,)
    
# Extract the eigenvector track
eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]
    
cvd = cooltools.expected_cis(clr=clr,nproc=nProcessor,view_df=view_df)
# calculate saddle matrix
interaction_sum, interaction_count =  cooltools.saddle(
  clr,cvd,eigenvector_track,'cis',
  n_bins=100,qrange=(0.025,0.975),view_df=view_df)
saddleplot_label(eigenvector_track,interaction_sum/interaction_count,100,title=base_name,qrange=(0.025,0.975),cbar_kws={'label':'average observed/expected contact frequency'});
plt.savefig(f'saddleplot_{base_name}.png', dpi=300)
plt.close()
    
    # save saddle matrix
pd.DataFrame(interaction_sum/interaction_count).to_csv(f'{base_name}_saddleMatrix_50k.tsv', index=False, sep='\t')
    
