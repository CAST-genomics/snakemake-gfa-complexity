import sys
import numpy as np
import pandas as pd
from pathlib import Path
import statsmodels.api as sm
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype

"""
Usage:
python manhattan_range.py [trait]
"""

AXIS_FONTSIZE = 6
TITLE_FONTSIZE = 5
AXIS_LABELSIZE = 2.5
LABEL_FONTSIZE = 4
TICK_FONTSIZE = 4
POINT_SIZE = 0.75

def main():
    trait = sys.argv[1]
    linear = trait+ ".assoc.sorted"
    bed_input = trait+ ".sig.filtered.bed"

    output = trait+".png"

    ids=tuple()
    orange_ids=tuple()
    orange_hids=False
    red_chr_ids=False
    label=True
    small=False
    titles=tuple()
    alpha=None

    """
    Create a manhattan plot from the results of a PLINK2 GWAS
    """
    plink_cols = {
    "#CHROM": "chromosome",
    "POS": "pos",
    "ID": "id",
    "P": "pval",
    }
    keep_cols = ["#CHROM", "POS", "ID", "P"]
    # create the plot

    plot_width = 30
    fig, ax = plt.subplots(
        1, 1, sharex=True, sharey=True, constrained_layout=True,
        figsize=(plot_width, plot_width*1/4)
    )


    # parse the .linear files
    max_pval = -1
    red_ids = ids

    idx = 0

    # first, handle cases where there may be more than one .linear file
    cur_ax = ax
    if cur_ax is tuple:
        cur_ax = cur_ax[idx]



    df = pd.read_csv(
        linear,
        sep="\t",
        header=0,
        usecols=keep_cols,
    ).rename(columns=plink_cols)
    df = df.sort_values("pos")
    pos_range = max(df["pos"]) - min(df["pos"])
    label_distance = pos_range/17
    # replace NaN with inf
    df["pval"] = df["pval"].fillna(np.inf)
    df["pval"] = df["pval"].fillna(np.inf)


    df["pval"] = df["pval"].replace(0, min(df[df["pval"] > 0]["pval"])/50)
    df['-log10(p)'] = -np.log10(df["pval"])
    # replace -infinity values with 0
    df['-log10(p)'].replace([-np.inf], 0, inplace=True)



    df['chromosome_int'] = df['chromosome'].str.extract('(\d+)')
    df['chromosome_cat'] = df['chromosome_int'].astype('category')

    #df.chromosome_cat = df.chromosome_cat.astype(
    #    CategoricalDtype(sorted(map(int, df.chromosome_cat.dtype.categories)), ordered=True)
    #)

    df = df.sort_values('chromosome_cat')

    # My input
    chr_lengths = [ 
    [1,247_249_719],
    [2,242_951_149],
    [3,199_501_827],
    [4,191_273_063],
    [5,180_857_866],
    [6,170_899_992],
    [7,158_821_424],
    [8,146_274_826],
    [9,140_273_252],
    [10,135_374_737],
    [11,134_452_384],
    [12,132_349_534],
    [13,114_142_980],
    [14,106_368_585],
    [15,100_338_915],
    [16,88_827_254],
    [17,78_774_742],
    [18,76_117_153],
    [19,63_811_651],
    [20,62_435_964],
    [21,46_944_323],
    [22,49_691_432] ]

    for i in range(0, len(chr_lengths)):
        chr_lengths[i][1] += chr_lengths[i - 1][1]

    df['chromosome_int'] = df['chromosome'].str.extract('(\d+)').astype(int)
    df["x_pos"] = [chr_lengths[x][1] for x in (df['chromosome_int']-1).replace(-1, 0)] + df["pos"]

    # TODO Bed stuff
    bed_cols = ["chromosome", "start", "end", "id", "Z"]


    bed_df = pd.read_csv(
        bed_input,
        sep="\t",
        names = bed_cols
    )
    bed_df['chromosome_int'] = bed_df['chromosome'].str.extract('(\d+)').astype(int)
    bed_df['chromosome_len'] = [chr_lengths[x][1] for x in (bed_df['chromosome_int']-1).replace(-1, 0)]
    bed_df["x_start"] = bed_df['chromosome_len'] + bed_df["start"]
    bed_df["x_end"] = bed_df['chromosome_len'] + bed_df["end"]

    for index, row in bed_df.iterrows():
        cur_ax.axvspan(row['x_start'], row['x_end'], color='red', alpha=0.3)

    # TODO Changed
    df[~df["id"].isin(red_ids + orange_ids)].plot(
        kind='scatter', x='x_pos', y='-log10(p)', ax=cur_ax,
        color='blue', alpha=0.5, s=10
    )

    # set title and perform cleanup/secondary tasks
    cur_ax.set_title(trait)
    cur_ax.set(xlabel=None, ylabel=None)
    max_val = cur_ax.get_ylim()[1]
    if max_pval < max_val:
        max_pval = max_val

    # if we haven't already flipped alpha to the log scale, do it now
    if alpha is not None and alpha > 0 and alpha < 0.5:
        alpha = -np.log10(alpha)

    # set the y-axis limit so that both axes have the same limit
    cur_ax.set_ylim(top=max_pval)
    if alpha is not None:
        cur_ax.axhline(y=alpha, color='red')

    fig.supxlabel('Chromosomal Position')
    fig.supylabel('$-log_{10} P-value$')

    plt.savefig(output)

if __name__ == "__main__":
    main()