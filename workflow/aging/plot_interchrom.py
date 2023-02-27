import pandas as pd
from matplotlib import pyplot as plt
from os.path import basename
import numpy as np

def plot_interchroms(sample_nominator, sample_denominator, output):
    """
    The function plots interchromosomal contact ratio. 
    Note that Input Hi-C maps should be downsampled.
    """
    i1 = pd.read_csv(sample_nominator, sep='\t', index_col=0)
    i2 = pd.read_csv(sample_denominator, sep='\t', index_col=0)
    name1 = basename(sample_nominator).split('_')[0]
    name2 = basename(sample_denominator).split('_')[0]
    chroms = i1.columns.tolist()

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), sharey=True, sharex=True)
    fig.suptitle('Interchromosomal contact ratio', fontsize=16, y=0.92, x=0.43)
    fig.tight_layout(pad=1.2)
    map1 = i1.values #sample_nominator
    map2 = i2.values #sample_denominator
    m = map1 / map2
    #np.fill_diagonal(m, np.nan)
    ax.set_title(name1 + '/' + name2 + '; average=%.3f' % np.nanmean(m), fontsize=16)
    im = ax.imshow(m, cmap='PuOr_r', vmax=1.2, vmin=0.8, extent=[-1,1,1,-1])

    tm_grid = np.arange(-1+1/len(chroms), 1, 2/len(chroms))
    ax.set_xticks(tm_grid)
    ax.set_yticks(tm_grid)
    ax.set_xticklabels(chroms, rotation=45, fontsize=11)
    ax.set_yticklabels(chroms, rotation=0, fontsize=13)
    ax.grid(False)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.1, 0.02, 0.75])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    plt.savefig(output, dpi=300, bbox_inches='tight', format='pdf')



#main
plot_interchroms(snakemake.input.sample_NOM, 
                 snakemake.input.sample_DEN, 
                 snakemake.output[0])