import pandas as pd
from matplotlib import ticker, pyplot as plt
import numpy as np
from os.path import basename

def plot_scaling_ratio(pairs, output, resolution, max_dist, window, column='count.sum', ignore_diags=2):
    fig, ax = plt.subplots(1, 1, figsize=(10,5))
    ax.set_title('Scaling ratio plot')
    ax.set_xlabel(f'Separation distance, b')
    ax.set_ylabel('Ratio')
    ax.axhline(y=1.0, ls='--', color='grey', lw=1.5)
    ax.xaxis.set_major_locator(ticker.LogLocator(base = 10.0, 
                                                 subs = np.arange(1.0, 10.0) * 0.1, 
                                                 numticks = 10))
    for pair in pairs:
        sample_nominator = pair[0]
        sample_denominator = pair[1]
        s1 = pd.read_csv(sample_nominator, sep='\t')
        s2 = pd.read_csv(sample_denominator, sep='\t')
        name1 = basename(sample_nominator).split('_')[0]
        name2 = basename(sample_denominator).split('_')[0]
        #prepare data
        agg_scaling = pd.DataFrame()
        s1 = s1.groupby('dist')[column].mean()
        s2 = s2.groupby('dist')[column].mean()
        s1[:ignore_diags] = np.nan
        s2[:ignore_diags] = np.nan
        agg_scaling['nom'] = s1
        agg_scaling['den'] = s2
        agg_scaling.index = s1.groupby('dist').mean().index * int(resolution)
        #plot data
        tmp = agg_scaling[agg_scaling.index <= int(max_dist)]
        smooth = tmp['nom'].rolling(window = int(window)).mean() / tmp['den'].rolling(window = int(window)).mean()
        ax.semilogx(smooth, label = f'{name1}/{name2}', lw=1.8)
    ax.legend()
    plt.savefig(output, dpi=300, bbox_inches='tight', format='pdf')

#main
pairs = list(zip(snakemake.input.nominators, snakemake.input.denominators))

plot_scaling_ratio(
    pairs=pairs,
    output=snakemake.output[0],
    resolution=snakemake.params.resolution, 
    max_dist=snakemake.params.max_dist, 
    window=snakemake.params.window,
)