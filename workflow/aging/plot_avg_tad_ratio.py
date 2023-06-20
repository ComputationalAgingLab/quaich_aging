from os.path import basename
import seaborn as sns
from matplotlib import pyplot as plt
from coolpuppy.lib.io import load_pileup_df
sns.set_theme(font_scale=2)

def plot_average_tad_ratio(pu_nom_path, pu_den_path,
                           name_nom, name_den, min_size,
                           output, 
                           ):
    #load pileups
    pu_nom = load_pileup_df(pu_nom_path)
    pu_den = load_pileup_df(pu_den_path)
    divv = pu_nom['data'][0] / pu_den['data'][0]
    #plot avg tad ratio
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    fig.suptitle(f'Average TADs change: min size={min_size}', 
                    fontsize=20, y=1.05)
    img = ax.imshow(divv, cmap='RdBu_r', vmax=1.20, vmin=0.8, interpolation='bicubic')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title(f'{name_nom}/{name_den}', fontsize=18)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.83, 0.15, 0.025, 0.7])
    fig.colorbar(img, cax=cbar_ax)

    for out in output:
    	plt.savefig(out, dpi=300, bbox_inches='tight')


#main
plot_average_tad_ratio(pu_nom_path=snakemake.input.pu_NOM, 
                       pu_den_path=snakemake.input.pu_DEN,
                       name_nom=snakemake.wildcards['sampleNOM'], 
                       name_den=snakemake.wildcards['sampleDEN'], 
                       min_size=snakemake.params.win,
                       output=snakemake.output, 
                        )
