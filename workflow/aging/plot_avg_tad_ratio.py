from bioframe import read_table
import cooler
import pandas as pd
from coolpuppy import coolpup
from os.path import basename
import seaborn as sns
from matplotlib import pyplot as plt
sns.set_theme(font_scale=2)

def plot_average_tad_ratio(clr_nom_path, clr_den_path,
                           tad_nom_path, tad_den_path,
                           exp_nom_path, exp_den_path,
                           output, view_path, 
                           resolution, min_size):
    resolution = int(resolution)
    clr_nom = cooler.Cooler(clr_nom_path + f'::/resolutions/{resolution}')
    clr_den = cooler.Cooler(clr_den_path + f'::/resolutions/{resolution}')
    tad_nom = read_table(tad_nom_path, schema='bed4')
    tad_den = read_table(tad_den_path, schema='bed4')
    exp_nom = pd.read_csv(exp_nom_path, sep='\t')
    exp_den = pd.read_csv(exp_den_path, sep='\t')
    view = read_table(view_path, schema='bed4')
    name_nom = basename(clr_nom_path).split('.')[0]
    name_den = basename(clr_den_path).split('.')[0]
    #create pileups
    cc_nom = coolpup.CoordCreator(tad_nom, resolution=resolution, 
                                features_format='bed', local=True, rescale_flank=0.5)
    cc_den = coolpup.CoordCreator(tad_den, resolution=resolution, 
                                features_format='bed', local=True, rescale_flank=0.5)
    pu_nom = coolpup.PileUpper(clr_nom, cc_nom, 
                                expected=exp_nom, 
                                rescale=True,
                                view_df=view, ignore_diags=0).pileupsWithControl()
    pu_den = coolpup.PileUpper(clr_den, cc_den, 
                                expected=exp_den, 
                                rescale=True,
                                view_df=view, ignore_diags=0).pileupsWithControl()
    #plot avg tad ratio
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    fig.suptitle(f'Average TADs change: min size={min_size}', 
                    fontsize=20, y=1.05)
    divv = pu_nom['data'][0] / pu_den['data'][0]
    img = ax.imshow(divv, cmap='RdBu_r', vmax=1.20, vmin=0.8, interpolation='bicubic')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title(f'{name_nom}/{name_den}', fontsize=18)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.83, 0.15, 0.025, 0.7])
    fig.colorbar(img, cax=cbar_ax)

    plt.savefig(output, format='pdf', dpi=300, bbox_inches='tight')


#main
plot_average_tad_ratio(clr_nom_path=snakemake.input.sample_NOM, 
                       clr_den_path=snakemake.input.sample_DEN,
                        tad_nom_path=snakemake.input.tad_NOM, 
                        tad_den_path=snakemake.input.tad_DEN,
                        exp_nom_path=snakemake.input.exp_NOM, 
                        exp_den_path=snakemake.input.exp_DEN,
                        output=snakemake.output[0], 
                        view_path=snakemake.input.view, 
                        resolution=snakemake.params.resolution, 
                        min_size=snakemake.params.win
                        )