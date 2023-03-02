import seaborn as sns
from os.path import basename
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def plot_eigenvector_correlation(vecs, vals, output, 
                                 resolution, fix_chr, fix_row='best'):
    """
    vecs - list of eigvec paths; 
    vals - list of eigval paths;
    fix_chr - choose chromosome name or `all`
    fix_row - choose eigvec name or `best` (choose best according to eigval)
    """
    keys = [basename(v).split('_')[0] for v in vecs]
    eigvec = {k:pd.read_csv(v, sep='\t') for k, v in zip(keys, vecs)}
    eigval = {k:pd.read_csv(v, sep='\t', index_col=0) for k, v in zip(keys, vals)}
    chroms = eigval[keys[0]].index.tolist()

    #combine eigenvectors from different sources
    comp_df = pd.DataFrame()
    for i, key in enumerate(keys):
        tmp = eigvec[key]
        if fix_chr == 'all':
            comp_vector = []
            for c in chroms:
                if fix_row == 'best':
                    main_eigvec = 'E' + str(np.argmax(eigval[key].loc[c].filter(regex='eigval*'))+1)
                else:
                    main_eigvec = fix_row
                comp_vector.append(tmp[tmp.chrom==c][main_eigvec].to_numpy())
            comp_df[key] = np.concatenate(comp_vector)
        elif (fix_chr != 'all') and (fix_row == 'best'):
            main_eigvec = 'E' + str(np.argmax(eigval[key].loc[fix_chr].filter(regex='eigval*'))+1)
            comp_vector = tmp[tmp.chrom==fix_chr][main_eigvec].to_numpy()
            comp_df[key] = comp_vector
        else:        
            comp_df[key] =  tmp[tmp['chrom']==fix_chr][fix_row]

    ax = sns.clustermap(comp_df.corr(), annot=True, figsize=(4,4))
    ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xmajorticklabels(), fontsize = 14)
    ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize = 14, rotation=0)
    ax.fig.suptitle(f'Compartment vectors cluster map, resolution={resolution}', x=0.5, y=1.1)
    plt.savefig(output, dpi=300, bbox_inches='tight', format='pdf')


#main
plot_eigenvector_correlation(vecs=snakemake.input.eigvecs, 
                             vals=snakemake.input.eigvals, 
                             output=snakemake.output[0], 
                             resolution=snakemake.params.resolution, 
                             fix_chr=snakemake.params.chrom, 
                             fix_row=snakemake.params.eig,)