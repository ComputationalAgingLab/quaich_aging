import cooler
import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def compute_cis_trans(clr_path, output_path=None, dropdiag=2, resolution=100000, stable=True, normalized=False):
    """
    This function returns a dictionary where keys are names of Hi-C matrices
    and values are matrices of cis-trans contact sums.
    
    hiclist :: [str] - list of paths to .mcool files
    dropdiag :: int - number of diags to drop (main diagonal counts from 1)
    resolution :: int - desired resolution for Hi-C matrix
    stable :: bool - whether to add 1 to contact sums matrix (for numerical 
                     stability at logarithm computation)
    normalized :: bool - whether to normalize matrix with a total sum of contacts
    """

    clr = cooler.Cooler(clr_path + f'::/resolutions/{resolution}')
    chromnames = clr.chromnames
    sizes = np.array(clr.chromsizes.values) / int(resolution) #compute normalization matrix
    norm = np.outer(sizes, sizes)

    pix = clr.pixels()[:]
    pix_no_diag = pix[pix.bin2_id >= pix.bin1_id + dropdiag] #delete m first diagonals
    bins = clr.bins()[:].drop(['start', 'end'], axis=1)
    cis_trans = np.asarray(
                pix_no_diag.merge(bins, left_on='bin1_id', right_index=True).\
                            merge(bins, left_on='bin2_id', right_index=True).\
                            groupby(['chrom_x', 'chrom_y']).sum().\
                            unstack(fill_value=0)['count']) #convert to square matrix

    cis_trans = cis_trans + cis_trans.T - np.diag(np.diag(cis_trans))
    if normalized:
        fullsum = np.nansum(cis_trans)
        if fullsum == 0.:
            raise ValueError('%s has zero contact sum.' % clr_path)
        cis_trans = cis_trans / fullsum * 2 * 1_000_000 #scale by constant
    if stable:
        cis_trans += 1
    result = cis_trans / norm
    result = pd.DataFrame(result, index=chromnames, columns=chromnames)
    if output_path is not None:
        result.to_csv(output_path, sep='\t')
    else:
        return result
    

#main
compute_cis_trans(  
                f'{snakemake.input[0]}', 
                f'{snakemake.output[0]}', 
                resolution=f'{snakemake.params.resolution}',
                dropdiag=2, stable=True, normalized=False)