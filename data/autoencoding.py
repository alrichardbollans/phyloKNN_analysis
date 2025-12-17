import os

import pandas as pd
from tqdm import tqdm

from data.helper_functions import number_of_simulation_iterations, input_data_path
from data.umapping import reduction_factor
from phyloAutoEncoder import autoencode_pairwise_distances


def main():
    for tag in tqdm(range(1, number_of_simulation_iterations + 1)):
        for case in ['ultrametric', 'with_extinct']:

            dir_path = os.path.join(input_data_path, 'simulations', case,'standard', str(tag))
            distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)
            unsup_model, unsupervised = autoencode_pairwise_distances(distances, reduction_factor, dir_path,plot=True)
            unsupervised.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny.csv'))

            dir_path = os.path.join(input_data_path, 'simulations', case, 'BiSSE', str(tag))
            distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)
            unsup_model, unsupervised = autoencode_pairwise_distances(distances, reduction_factor, dir_path,plot=True)
            unsupervised.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny.csv'))

            dir_path = os.path.join(input_data_path, 'simulations', case, 'HiSSE', str(tag))
            distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)
            unsup_model, unsupervised = autoencode_pairwise_distances(distances, reduction_factor, dir_path,plot=True)
            unsupervised.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny.csv'))

        dir_path = os.path.join(input_data_path, 'real_data', 'ultrametric','Clonality', str(tag))
        distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)
        unsup_model, unsupervised = autoencode_pairwise_distances(distances, reduction_factor, dir_path,plot=True)
        unsupervised.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny.csv'))

        dir_path = os.path.join(input_data_path, 'real_data', 'ultrametric','Seed Mass', str(tag))
        distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)
        unsup_model, unsupervised = autoencode_pairwise_distances(distances, reduction_factor, dir_path,plot=True)
        unsupervised.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny.csv'))


if __name__ == '__main__':
    main()
