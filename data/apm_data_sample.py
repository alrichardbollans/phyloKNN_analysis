import os

import pandas as pd
from tqdm import tqdm

from data.helper_functions import input_data_path
from data.umapping import reduction_factor, unsupervised_umap_wrapper, embed_distances_with_umap
from phyloAutoEncoder import autoencode_pairwise_distances
number_of_apm_folds=10

def main():
    for tag in tqdm(range(1, number_of_apm_folds + 1)):
        dir_path = os.path.join(input_data_path, 'my_apm_data', 'ultrametric','APM', str(tag))
        unsupervised_umap_wrapper(dir_path)

        distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)
        unsup_model, unsupervised_df = autoencode_pairwise_distances(distances, reduction_factor, dir_path, plot=True)
        unsupervised_df.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny.csv'))

def with_full_tree():
    db_path = os.environ.get('KEWDROPBOXPATH')
    # Data from APM Traits
    binary_data = pd.read_csv(os.path.join(db_path, 'ApmTraits', 'apm_activity', 'outputs', 'compiled_extraction_apm_data.csv'))
    number_of_output_features = int(len(binary_data)*reduction_factor)

    repo_path = os.environ.get('KEWSCRATCHPATH')
    dir_path = os.path.join(repo_path, 'gentianales_trees', 'WCVP_12', 'Uphy', 'outputs',
                                        'Species')
    tree_path = os.path.join(dir_path, 'species_distances.csv')
    distances = pd.read_csv(tree_path, index_col=0)
    reduction_fraction = number_of_output_features/len(distances.columns)

    umap_embedding = embed_distances_with_umap(distances, reduction_fraction)
    umap_embedding.to_csv(os.path.join(dir_path, 'umap_unsupervised_embedding_full_tree.csv'))

    unsup_model, unsupervised_df = autoencode_pairwise_distances(distances, reduction_fraction, dir_path, plot=True)
    unsupervised_df.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny_full_tree.csv'))


if __name__ == '__main__':
    main()
    with_full_tree()