import os.path

import pandas as pd
import umap
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

from data.helper_functions import number_of_simulation_iterations, input_data_path, reduction_factor


def embed_distances_with_umap(distances: pd.DataFrame, reduction_fraction=reduction_factor):
    scaled_penguin_data = StandardScaler().fit_transform(distances)

    reducer = umap.UMAP(n_components=int(len(distances.columns) * reduction_fraction))

    embedding = reducer.fit_transform(scaled_penguin_data)
    umap_embedding = pd.DataFrame(embedding, index=distances.index)
    return umap_embedding


def unsupervised_umap_wrapper(dir_path: str, reduction_fraction=reduction_factor):
    distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)

    umap_embedding = embed_distances_with_umap(distances, reduction_fraction)
    umap_embedding.to_csv(os.path.join(dir_path, 'umap_unsupervised_embedding.csv'))


def main():
    for tag in tqdm(range(1, number_of_simulation_iterations + 1)):
        for case in ['ultrametric', 'with_extinct']:
            dir_path = os.path.join(input_data_path, 'simulations', case,'standard', str(tag))
            unsupervised_umap_wrapper(dir_path)

            dir_path = os.path.join(input_data_path, 'simulations', case, 'BiSSE', str(tag))
            unsupervised_umap_wrapper(dir_path)

            dir_path = os.path.join(input_data_path, 'simulations', case, 'HiSSE', str(tag))
            unsupervised_umap_wrapper(dir_path)

        dir_path = os.path.join(input_data_path, 'real_data', 'ultrametric','Clonality', str(tag))
        unsupervised_umap_wrapper(dir_path)

        dir_path = os.path.join(input_data_path, 'real_data', 'ultrametric','Seed Mass', str(tag))
        unsupervised_umap_wrapper(dir_path)


if __name__ == '__main__':
    main()
