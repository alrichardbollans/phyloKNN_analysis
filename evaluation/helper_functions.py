import os
import pathlib

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from data.helper_functions import number_of_simulation_iterations

bin_model_names = ['corHMM',  'phylnn_raw', 'phylnn_fill_means',
                   'logit_eigenvecs', 'logit_umap', 'logit_umap_supervised', 'logit_autoencoded',
                   'logit_autoenc_supervised',
                   'xgb_eigenvecs', 'xgb_umap', 'xgb_umap_supervised', 'xgb_autoencoded', 'xgb_autoenc_supervised']
bin_model_names.remove('phylnn_raw')
cont_model_names = ['phylopars',  'phylnn_raw', 'phylnn_fill_means',
                    'linear_eigenvecs', 'linear_umap', 'linear_autoencoded',
                    'xgb_eigenvecs', 'xgb_umap', 'xgb_autoencoded']
cont_model_names.remove('phylnn_raw')
rename_models_and_ev_models = {'phylnn_fill_means': 'phylokNN', 'logit_eigenvecs': 'Eigenvec (L)',
                               'logit_umap': 'UMAP (L)',
                               'logit_umap_supervised': 'UMAP* (L)', 'logit_autoencoded': 'Autoenc (L)',
                               'logit_autoenc_supervised': 'Autoenc* (L)',
                               'xgb_eigenvecs': 'Eigenvec (XGB)', 'xgb_umap': 'UMAP (XGB)',
                               'xgb_umap_supervised': 'UMAP* (XGB)', 'xgb_autoencoded': 'Autoenc (XGB)',
                               'xgb_autoenc_supervised': 'Autoenc* (XGB)',
                               'linear_eigenvecs': 'Eigenvec (L)', 'linear_umap': 'UMAP (L)',
                               'linear_umap_supervised': 'UMAP* (L)', 'linear_autoencoded': 'Autoenc (L)',
                               'linear_autoenc_supervised': 'Autoenc* (L)',
                               'logit_eigenvecs_full_tree': 'Eigenvec (Full tree) (L)',
                               'logit_umap_full_tree': 'UMAP (Full tree) (L)', 'logit_autoencoded_full_tree': 'Autoenc (Full tree) (L)',
                               'xgb_eigenvecs_full_tree': 'Eigenvec (Full tree) (XGB)', 'xgb_umap_full_tree': 'UMAP (Full tree) (XGB)',
                               'xgb_autoencoded_full_tree': 'Autoenc (Full tree) (XGB)'
                               }

binary_model_order = ['corHMM', 'phylokNN', 'Eigenvec (L)', 'Eigenvec (XGB)', 'UMAP (L)', 'UMAP* (L)',
                      'UMAP (XGB)', 'UMAP* (XGB)',
                      'Autoenc (L)', 'Autoenc* (L)', 'Autoenc (XGB)', 'Autoenc* (XGB)']

continuous_model_order = ['phylopars', 'phylokNN', 'Eigenvec (L)', 'Eigenvec (XGB)', 'UMAP (L)',
                          'UMAP (XGB)',
                          'Autoenc (L)', 'Autoenc (XGB)']


def check_scales():
    # This is just a sanity check to check the scaling
    EB_values = []
    BMT_values = []
    standard_values = []

    for tag in range(1, number_of_simulation_iterations + 1):
        tag = str(tag)
        eb_val = pd.read_csv(
            os.path.join('..', 'data', 'non_standard_simulations', 'EB', 'continuous', tag, 'ground_truth.csv'))[
            'trait_EB_scaled'].tolist()
        EB_values += eb_val

        bmt_val = pd.read_csv(
            os.path.join('..', 'data', 'non_standard_simulations', 'BMT', 'continuous', tag, 'ground_truth.csv'))[
            'trait_BM_trend_scaled'].tolist()
        BMT_values += bmt_val

        standard_val = pd.read_csv(os.path.join('..', 'data', 'simulations', 'continuous', tag, 'ground_truth.csv'))[
            'F1.1/1'].tolist()
        standard_values += standard_val
    print('EB values:', np.mean(EB_values))
    print('BMT values:', np.mean(BMT_values))
    print('Standard values:', np.mean(standard_values))
    plot_df = pd.DataFrame({'standard_values': standard_values, 'EB_values': EB_values, 'BMT_values': BMT_values})
    sns.boxplot(plot_df)
    plt.show()


def output_df(df, bin_or_cont, out_dir, group: str = 'Ev Model'):
    for group_member in df[group].unique():
        ev_df = df[df[group] == group_member]
        filename = "".join(x for x in group_member if x.isalnum())
        out_path = os.path.join(out_dir, bin_or_cont)
        pathlib.Path(out_path).mkdir(exist_ok=True, parents=True)
        ev_df.describe(include='all').to_csv(os.path.join(out_path, f'{filename}_mean_results.csv'))