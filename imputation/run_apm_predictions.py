import os

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, make_scorer
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
from xgboost import XGBClassifier

from data.apm_data_sample import number_of_apm_folds
from imputation.helper_functions import get_prediction_data_paths, phylnn_predict
from imputation.run_encodings_predictions import get_umap_data, add_y_to_data, get_eigenvectors, get_autoencoded_data, logit_init_kwargs, \
    logit_grid_search_params, fit_and_output, xgb_clf_init_kwargs, xgb_clf_grid_search_params, get_semi_supervised_umap_data, \
    get_semi_supervised_autoencoded_data
from phylokNN import nan_safe_metric_wrapper


def run_predictions():
    for iteration in tqdm(range(1, number_of_apm_folds + 1)):
        m = 'mcar'
        bin_or_cont = 'binary'

        real_or_sim = 'my_apm_data'
        average_precision_score_nan_safe = nan_safe_metric_wrapper(average_precision_score)
        _scorer = make_scorer(average_precision_score_nan_safe, greater_is_better=True, response_method='predict_proba')
        case = 'ultrametric'
        ev_model = 'APM'
        phylnn_predict(case, ev_model, iteration, m, _scorer)

        umap_X = get_umap_data(case, ev_model, iteration)
        umap_df, umap_encoding_vars, umap_target_name = add_y_to_data(umap_X, case, ev_model, iteration, m)

        eigen_X = get_eigenvectors(case, ev_model, iteration)
        eigen_df, eigen_encoding_vars, eigen_target_name = add_y_to_data(eigen_X, case, ev_model, iteration, m)

        autoenc_X = get_autoencoded_data(case, ev_model, iteration)
        autoenc_df, autoenc_encoding_vars, autoenc_target_name = add_y_to_data(autoenc_X, case, ev_model, iteration, m)
        out_dir = get_prediction_data_paths(case, ev_model, iteration, m)

        if bin_or_cont == 'binary':
            clf_instance = LogisticRegression(**logit_init_kwargs)
            fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_umap', umap_df, umap_encoding_vars, umap_target_name,
                           bin_or_cont, scorer='average_precision')

            clf_instance = LogisticRegression(**logit_init_kwargs)
            fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_eigenvecs', eigen_df, eigen_encoding_vars,
                           eigen_target_name, bin_or_cont, scorer='average_precision')

            clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
            fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_umap', umap_df, umap_encoding_vars, umap_target_name,
                           bin_or_cont, scorer='average_precision')

            clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
            fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_eigenvecs', eigen_df, eigen_encoding_vars,
                           eigen_target_name, bin_or_cont, scorer='average_precision')
            # ### autoencoder
            clf_instance = LogisticRegression(**logit_init_kwargs)
            fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_autoencoded', autoenc_df, autoenc_encoding_vars,
                           eigen_target_name, bin_or_cont, scorer='average_precision')

            clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
            fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_autoencoded', autoenc_df, autoenc_encoding_vars,
                           umap_target_name,
                           bin_or_cont, scorer='average_precision')
            # # ### Semisupervised cases not included. The general simulation analysis shows they dont provide improvement
            # # ### Semisupervised umap
            semi_supervised_umap_df, semi_sup_umap_encoding_vars, semi_sup_umap_target_name = get_semi_supervised_umap_data(case,
                                                                                                                            ev_model,
                                                                                                                            iteration, m)
            clf_instance = LogisticRegression(**logit_init_kwargs)
            fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_umap_supervised', semi_supervised_umap_df,
                           semi_sup_umap_encoding_vars,
                           semi_sup_umap_target_name, bin_or_cont)
            clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
            fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_umap_supervised', semi_supervised_umap_df,
                           semi_sup_umap_encoding_vars,
                           semi_sup_umap_target_name, bin_or_cont)

            ## Semisupervised autoenc
            semi_supervised_autoenc_df, semi_sup_autoenc_encoding_vars, semi_sup_autoenc_target_name = get_semi_supervised_autoencoded_data(
                case, ev_model,
                iteration, m)

            clf_instance = LogisticRegression(**logit_init_kwargs)
            fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_autoenc_supervised', semi_supervised_autoenc_df,
                           semi_sup_autoenc_encoding_vars,
                           semi_sup_autoenc_target_name, bin_or_cont)
            clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
            fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_autoenc_supervised', semi_supervised_autoenc_df,
                           semi_sup_autoenc_encoding_vars,
                           semi_sup_autoenc_target_name, bin_or_cont)

def run_predictions_for_encodings_with_full_tree():
    repo_path = os.environ.get('KEWSCRATCHPATH')
    dir_path = os.path.join(repo_path, 'gentianales_trees', 'WCVP_12', 'Uphy', 'outputs',
                            'Species')
    m = 'mcar'
    bin_or_cont = 'binary'
    case = 'ultrametric'
    ev_model = 'APM'
    ## Umap data
    umap_X = pd.read_csv(os.path.join(dir_path, 'umap_unsupervised_embedding_full_tree.csv'), index_col=0)
    ## Scale the data
    umap_X = pd.DataFrame(StandardScaler().fit_transform(umap_X), index=umap_X.index)


    ## Eigenvecs
    eigen_X = pd.read_csv(os.path.join(dir_path, 'all_eigenvectors_full_tree.csv'), index_col=0)

    broken_stick_params = pd.read_csv(os.path.join(dir_path, 'broken_stick_parameters.csv'), index_col=0)
    num_cols_to_use = broken_stick_params['broken_stick_number'].iloc[0]
    eigen_X = eigen_X.iloc[:, : num_cols_to_use]
    ## Scale the data
    eigen_X = pd.DataFrame(StandardScaler().fit_transform(eigen_X), index=eigen_X.index)

    ## Autoencoded data
    autoenc_X = pd.read_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny_full_tree.csv'), index_col=0)
    ## Scale the data
    autoenc_X = pd.DataFrame(StandardScaler().fit_transform(autoenc_X), index=autoenc_X.index)
    for iteration in tqdm(range(1, number_of_apm_folds + 1)):

        umap_df, umap_encoding_vars, umap_target_name = add_y_to_data(umap_X, case, ev_model, iteration, m)

        eigen_df, eigen_encoding_vars, eigen_target_name = add_y_to_data(eigen_X, case, ev_model, iteration, m)

        autoenc_df, autoenc_encoding_vars, autoenc_target_name = add_y_to_data(autoenc_X, case, ev_model, iteration, m)
        out_dir = get_prediction_data_paths(case, ev_model, iteration, m)

        raise NotImplementedError('Just use with best selected method')

        if bin_or_cont == 'binary':
            clf_instance = LogisticRegression(**logit_init_kwargs)
            fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_umap_full_tree', umap_df, umap_encoding_vars, umap_target_name,
                           bin_or_cont, scorer='average_precision')

            clf_instance = LogisticRegression(**logit_init_kwargs)
            fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_eigenvecs_full_tree', eigen_df, eigen_encoding_vars,
                           eigen_target_name, bin_or_cont, scorer='average_precision')

            clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
            fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_umap_full_tree', umap_df, umap_encoding_vars, umap_target_name,
                           bin_or_cont, scorer='average_precision')

            clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
            fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_eigenvecs_full_tree', eigen_df, eigen_encoding_vars,
                           eigen_target_name, bin_or_cont, scorer='average_precision')
            # # ### autoencoder
            clf_instance = LogisticRegression(**logit_init_kwargs)
            fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_autoencoded_full_tree', autoenc_df, autoenc_encoding_vars,
                           umap_target_name, bin_or_cont, scorer='average_precision')

            clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
            fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_autoencoded_full_tree', autoenc_df, autoenc_encoding_vars,
                           umap_target_name,
                           bin_or_cont, scorer='average_precision')


if __name__ == '__main__':
    run_predictions()
    # run_predictions_for_encodings_with_full_tree()
