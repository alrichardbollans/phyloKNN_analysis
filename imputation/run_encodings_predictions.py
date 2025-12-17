import os
import sys


import pandas as pd
import umap
from phyloAutoEncoder import autoencode_pairwise_distances
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import brier_score_loss, mean_absolute_error
from sklearn.model_selection import GridSearchCV, StratifiedKFold, KFold
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
from xgboost import XGBClassifier, XGBRegressor

sys.path.append('../..')
from data.helper_functions import number_of_simulation_iterations, reduction_factor, simulation_types
from imputation.helper_functions import missingness_types, get_input_data_paths, \
    get_prediction_data_paths, n_split_for_nested_cv, get_bin_or_cont_from_ev_model


xgb_clf_init_kwargs = {'eval_metric': brier_score_loss}
xgb_clf_grid_search_params = {'max_depth': [1, 3, 6, 10], 'learning_rate': [0.01, 0.1, 0.3], 'subsample': [0.5, 0.8, 1], 'gamma': [0, 0.1, 1],
                              'max_delta_step': [0, 1]}

logit_init_kwargs = {}
logit_grid_search_params = {'C': [0.001, 0.01, 0.1, 1, 10],
                            'class_weight': ['balanced', None, {0: 0.1, 1: 0.9}, {0: 0.15, 1: 0.85}, {0: 0.3, 1: 0.7}]}

xgb_rgr_init_kwargs = {'eval_metric': mean_absolute_error}
xgb_rgr_grid_search_params = xgb_clf_grid_search_params.copy()

linear_init_kwargs = {}
linear_grid_search_params = {}


def fit_and_output(clf_instance, grid_search_param_grid, out_dir, model_name, full_df, encoding_vars, target_name, bin_or_cont, scorer=None):
    assert target_name not in encoding_vars
    X_test = full_df[full_df[target_name].isna()]
    X_train = full_df[~full_df[target_name].isna()]

    values = X_train[target_name].unique().tolist()
    if len(values) == 1:
        out_df = X_test.copy()
        if values[0] == 1:
            out_df[0] = 0
            out_df[1] = 1
            out_df = out_df[[0, 1]]
        elif values[0] == 0:
            out_df[0] = 1
            out_df[1] = 0
            out_df = out_df[[0, 1]]
        else:
            raise ValueError(f'Unexpected unique values: {values}')
    else:
        if scorer is None:
            if bin_or_cont == 'binary':
                scorer = 'neg_brier_score'
            elif bin_or_cont == 'continuous':
                scorer = 'neg_mean_absolute_error'
        if bin_or_cont == 'binary':
            cv = StratifiedKFold(n_splits=n_split_for_nested_cv, shuffle=True, random_state=53)
        elif bin_or_cont == 'continuous':
            cv = KFold(n_splits=n_split_for_nested_cv, shuffle=True, random_state=53)
        gs = GridSearchCV(
            estimator=clf_instance,
            param_grid=grid_search_param_grid,
            cv=cv,
            n_jobs=-1,
            scoring=scorer,
            verbose=0,
            refit=True
        )

        fitted_gs = gs.fit(X_train[encoding_vars], X_train[target_name])
        if bin_or_cont == 'binary':

            y_pred = fitted_gs.predict_proba(X_test[encoding_vars])
            out_df = pd.DataFrame(y_pred, index=X_test.index, columns=[0, 1])

        elif bin_or_cont == 'continuous':
            y_pred = fitted_gs.predict(X_test[encoding_vars])
            out_df = pd.DataFrame(y_pred, index=X_test.index, columns=['estimate'])

        out_df.index.name = 'accepted_species'

    out_df.to_csv(os.path.join(out_dir, f'{model_name}.csv'))


def add_y_to_data(X, case, ev_model, iteration, missingness):
    treepath, value_path = get_input_data_paths(case, ev_model, iteration)
    y = pd.read_csv(os.path.join(value_path, f'{missingness}_values.csv'), index_col=0)

    target_name = y.columns.to_list()[0]

    df = pd.merge(y, X, left_index=True, right_index=True)  # do this way to preserve ordering of y
    assert len(df) == len(y)
    df[target_name] = df.pop(target_name)  # then move column to end
    encoding_vars = X.columns.to_list()

    return df, encoding_vars, target_name


def get_umap_data(case, ev_model, iteration):
    treepath, value_path = get_input_data_paths(case, ev_model, iteration)

    X = pd.read_csv(os.path.join(treepath, 'umap_unsupervised_embedding.csv'), index_col=0)
    ## Scale the data
    X = pd.DataFrame(StandardScaler().fit_transform(X), index=X.index)

    return X


def get_semi_supervised_umap_data(case, ev_model, iteration, missingness):
    treepath, value_path = get_input_data_paths(case, ev_model, iteration)

    distances = pd.read_csv(os.path.join(treepath, 'tree_distances.csv'), index_col=0)

    full_df, encoding_vars, target_name = add_y_to_data(distances, case, ev_model, iteration, missingness)

    scaled_penguin_data = StandardScaler().fit_transform(full_df[encoding_vars])
    fitter = umap.UMAP(n_components=int(len(distances.columns) * reduction_factor))

    fitter.fit(scaled_penguin_data, y=full_df[target_name].fillna(-1, inplace=False).values)  # give NaNs a label of -1
    embedding = fitter.transform(scaled_penguin_data)

    X = pd.DataFrame(embedding, index=full_df.index)
    out_df = pd.DataFrame(StandardScaler().fit_transform(X), index=X.index, columns=X.columns)
    out_df[target_name] = full_df[target_name]

    encoding_vars = [c for c in out_df.columns if c != target_name]
    # tsne_plot(out_df, encoding_vars, target_name,
    #           os.path.join(get_prediction_data_paths(real_or_sim, bin_or_cont, iteration, missingness), 'supervised_umap.png'))
    return out_df, encoding_vars, target_name


def get_autoencoded_data(case, ev_model, iteration):
    treepath, value_path = get_input_data_paths(case, ev_model, iteration)

    X = pd.read_csv(os.path.join(treepath, 'unsupervised_autoencoded_phylogeny.csv'), index_col=0)
    ## Scale the data
    X = pd.DataFrame(StandardScaler().fit_transform(X), index=X.index)

    return X


def get_semi_supervised_autoencoded_data(case, ev_model, iteration, missingness):
    '''
    show majority samples to autoencoder to learn latent space. Idea is it learns a latent space that captures the underlying patterns of normal behavior.
    Then, when new data is projected into this latent space, deviations from normal patterns can be more easily separated using a supervised classifier.

    From approach in Gerych, Walter, Emmanuel Agu, and Elke Rundensteiner.
    ‘Classifying Depression in Imbalanced Datasets Using an Autoencoder- Based Anomaly Detection Approach’.
    In 2019 IEEE 13th International Conference on Semantic Computing (ICSC), 124–27. Newport Beach, CA, USA: IEEE, 2019.
    https://doi.org/10.1109/ICOSC.2019.8665535.

    It may not be the best way to use supervised autoencoders but it may improve on the unsupervised case.

    '''
    treepath, value_path = get_input_data_paths(case, ev_model, iteration)

    distances = pd.read_csv(os.path.join(treepath, 'tree_distances.csv'), index_col=0)

    full_df, distance_vars, target_name = add_y_to_data(distances, case, ev_model, iteration, missingness)

    majority_class = full_df[target_name].mode()[0]
    train_data = full_df[full_df[target_name] == majority_class][distance_vars]

    encoder_model, encoded_train = autoencode_pairwise_distances(train_data, reduction_fraction=reduction_factor)
    embedding = encoder_model.predict(full_df[distance_vars])

    X = pd.DataFrame(embedding, index=full_df.index)
    out_df = pd.DataFrame(StandardScaler().fit_transform(X), index=X.index, columns=X.columns)
    out_df[target_name] = full_df[target_name]

    encoding_vars = [c for c in out_df.columns if c != target_name]
    # tsne_plot(out_df, encoding_vars, target_name,
    #           os.path.join(get_prediction_data_paths(real_or_sim, bin_or_cont, iteration, missingness), 'supervised_autoenc.png'))
    return out_df, encoding_vars, target_name


def get_eigenvectors(case, ev_model, iteration):
    treepath, value_path = get_input_data_paths(case, ev_model, iteration)
    X = pd.read_csv(os.path.join(treepath, 'all_eigenvectors.csv'), index_col=0)

    broken_stick_params = pd.read_csv(os.path.join(treepath, 'broken_stick_parameters.csv'), index_col=0)
    num_cols_to_use = broken_stick_params['broken_stick_number'].iloc[0]
    X = X.iloc[:, : num_cols_to_use]
    ## Scale the data
    X = pd.DataFrame(StandardScaler().fit_transform(X), index=X.index)
    return X


def run_predictions():
    for iteration in tqdm(range(1, number_of_simulation_iterations + 1)):
        for m in missingness_types:
            for case in ['ultrametric', 'with_extinct']:
                for ev_model in simulation_types['binary'] + simulation_types['continuous']:
                    if ev_model in ['Seed Mass', 'Clonality'] and case == 'with_extinct':
                        continue

                    umap_X = get_umap_data(case, ev_model, iteration)
                    umap_df, umap_encoding_vars, umap_target_name = add_y_to_data(umap_X, case, ev_model, iteration, m)

                    eigen_X = get_eigenvectors(case, ev_model, iteration)
                    eigen_df, eigen_encoding_vars, eigen_target_name = add_y_to_data(eigen_X, case, ev_model, iteration, m)

                    autoenc_X = get_autoencoded_data(case, ev_model, iteration)
                    autoenc_df, autoenc_encoding_vars, autoenc_target_name = add_y_to_data(autoenc_X, case, ev_model, iteration, m)
                    out_dir = get_prediction_data_paths(case, ev_model, iteration, m)

                    bin_or_cont = get_bin_or_cont_from_ev_model(ev_model)
                    if bin_or_cont == 'binary':
                        # Compare logistic regression and XGBoost models i.e. for modelling simpler relationships and complex relationships
                        clf_instance = LogisticRegression(**logit_init_kwargs)
                        fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_umap', umap_df, umap_encoding_vars, umap_target_name,
                                       bin_or_cont)

                        clf_instance = LogisticRegression(**logit_init_kwargs)
                        fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_eigenvecs', eigen_df, eigen_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_umap', umap_df, umap_encoding_vars, umap_target_name,
                                       bin_or_cont)

                        clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_eigenvecs', eigen_df, eigen_encoding_vars,
                                       eigen_target_name, bin_or_cont)
                        #
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

                        # ### autoencoder
                        clf_instance = LogisticRegression(**logit_init_kwargs)
                        fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_autoencoded', autoenc_df, autoenc_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_autoencoded', autoenc_df, autoenc_encoding_vars,
                                       umap_target_name,
                                       bin_or_cont)

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

                    elif bin_or_cont == 'continuous':
                        clf_instance = LinearRegression(**linear_init_kwargs)
                        fit_and_output(clf_instance, linear_grid_search_params, out_dir, 'linear_umap', umap_df, umap_encoding_vars, umap_target_name,
                                       bin_or_cont)

                        clf_instance = LinearRegression(**linear_init_kwargs)
                        fit_and_output(clf_instance, linear_grid_search_params, out_dir, 'linear_eigenvecs', eigen_df, eigen_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        clf_instance = XGBRegressor(**xgb_rgr_init_kwargs)
                        fit_and_output(clf_instance, xgb_rgr_grid_search_params, out_dir, 'xgb_umap', umap_df, umap_encoding_vars, umap_target_name,
                                       bin_or_cont)

                        clf_instance = XGBRegressor(**xgb_rgr_init_kwargs)
                        fit_and_output(clf_instance, xgb_rgr_grid_search_params, out_dir, 'xgb_eigenvecs', eigen_df, eigen_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        ### autoencoder
                        clf_instance = LinearRegression(**linear_init_kwargs)
                        fit_and_output(clf_instance, linear_grid_search_params, out_dir, 'linear_autoencoded', autoenc_df, autoenc_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        clf_instance = XGBRegressor(**xgb_rgr_init_kwargs)
                        fit_and_output(clf_instance, xgb_rgr_grid_search_params, out_dir, 'xgb_autoencoded', autoenc_df, autoenc_encoding_vars,
                                       umap_target_name,
                                       bin_or_cont)


if __name__ == '__main__':
    run_predictions()
