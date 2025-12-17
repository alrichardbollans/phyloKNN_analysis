import os.path
import pathlib
import pickle
import sys

import pandas as pd
from sklearn.metrics import make_scorer, mean_absolute_error, brier_score_loss
from sklearn.model_selection import KFold
sys.path.append('../..')
from data.helper_functions import input_data_path, simulation_types
from phylokNN import nan_safe_metric_wrapper, phyloNN_bayes_opt, PhylNearestNeighbours

repo_path = os.environ.get('KEWSCRATCHPATH')
prediction_path = os.path.join(repo_path, 'phyloKNN', 'analysis', 'imputation')

missingness_types = ['mcar', 'phyloNa']

n_split_for_nested_cv = 5

def get_iteration_path_from_base(base: str, case: str, ev_model: str, iteration: int):
    if ev_model in ['BiSSE', 'HiSSE']:
        basepath = os.path.join(base,'simulations', case, ev_model)

    elif ev_model in ['Clonality', 'Seed Mass']:
        basepath = os.path.join(base, 'real_data', case, ev_model)
    elif ev_model in ['APM']:
        basepath = os.path.join(base, 'my_apm_data', case, ev_model)

    else:
        basepath = os.path.join(base, 'simulations', case, 'standard')

    treepath = os.path.join(basepath, str(iteration))
    value_path = os.path.join(basepath, str(iteration), ev_model)

    return treepath, value_path


def get_input_data_paths(case: str, ev_model: str, iteration: int):
    return get_iteration_path_from_base(input_data_path, case, ev_model, iteration)


def get_prediction_data_paths(case: str, ev_model: str, iteration: int, missingness_type: str):
    return os.path.join(prediction_path, case, ev_model, str(iteration), missingness_type)


def check_data(ground_truth, missing_values):
    assert len(ground_truth) == len(missing_values)
    assert len(ground_truth.columns) == len(missing_values.columns)

    def check_df(df):
        index_name = df.columns[0]
        assert index_name == 'accepted_species'

        assert len(df.columns) == 2

    check_df(ground_truth)
    check_df(missing_values)

    target_name = missing_values.columns[1]

    mcar_nans = missing_values[missing_values[target_name].isna()]
    assert len(mcar_nans) > 1
    assert len(mcar_nans) < len(missing_values)


def get_bin_or_cont_from_ev_model(ev_model: str):
    if ev_model in simulation_types['binary'] or ev_model=='APM':
        return 'binary'
    elif ev_model in simulation_types['continuous']:
        return 'continuous'
    else:
        raise ValueError(f'Unknown data type {ev_model}')

def phylnn_predict(case: str, ev_model: str, iteration: int, missing_type: str, val_scorer=None):
    treepath, value_path = get_input_data_paths(case, ev_model, iteration)
    ground_truth = pd.read_csv(os.path.join(value_path, 'ground_truth.csv'))
    missing_values = pd.read_csv(os.path.join(value_path, f'{missing_type}_values.csv'))
    check_data(ground_truth, missing_values)

    bin_or_cont = get_bin_or_cont_from_ev_model(ev_model)

    distance_csv = os.path.join(treepath, 'tree_distances.csv')
    if val_scorer is None:
        if bin_or_cont == 'continuous':
            val_scorer = make_scorer(nan_safe_metric_wrapper(mean_absolute_error), greater_is_better=False)
        elif bin_or_cont == 'binary':
            val_scorer = make_scorer(nan_safe_metric_wrapper(brier_score_loss), greater_is_better=False, response_method='predict_proba')
        else:
            raise ValueError(f'Unknown data type {bin_or_cont}')
    if bin_or_cont == 'continuous':
        clf = False
    elif bin_or_cont == 'binary':
        clf = True
    else:
        raise ValueError(f'Unknown data type {bin_or_cont}')
    out_dir = get_prediction_data_paths(case, ev_model, iteration, missing_type)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    distance_df = pd.read_csv(distance_csv, index_col=0)
    # reduce size of distance dataframe to match trait data
    cols_for_distance_df = [c for c in distance_df.columns if c in missing_values[missing_values.columns[0]].values]
    distance_df = distance_df[cols_for_distance_df]
    distance_df = distance_df.loc[cols_for_distance_df]

    njobs = -1
    verbose = 0

    target_name = missing_values.columns[1]
    train = missing_values[~missing_values[target_name].isna()]

    best_ratio, best_kappa = phyloNN_bayes_opt(
        distance_df,
        clf=clf,
        scorer=val_scorer, cv=KFold(n_splits=n_split_for_nested_cv, shuffle=True, random_state=42), X=train, y=train[target_name].values, njobs=njobs,
        verbose=verbose)

    best_phyln_fill_means = PhylNearestNeighbours(distance_df, clf, ratio_max_branch_length=best_ratio, kappa=best_kappa,
                                                  fill_in_unknowns_with_mean=True)
    best_phyln_fill_means.fit(train, train[target_name].values)

    pickle.dump(best_phyln_fill_means, open(os.path.join(out_dir, 'phylnn_fill_means_hparams.pkl'), 'wb'))

    test = missing_values[missing_values[target_name].isna()]
    assert len(set(train.index).intersection(set(test.index))) == 0
    if clf:
        prediction_fill_means = best_phyln_fill_means.predict_proba(test)
    else:
        prediction_fill_means = best_phyln_fill_means.predict(test)

    pd.DataFrame(prediction_fill_means, index=test[test.columns[0]]).to_csv(os.path.join(out_dir, 'phylnn_fill_means.csv'))

    ## Without filling means
    best_phyln_raw = PhylNearestNeighbours(distance_df, clf, ratio_max_branch_length=best_ratio, kappa=best_kappa,
                                           fill_in_unknowns_with_mean=False)
    best_phyln_raw.fit(train, train[target_name].values)

    pickle.dump(best_phyln_raw, open(os.path.join(out_dir, 'phylnn_raw_hparams.pkl'), 'wb'))

    if clf:
        prediction_raw = best_phyln_raw.predict_proba(test)
    else:
        prediction_raw = best_phyln_raw.predict(test)

    pd.DataFrame(prediction_raw, index=test[test.columns[0]]).to_csv(os.path.join(out_dir, 'phylnn_raw.csv'))
