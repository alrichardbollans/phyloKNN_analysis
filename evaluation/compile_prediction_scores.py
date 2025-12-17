import os
import pathlib

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.metrics import brier_score_loss, mean_absolute_error, average_precision_score

from data.helper_functions import number_of_simulation_iterations, simulation_types
from evaluation.helper_functions import bin_model_names, cont_model_names, rename_models_and_ev_models, \
    binary_model_order, continuous_model_order
from imputation.helper_functions import get_input_data_paths, get_prediction_data_paths, missingness_types, get_bin_or_cont_from_ev_model

_iterations_still_to_run = []


def check_prediction_data(dfs: list[pd.DataFrame], ground_truth: pd.DataFrame, missing_values: pd.DataFrame):
    assert missing_values.shape[0] == ground_truth.shape[0]
    assert missing_values.shape[1] == ground_truth.shape[1]
    assert ground_truth.columns[0] == missing_values.columns[0]

    assert len(ground_truth.columns) == 1

    gt_target_name = ground_truth.columns[0]
    missing_values = missing_values[missing_values[gt_target_name].isna()]

    for df in dfs:
        try:
            pd.testing.assert_index_equal(missing_values.index, df.index, check_names=False)

        except AssertionError:
            test_species = missing_values.index.tolist()
            pred_species = df.index.tolist()
            issues = [c for c in test_species if c not in pred_species]
            if issues != ['×_Staparesia_meintjesii', '×_Stapvalia_oskopensis']:
                # print(issues)
                # print('##############')
                # print([c for c in pred_species if c not in test_species])
                raise AssertionError(f'Model issue {df.columns[0]}')


def get_model_names(bin_or_cont):
    if bin_or_cont == 'binary':

        return bin_model_names
    elif bin_or_cont == 'continuous':
        return cont_model_names
    elif bin_or_cont == 'both':
        return list(set(bin_model_names + cont_model_names))
    else:
        raise ValueError(f'Unknown data type {bin_or_cont}')


def get_model_predictions(case: str, ev_model: str, iteration: int, missing_type: str, ground_truth, missing_values):
    imputation_path = get_prediction_data_paths(case, ev_model, iteration, missing_type)
    bin_or_cont = get_bin_or_cont_from_ev_model(ev_model)
    model_names = get_model_names(bin_or_cont)
    dfs = []
    for model_name in model_names:
        try:
            model_df = pd.read_csv(os.path.join(imputation_path, f'{model_name}.csv'), index_col=0)
            if bin_or_cont == 'binary':
                assert len(model_df.columns) == 2
                model_df = model_df[['1']]
            elif bin_or_cont == 'continuous':
                assert len(model_df.columns) == 1
            model_df.columns = [model_name]
            dfs.append(model_df)
        except FileNotFoundError:
            print(f'Missing {model_name} for {ev_model} {bin_or_cont} {iteration} {missing_type}')
            _iterations_still_to_run.append(iteration)
    check_prediction_data(dfs, ground_truth, missing_values)

    try:
        check_prediction_data(dfs, ground_truth, missing_values)
    except AssertionError as m:
        raise AssertionError(f'{m}. Issue with {ev_model}: str, {bin_or_cont}: str, {iteration}: int, {missing_type}. ')
    return dfs


def evaluate_model_outputs_for_specific_case(case: str, ev_model: str, iteration: int, missing_type: str, scorer=None):
    treepath, data_path = get_input_data_paths(case, ev_model, iteration)
    out_dict = {}

    # Get some information about test scenario
    bin_or_cont = get_bin_or_cont_from_ev_model(ev_model)

    out_dict['Ev Model'] = ev_model
    out_dict['Missing Type'] = missing_type
    out_dict['Tree Type'] = case
    out_dict['Variable Type'] = bin_or_cont

    if bin_or_cont == 'binary':
        signal_df = pd.read_csv(os.path.join(data_path, 'phylogenetic_signal_results_D.csv'), index_col=0)

    else:
        signal_df = pd.read_csv(os.path.join(data_path, 'phylogenetic_signal_results_lambda.csv'), index_col=0)
    signal = signal_df['value'].iloc[0]
    out_dict['Signal'] = signal
    # Compile predictions on test data with ground truth
    ground_truth = pd.read_csv(os.path.join(data_path, 'ground_truth.csv'), index_col=0)
    assert len(ground_truth.columns) == 1
    missing_values = pd.read_csv(os.path.join(data_path, f'{missing_type}_values.csv'), index_col=0)

    dfs = get_model_predictions(case, ev_model, iteration, missing_type, ground_truth, missing_values)

    gt_target_name = ground_truth.columns[0]
    missing_values = missing_values[missing_values[gt_target_name].isna()]
    test_ground_truths = ground_truth.loc[missing_values.index]
    full_df = test_ground_truths

    # Compile model predictions
    for df in dfs:
        full_df = pd.merge(full_df, df, left_index=True, right_index=True)

    # Check some potential issues
    issues = full_df[full_df['phylnn_fill_means'].isna()]
    assert len(issues) == 0

    # if drop_nans:
    #     full_df = full_df[full_df['phylnn_raw'].notna()]
    #     pd.testing.assert_series_equal(full_df['phylnn_raw'], full_df['phylnn_fill_means'], check_names=False)
    #     full_df = full_df.drop('phylnn_fill_means', axis=1)
    # else:
    #     if 'phylnn_raw' in full_df:
    #         full_df = full_df.drop('phylnn_raw', axis=1)

    # Get scores for model names
    model_names = get_model_names(bin_or_cont)
    if len(full_df) > 0:
        for model_name in model_names:
            if model_name in full_df.columns:
                if scorer is None:
                    if bin_or_cont == 'binary':
                        score = brier_score_loss(full_df[gt_target_name], full_df[model_name])
                    elif bin_or_cont == 'continuous':
                        score = mean_absolute_error(full_df[gt_target_name], full_df[model_name])
                else:
                    print(f'Using {scorer.__name__} for {model_name}')
                    score = scorer(full_df[gt_target_name], full_df[model_name])
                out_dict[model_name] = score
    return out_dict


def plot_results(df, model_names, out_dir, scorer_label):
    plot_df = df[model_names]
    sns.boxplot(data=plot_df)
    plt.xticks(rotation=30, ha='right')
    plt.ylabel(scorer_label)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'violin_plot.jpg'), dpi=300)
    plt.clf()
    plt.close()


def output_results_from_df(full_df: pd.DataFrame, out_dir: str, ev_model: str, scorer_label='Loss'):
    full_df = full_df.reset_index(drop=True)
    pathlib.Path(out_dir).mkdir(exist_ok=True, parents=True)

    full_df.to_csv(os.path.join(out_dir, f'results.csv'))
    full_df.describe(include='all').to_csv(os.path.join(out_dir, f'results_summary.csv'))
    bin_or_cont = get_bin_or_cont_from_ev_model(ev_model)
    model_names = get_model_names(bin_or_cont)

    plot_results(full_df, [c for c in model_names if c in full_df.columns], out_dir, scorer_label=scorer_label)


def collate_simulation_outputs(ev_model: str,
                               range_to_eval: int = number_of_simulation_iterations, scorer=None, out_dir=None):
    full_df = pd.DataFrame()
    for tag in range(1, range_to_eval + 1):
        for case in ['ultrametric', 'with_extinct']:
            if ev_model in ['Seed Mass', 'Clonality', 'APM'] and case == 'with_extinct':
                continue

            for missing_type in missingness_types:
                if ev_model == 'APM' and missing_type == 'phyloNa':
                    continue
                run_dict = evaluate_model_outputs_for_specific_case(case, ev_model, tag, missing_type, scorer=scorer)
                run_df = pd.DataFrame(run_dict, index=[tag])
                full_df = pd.concat([full_df, run_df])
    if out_dir is None:
        out_dir = os.path.join('compiled_score_outputs', ev_model)
    output_results_from_df(full_df, out_dir, ev_model, scorer_label=scorer.__name__ if scorer is not None else 'Loss')

    return full_df


def read_all_results(output_too=False):
    binary_df = pd.DataFrame()
    for ev_model in simulation_types['binary']:
        ev_model_df = pd.read_csv(os.path.join('compiled_score_outputs', ev_model, 'results.csv'))
        binary_df = pd.concat([binary_df, ev_model_df])

    continuous_df = pd.DataFrame()
    for ev_model in simulation_types['continuous']:
        ev_model_df = pd.read_csv(os.path.join('compiled_score_outputs', ev_model, 'results.csv'))
        continuous_df = pd.concat([continuous_df, ev_model_df])

    if output_too:
        binary_df.describe(include='all').to_csv(os.path.join('compiled_score_outputs', 'binary_results_summary.csv'))
        continuous_df.describe(include='all').to_csv(os.path.join('compiled_score_outputs', 'continuous_results_summary.csv'))
        sns.set_theme()
        err_kws = {"color": ".2", "linewidth": 1.5}
        p_df = pd.melt(binary_df, value_vars=bin_model_names, var_name='Model', value_name='Brier Score')
        p_df['Model'] = p_df['Model'].map(rename_models_and_ev_models).fillna(p_df['Model'])
        g = sns.barplot(p_df, x='Model', y='Brier Score', order=binary_model_order,
                        capsize=.4, err_kws=err_kws)
        g.set_xticklabels(g.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        plt.tight_layout()
        plt.savefig(os.path.join('compiled_score_outputs', 'binary_means.jpg'), dpi=300)
        plt.close()

        p_df = pd.melt(continuous_df, value_vars=cont_model_names, var_name='Model', value_name='Mean Absolute Error')
        p_df['Model'] = p_df['Model'].map(rename_models_and_ev_models).fillna(p_df['Model'])
        g = sns.barplot(p_df, x='Model', y='Mean Absolute Error', order=continuous_model_order,
                        capsize=.4,
                        err_kws=err_kws,
                        )
        g.set_xticklabels(g.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        plt.tight_layout()
        plt.savefig(os.path.join('compiled_score_outputs', 'continuous_means.jpg'), dpi=300)
        plt.close()

    return binary_df, continuous_df


def main():
    # collate_simulation_outputs('APM', range_to_eval=10, out_dir='APM_outputs', scorer=average_precision_score)
    for ev_model in simulation_types['binary'] + simulation_types['continuous']:
        collate_simulation_outputs(ev_model)

    print(set(_iterations_still_to_run))
    read_all_results(output_too=True)


if __name__ == '__main__':
    main()
