import os
import pathlib
from itertools import combinations

import pandas as pd
from scipy.stats import ttest_rel

from evaluation.compile_prediction_scores import read_all_results
from evaluation.helper_functions import bin_model_names, cont_model_names

_alpha = 0.05


def ttests(full_df: pd.DataFrame, model_names, out_dir):
    pathlib.Path(out_dir).mkdir(exist_ok=True, parents=True)

    ttest_dict = {}
    for model_name2 in model_names:
        if 'corHMM' in model_names:
            model_name1 = 'corHMM'
        if 'phylopars'in model_names:
            model_name1 = 'phylopars'
        if model_name1 != model_name2:
            if model_name1 in full_df.columns and model_name2 in full_df.columns:
                t_stat, p_value = ttest_rel(full_df[model_name1], full_df[model_name2], nan_policy='omit')
                ttest_dict[f'{model_name1}_{model_name2}'] = [t_stat, p_value]
    # Convert to a DataFrame
    ttest_df = pd.DataFrame(ttest_dict, index=['stat', 'p value'])

    # Save to CSV
    reformatted_ttest_df = []
    for c in ttest_df.columns:
        # if test_results[c].loc['p value']<_alpha:
        reformatted_ttest_df.append([c, ttest_df[c].loc['stat'], ttest_df[c].loc['p value']])

    # raw_test_results = pd.read_csv(os.path.join(eval_dir, 'raw_ttest_results.csv'), index_col=0)
    # for c in raw_test_results.columns:
    #     # if raw_test_results[c].loc['p value']<_alpha:
    #         significant_cases.append([c, raw_test_results[c].loc['stat'],raw_test_results[c].loc['p value']])
    reformatted_ttest_df = pd.DataFrame(reformatted_ttest_df, columns=['models', 'stat', 'p value'])

    corrected_df = hochberg_correction(reformatted_ttest_df, p_value_col='p value')
    corrected_df.to_csv(os.path.join(out_dir, 'ttest_summary.csv'))


def hochberg_correction(df: pd.DataFrame, p_value_col: str):
    # Hochberg Method
    # Yosef Hochberg, ‘A Sharper Bonferroni Procedure for Multiple Tests of Significance’, Biometrika 75, no. 4 (1988): 800–802,
    # https://doi.org/10.1093/biomet/75.4.800.

    # The adjustment is the same as holm-bonferroni, but the usage is slightly different.
    new_df = df.sort_values(by=p_value_col)
    n = len(new_df.index)
    new_df.reset_index(inplace=True, drop=True)

    new_df['hochberg_adjusted_p_value'] = new_df.apply(lambda x: x[p_value_col] * (n - x.name), axis=1)
    new_df = new_df.sort_values(by='hochberg_adjusted_p_value')

    return new_df


def summarise_ttests():
    bin_df, cont_df = read_all_results()
    ttests(bin_df, bin_model_names, os.path.join('ttest', 'binary'))
    ttests(cont_df, cont_model_names, os.path.join('ttest', 'continuous'))

    for ev_model in bin_df['Ev Model']:
        ttests(bin_df[bin_df['Ev Model']==ev_model], bin_model_names, os.path.join('ttest','binary', ev_model))
    for ev_model in cont_df['Ev Model']:
        ttests(cont_df[cont_df['Ev Model']==ev_model], cont_model_names, os.path.join('ttest','continuous', ev_model))

if __name__ == '__main__':
    summarise_ttests()
