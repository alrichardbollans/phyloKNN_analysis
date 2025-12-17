import os
import pathlib

import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import spearmanr

from evaluation.compare_types import plot_binary_and_continuous_cases
from evaluation.compile_prediction_scores import read_all_results
import seaborn as sns

from evaluation.do_ttests_for_results import ttests, hochberg_correction
from evaluation.helper_functions import bin_model_names, cont_model_names, rename_models_and_ev_models


def plot_signal_performance(df, model, bin_or_cont, out_dir, drop_outlier=False):
    sns.set_theme()
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    if drop_outlier:
        # For readability drop this outlier
        if bin_or_cont == 'binary':
            df = df[df['Signal'] <= 30]
        else:
            raise ValueError('Only binary models have outliers')
    df = df.rename(columns=rename_models_and_ev_models)
    if model in rename_models_and_ev_models:
        model1_renamed = rename_models_and_ev_models[model]
    else:
        model1_renamed = model

    plt.figure(figsize=(12, 8))
    sns.regplot(data=df, x='Signal', y=model1_renamed)
    if bin_or_cont == 'binary':

        plt.ylim(0, 1)
        plt.ylabel('Brier Score')
        plt.xlabel("Fritz and Purvis' D")
    else:

        plt.ylabel('Mean Absolute Error')
        plt.xlabel("Pagel's λ")
    plt.savefig(os.path.join(out_dir, f'{model}_{str(drop_outlier)}.jpg'), dpi=300)
    plt.close()

    plt.figure(figsize=(12, 8))
    sns.regplot(data=df[df['Missing Type'] == 'mcar'], x='Signal', y=model1_renamed)
    if bin_or_cont == 'binary':

        plt.ylim(0, 1)
        plt.ylabel('Brier Score')
        plt.xlabel("Fritz and Purvis' D")
    else:

        plt.ylabel('Mean Absolute Error')
        plt.xlabel("Pagel's λ")
    plt.savefig(os.path.join(out_dir, f'{model}_mcar_{str(drop_outlier)}.jpg'), dpi=300)
    plt.close()

    plt.figure(figsize=(12, 8))
    sns.regplot(data=df[df['Missing Type'] == 'phyloNa'], x='Signal', y=model1_renamed)
    if bin_or_cont == 'binary':

        plt.ylim(0, 1)
        plt.ylabel('Brier Score')
        plt.xlabel("Fritz and Purvis' D")
    else:

        plt.ylabel('Mean Absolute Error')
        plt.xlabel("Pagel's λ")
    plt.savefig(os.path.join(out_dir, f'{model}_phyloNa_{str(drop_outlier)}.jpg'), dpi=300)
    plt.close()


def plot_signal_performance_two_models(df, model1, model2, bin_or_cont, out_dir, drop_outlier=False):
    sns.set_theme()

    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    if drop_outlier:
        # For readability drop this outlier
        if bin_or_cont == 'binary':
            df = df[df['Signal'] <= 30]
        else:
            raise ValueError('Only binary models have outliers')


    df = df.rename(columns=rename_models_and_ev_models)
    if model1 in rename_models_and_ev_models:
        model1_renamed = rename_models_and_ev_models[model1]
    else:
        model1_renamed = model1
    if model2 in rename_models_and_ev_models:
        model2_renamed = rename_models_and_ev_models[model2]
    else:
        model2_renamed = model2
    sns.regplot(data=df, x='Signal', y=model1_renamed, label=model1_renamed, color='b',
                scatter_kws={'alpha':0.6})
    sns.regplot(data=df, x='Signal', y=model2_renamed, label=model2_renamed, color='orange',
                scatter_kws={'alpha':0.6})
    plt.legend()

    if bin_or_cont == 'binary':

        plt.ylim(-0.1, 1)
        plt.ylabel('Brier Score')
        plt.xlabel("Fritz and Purvis' D")
    else:

        plt.ylabel('Mean Absolute Error')
        plt.xlabel("Pagel's λ")

    plt.savefig(os.path.join(out_dir, f'{model1}_{model2}_{str(drop_outlier)}.jpg'), dpi=300)
    plt.close()

def main():
    out_dir = 'signal_analysis'
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    bin_df, cont_df = read_all_results()
    # # plot signal distribution
    sns.displot(bin_df.reset_index(drop=True), x='Signal', kind='kde')
    plt.savefig(os.path.join(out_dir, f'D.jpg'), dpi=300)
    plt.close()

    sns.displot(cont_df, x='Signal', kind='kde')
    plt.savefig(os.path.join(out_dir, f'lambda.jpg'), dpi=300)
    plt.close()

    # Calculate correlations
    spearman_values = []
    spearman_indices = []
    for model in bin_model_names:
        model_correlation = spearmanr(bin_df['Signal'], bin_df[model], nan_policy='omit')
        spearman_values.append([model_correlation.statistic, model_correlation.pvalue])
        spearman_indices.append(f'{model}_bin')
    for model in cont_model_names:
        model_correlation = spearmanr(cont_df['Signal'], cont_df[model], nan_policy='omit')
        spearman_values.append([model_correlation.statistic, model_correlation.pvalue])
        spearman_indices.append(f'{model}_cont')
    spearmanr_df = pd.DataFrame(spearman_values, index=spearman_indices,
                                columns=['spearman_values', 'pvalues'])
    h_sp_df = hochberg_correction(spearmanr_df, 'pvalues')
    h_sp_df.to_csv(os.path.join(out_dir, 'spearman_correlations.csv'))
    ## plot performance against signal for each model
    for model in bin_model_names:
        plot_signal_performance(bin_df, model, 'binary', os.path.join(out_dir, 'binary models', model))
        plot_signal_performance(bin_df, model, 'binary', os.path.join(out_dir, 'binary models', model), drop_outlier=True)

        if model != 'corHMM':
            # plot performance against signal for corHMM and a ML model
            plot_signal_performance_two_models(bin_df, 'corHMM', model, 'binary', os.path.join(out_dir, 'models vs corHMM'))
            plot_signal_performance_two_models(bin_df, 'corHMM', model, 'binary', os.path.join(out_dir, 'models vs corHMM'), drop_outlier=True)

    for model in cont_model_names:
        plot_signal_performance(cont_df, model, 'continuous', os.path.join(out_dir, 'continuous models', model))

        if model != 'phylopars':
            # plot performance against signal for corHMM and a ML model
            plot_signal_performance_two_models(cont_df, 'phylopars', model, 'continuous', os.path.join(out_dir, 'models vs phylopars'))

    # #  phylogenetic imputation may perform poorly when lambda is less than 0.6 (Molina-Venegas et al., 2018)
    # ## compare models with signal above 0.6 and signal below 0.6
    binary_above_06_df = bin_df[bin_df['Signal'] <= 0.5]
    continuous_above_06_df = cont_df[cont_df['Signal'] >= 0.6]
    plot_binary_and_continuous_cases(binary_above_06_df, continuous_above_06_df, 'Ev Model',
                                     os.path.join(out_dir, 'strong signal'))
    binary_above_06_df.describe(include='all').to_csv(os.path.join(out_dir, 'strong signal', 'binary_results_summary.csv'))
    continuous_above_06_df.describe(include='all').to_csv(os.path.join(out_dir, 'strong signal', 'continuous_results_summary.csv'))

    binary_below_06_df = bin_df[bin_df['Signal'] > 0.5]
    continuous_below_06_df = cont_df[cont_df['Signal'] < 0.6]
    plot_binary_and_continuous_cases(binary_below_06_df, continuous_below_06_df, 'Ev Model', os.path.join(out_dir, 'weak signal'))
    binary_below_06_df.describe(include='all').to_csv(os.path.join(out_dir, 'weak signal', 'binary_results_summary.csv'))
    continuous_below_06_df.describe(include='all').to_csv(os.path.join(out_dir, 'weak signal', 'continuous_results_summary.csv'))

    ttests(binary_below_06_df, bin_model_names, os.path.join(out_dir, 'weak signal', 'binary'))
    ttests(continuous_below_06_df, cont_model_names, os.path.join(out_dir, 'weak signal', 'continuous'))



if __name__ == '__main__':
    main()
