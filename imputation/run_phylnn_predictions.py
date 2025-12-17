import sys

from tqdm import tqdm
sys.path.append('../..')
from data.helper_functions import number_of_simulation_iterations, simulation_types
from imputation.helper_functions import phylnn_predict, missingness_types


def main():
    for tag in tqdm(range(1, number_of_simulation_iterations + 1)):
        for m in missingness_types:
            for case in ['ultrametric', 'with_extinct']:
                for ev_model in simulation_types['binary'] + simulation_types['continuous']:
                    if ev_model in ['Seed Mass', 'Clonality'] and case == 'with_extinct':
                        continue
                    phylnn_predict(case, ev_model, tag, m)


if __name__ == '__main__':
    main()
