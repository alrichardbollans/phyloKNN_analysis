import os.path

repo_path = os.environ.get('KEWSCRATCHPATH')
input_data_path = os.path.join(repo_path, 'phyloKNN_analysis', 'data')
prediction_path = os.path.join(repo_path, 'phyloKNN_analysis', 'imputation')

number_of_simulation_iterations = 100

reduction_factor = 0.1

simulation_types = {'binary':['ER', 'ARD', 'BiSSE', 'HiSSE', 'bBMT', 'Clonality'],
                    'continuous':['BM', 'OU', 'EB', 'LB', 'BMT', 'Seed Mass']}


