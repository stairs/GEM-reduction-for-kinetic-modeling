from pandas import read_csv
from scipy.io import loadmat


def load_reference_fccs(path: str):
    result = read_csv(path, delimiter='\t', index_col=0)
    result.columns = result.columns.str[1:-1]
    result.index = result.index.str[1:-1]
    result = result.drop(columns='Summation Error')
    return result


def load_sampled_fccs(path: str):
    data = loadmat(path)
    return data['CJ_rec']

def load_sampled_parameters(path: str):
    data = loadmat(path)
    return data['Parameters']