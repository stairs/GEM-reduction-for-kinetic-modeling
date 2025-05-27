from typing import List
import numpy as np
import pandas as pd


def exclude_reactions(fcc_samples, all_rxns: List[str], exclude_rxns: List[str]):
    indices = [index for (index, rxn) in enumerate(all_rxns) if rxn not in exclude_rxns]
    return fcc_samples[np.ix_(indices, indices)], all_rxns[indices]


def exclude_causes(fcc_samples, fcc_ref: pd.DataFrame, exclude_rxns: List[str]):
    causes = fcc_ref.columns
    fcc_ref.skew()
    indices = [index for (index, rxn) in enumerate(causes) if rxn not in exclude_rxns]
    return fcc_samples[:, indices], fcc_ref.iloc[:, indices]
