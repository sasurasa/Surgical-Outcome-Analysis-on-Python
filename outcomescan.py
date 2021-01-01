import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency


def chi_pv(data, outcome, factor):
	table = pd.crosstab(data[outcome], data[factor])
	c, p, dof, expected = chi2_contingency(table)
	return p

def outcomescan_2(data, outcome):
    d = {}
    for factor in data.columns.values.tolist():
        if data[factor].nunique() > 5:
            continue
        else:
            pv=chi_pv(data, outcome, factor)
            d[factor] = pv
    print(d)