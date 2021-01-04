#!/usr/bin/python3

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

#Using Chi-square to scan for categorical data (less than 5 strata) that is associated with the outcome
#Return a dictionary consisting of factor:p-value

def chi_pv(data, outcome, factor):
	table = pd.crosstab(data[outcome], data[factor])
	c, p, dof, expected = chi2_contingency(table)
	return p

def outcomescan_2(data, outcome):
    d = {}
    for factor in data.columns.values.tolist():
        if data[factor].nunique() > 5:
            continue
        elif factor == outcome:
            continue
        else:
            pv=chi_pv(data, outcome, factor)
            d[factor] = pv
    df = pd.DataFrame(d.items(), columns=['variable', 'Chisquare p-value'])
    print(df)
	
path = input('Enter data file with path:')
outcome = input('Enter surgical outcome to be analysed against:')


data = pd.read_excel(path, sheet_name='Sheet1', parse_dates=True, engine='openpyxl')

outcomescan_2(data, outcome)
