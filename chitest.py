#Setting Chi-square test for association between a certain factor with surgical outcome

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

#data: python Dataframe
#outcome: binary data
#factor: should be binary as well
path = 'Desktop/..'
data = pd.read_excel(path, sheet_name='Sheet1', parse_dates = True)

#Input: data, outcome and the factor to be tested
#Output: table dimension and Chi-square p-value

def chi_p(data, outcome, factor):
	table = pd.crosstab(data[outcome], data[factor])
	c, p, dof, expected = chi2_contingency(table)
	print(table.shape)
	print(p)
