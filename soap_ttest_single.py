#!/usr/bin/python3

#Independent T-test (unpaired T-test, Non-repeated measure)
#Performing test for equality between 2-group mean (H0) with an assuption that the 2 groups have normal distribution of values. If the H0 is rejected (p-value < 0.05), the alternative hypothesis is that the mean of one group is either greater or lesser than the other group (2-tail test).
#For 2 sets of data with homogeniety in their distribution (equal variance).
#All factors to be analysed against the outcome must be binary and coded with 0 and 1 only.

import pandas as pd
import scipy
from scipy import stats
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

path = input('Enter data file with path:')
data = pd.read_excel(path, sheet_name='Sheet1', parse_dates=True, engine='openpyxl')
outcome = input('Enter surgical outcome to be analysed against:')
var = input('Variable to be tested, group by outcome:')

def t_test_multi(data, outcome, var):
	if data[outcome].nunique() != 2:
		print('The outcome', outcome, 'is non-binary')
	else:
		var_by_outcome = data.groupby(outcome)[var].describe()
		print(var,'\n', var_by_outcome)
		cat1 = data[data[outcome] == 0]
		cat2 = data[data[outcome] == 1]
		print('-----------------------------------------------------------------\n')
		print(stats.ttest_ind(cat1[var].dropna(), cat2[var].dropna()))
		print(stats.mannwhitneyu(cat1[var].dropna(), cat2[var].dropna()))
			
		plt.boxplot([cat1[var].dropna(), cat2[var].dropna()], showmeans = True)
		plt.show()

t_test_multi(data, outcome, var)






