import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ttest_ind

path = input('Enter data file with path:')
data = pd.read_excel(path, sheet_name='Sheet1', parse_dates=True, engine='openpyxl')

var = input('Continuous variable to be tested:')
bivar = list(map(str, input('Enter a list of binary variables to be used for comparisons of var,separated by a space: ').split()))

def tmbi(data, bivar, var):
	edb = []
	for i in range(len(bivar)):
		eb = data.groupby(bivar[i])[var].mean()
		edb.append(eb)
		cat1 = data[data[bivar[i]] == 0]
		cat2 = data[data[bivar[i]] == 1]
		print('T-test for '+ var + ' by ' + bivar[i])
		print(stats.ttest_ind(cat1[var].dropna(), cat2[var].dropna()))
	a = pd.concat(edb, axis =1)
	a.columns = bivar
	a = a.add_prefix(var+'_')
	
	print ('=============================================================================')
	print('Mean '+ var + ' by each binary factors')
	print(a)
	

def mmbi(data, bivar, var):
	edb2 = []	
	for i in range(len(bivar)):
		eb2 = data.groupby(bivar[i])[var].median()
		edb2.append(eb2)
		cat1 = data[data[bivar[i]] == 0]
		cat2 = data[data[bivar[i]] == 1]
		print('Ranksum-test for '+ var + ' by ' + bivar[i])
		print(stats.mannwhitneyu(cat1[var].dropna(), cat2[var].dropna()))
	b = pd.concat(edb2, axis =1)
	b.columns = bivar
	b = b.add_prefix(var+'_')

	print ('=============================================================================')
	print('Median '+ var + ' by each binary factors')
	print(b)

tmbi(data, bivar, var)
print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
mmbi(data, bivar, var)
	
	
	