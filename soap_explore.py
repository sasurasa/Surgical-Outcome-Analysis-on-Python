#!/usr/bin/python3

#Importing Excel file into a pandas' dataframe format

import pandas as pd

path = input('Enter data file with path:')
outcome = input('Enter surgical outcome to be analysed against:')


data = pd.read_excel(path, sheet_name='Sheet1', parse_dates=True, engine='openpyxl', index_col=outcome)

def soap_explore(data):
	size = data.size
	dimension = data.shape
	variables = data.columns.values.tolist()
	print('===============================================================================================================')
	print('\nThe dataframe has',size, 'cells, with',dimension,'(row x column) dimension.\n')
	print('All the variables include;',variables,'\n')
	print('===============================================================================================================')
	print('Types and numbers of each variables are;')
	print(data.info())

soap_explore(data)

from scipy import stats
from scipy.stats import shapiro


def shapif(data):
    d1 = {}
    for factor in data.columns.values.tolist():
        if data[factor].dtypes == 'float64':
        	x = data[factor].notna()
        	s = stats.shapiro(x)
        	d1[factor] = s
    df1 = pd.DataFrame(d1.items(), columns=['variable', 'Shapiro-Wilk result'])
    print(df1)

def shapii(data):
    d2 = {}
    for factor in data.columns.values.tolist():
        if data[factor].dtypes == 'int64':
        	x = data[factor]
        	s = stats.shapiro(x)
        	d2[factor] = s
    df2 = pd.DataFrame(d2.items(), columns=['variable', 'Shapiro-Wilk result'])
    print(df2)

shapif(data)
print('\n')
shapii(data)
