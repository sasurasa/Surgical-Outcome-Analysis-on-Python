#!/usr/bin/python3

#Importing Excel file into a pandas' dataframe format

import pandas as pd

path = input('Enter data file with path:')
path = path + ".xlsx"
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