import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)

path = input('Enter data file with path:')
outcome = input('Enter surgical outcome to be analysed against:')
data = pd.read_excel(path, sheet_name='Sheet1', parse_dates=True, engine='openpyxl')

def describe(data):
	dat = []
	for var in data.columns.values.tolist():
		if data[var].nunique() > 2 and data[var].dtypes != 'O':
			des = (data[var].describe())
			dat.append(des)
	dsc = pd.concat(dat, axis=1)
	print(dsc.fillna('-'))				

def countbi(data):
	dat = []
	for var in data.columns.values.tolist():
		if data[var].nunique() == 2:
			des = data.groupby(var)[var].count()
			dat.append(des)
	dsc = pd.concat(dat, axis = 1)
	sum = dsc.aggregate('sum', axis = 0)
	dsc = dsc.append(sum, ignore_index=True)
	lst_a = dsc.index.tolist()
	lst_a[-1] = 'vertisum'
	dsc.index = lst_a

	cols = dsc.columns.tolist()
	for i in cols:
		dsc['%'+i] = dsc[i]/dsc[i].vertisum*100
	print(dsc.fillna('-'))


	

def mean_by_outcome(data):
	mbo = []
	for var in data.columns.values.tolist():
		if data[var].nunique() > 2 and data[var].dtypes != 'O':
			mbo.append(data.groupby(outcome).mean()[var])
	des = pd.concat(mbo, axis = 1)
	print(des)

print('General summary statistics:\n')	
describe(data)
print('Count binary dara:\n')
countbi(data)
print ('Mean by outcome group \n')
mean_by_outcome(data)