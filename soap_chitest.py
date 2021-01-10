#!/usr/bin/python3


#Setting Chi-square test for association between a certain factor with surgical outcome

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

#data: python Dataframe
#outcome: binary data
#factor: should be binary as well
path = input('Excel path: ')
data = pd.read_excel(path, sheet_name='Sheet1', parse_dates = True, engine='openpyxl')
var_b = input('Outcome to be analysed: ')
var_a = input('Variable to be tested against outcome: ')

def chi_p(data, var_b, var_a):
	table = pd.crosstab(data[var_b], data[var_a])
	c, p, dof, expected = chi2_contingency(table)
	print('Data dimension: ', table.shape)
	print('Chisquare p-value =',p)
print('=================================================================================\n')
chi_p(data, var_b, var_a)
print('=================================================================================\n')
twosub = data[[var_a, var_b]]
var_b_list = data[var_b].unique().tolist()
all_col_list = []
for i in var_b_list:
	Bi = twosub[data[var_b] == i]
	Ci = Bi.groupby(var_a).count()
	dict = {Ci.columns[0]:Ci.columns[0]+str(i)}
	Ci = Ci.rename(columns=dict)
	all_col_list.append(Ci)

d = pd.concat(all_col_list, axis=1, join='inner')
sum = d.aggregate('sum', axis = 1)
e = pd.concat([d, sum], axis=1, join='inner')
e = e.rename(columns={0:'horizonsum'})
sum = e.aggregate('sum', axis = 0)
e = e.append(sum, ignore_index=True)
lst_a = data[var_a].unique()
lst_a.sort()
lst_a = lst_a.tolist()
lst_a.append('vertisum')
e[var_a] = lst_a

cols = e.columns.tolist()
cols = cols[-1:] + cols[:-1]
e = e[cols]

for i in cols:
	if i == 'horizonsum':
		continue
	if i == var_a:
		continue
	else:
		e['%'+i] = e[i]/e['horizonsum']*100

print(e)print(e)
print('\n')

import matplotlib.pyplot as plt

f = e.columns.tolist()
var = e.iloc[0:-1][f[0]]
outc = e.iloc[0:-1][f[-1]]
fig = plt.figure(figsize =(10, 5))
plt.bar(var,outc, label=var.name)
plt.xlabel(var.name)
plt.xticks(range(var.min(),var.max()+1))
plt.ylabel(outc.name)
plt.title('Chi-square p-value: '+ str(p))
plt.show()



