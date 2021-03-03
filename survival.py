#Survival analysis

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import shapiro
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

df = soapsheetin(path)


import sksurv
from sksurv.nonparametric import kaplan_meier_estimator
import lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

def single_kmc(data, status, interval):

	sta = []
	for i in data[status]:
    	b = bool(i)
    	sta.append(b)

	interv = data[interval]

	time, survival_prob = kaplan_meier_estimator(sta, interv)
	plt.step(time, survival_prob)
	plt.ylabel("est. probability of survival")
	plt.xlabel("time $(days)$")



def compare_kmc(data, factor, status, interval):
	f = data[factor].drop_duplicates().tolist()
	f.sort()
	
	for i in f:
		group_i = data[data[factor] == i]
		interv_i = group_i[interval]
		interv_i.reset_index(drop=True)
		interv_i = interv_i.tolist()
		interv_i.append(0)
	
		sta_i = []
		for j in group_i[status]:
			c = bool(j)
			sta_i.append(c)
		sta_i.append(False)

	
	

		time, survival_prob = kaplan_meier_estimator(sta_i, interv_i)
		plt.step(time, survival_prob, where='post', label=str(factor)+'=%s' % i)
	plt.ylabel("est. probability of survival")
	plt.xlabel("time $(days)$")
	plt.legend(loc='best')

	if data[factor].nunique() != 2:
		print('The factor', factor, 'is non-binary, so the Logrank statistics is not calculated')
	else:
		f = data[factor].drop_duplicates().tolist()
		f.sort()
		time = []
		censor = []
		
		for i in f:
			interv_i = []
			group_i = df[df[factor] == i]
			for j in group_i[interval]:
				interv_i.append(j)
			time.append(interv_i)
				
			censor_i = []
			for k in group_i[status]:
				censor_i.append(k)
			censor.append(censor_i)
		T = time[0]
		T1 = time[1]
		E = censor[0]
		E1 = censor[1]

		results = logrank_test(T, T1, E, E1)
		results.print_summary()


			




	
