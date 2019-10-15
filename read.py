with open("equations/stage0_xs.txt", "r") as file:
	rows = []
	for row in file:
		rows.append(row)

import numpy
import parameters_config_EFT
from parameters_config_EFT import PARAMS
A_matrix=numpy.zeros(shape=(len(PARAMS), len(PARAMS)))
B_matrix=numpy.zeros(shape=(len(PARAMS), len(PARAMS), len(PARAMS)))

EFT_param_list=[]
EFT_param_names=[]

for i in PARAMS.keys():
	EFT_param_list.append(i.split('_'))
	EFT_param_names.append(i.split('_')[0])


scaling_functions_decomposed=[]
for j in rows:
	scaling_functions_decomposed.append(j.split('+','-'))




	


	
