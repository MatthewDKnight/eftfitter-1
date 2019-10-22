#with open("equations/stage0_xs.txt", "r") as file:
#	rows = []
#	for row in file:
#		rows.append(row)

import numpy as np
import parameters_config_EFT
from parameters_config_EFT import PARAMS
import string
import re

def splitRow(string):
	"""
	Will split a string into terms that start with '+' or '-'.
	"""
	pieces = []
	i = 0
	for j in range(len(string)):
		if string[j] == "+" or string[j] == "-":
			pieces.append(string[i:j])
			i = j
	pieces.append(string[i:].strip("\n"))
	return pieces

def decomposeTerm(term):
	"""
	Takes in a term like like '-0.682 * cT ' or '+ 10.8796 * cWW * cWW '
	and seperates it into the coefficient and the corresponding EFT
	parameter(s).
	"""
	term = term.replace("*", "") #remove *

	#find where the first parameter by finding first lowercase letter in string
	for i in range(len(term)):
		if ord(term[i])>97 and ord(term[i])<122: #use ascii
			first_half = term[:i] #contains coefficient
			second_half = term[i:] #contains term(s)
			break

	#for first_half, remove spaces and then float to give coefficient
	first_half = first_half.replace(" ", "")
	coefficient = float(first_half) 
	
	#make list of EFT parameters in term
	params = second_half.split(" ", 1)
	if params[1] == "": #if only one term
		del params[1]
	else:
		params[1] = params[1].replace(" ", "")
	return coefficient, params


def new_get_x(vals=np.zeros(shape=(39, 2)).tolist(), MINDEX=0):
	#read in text file, row by r
	with open("equations/stage1_1_xs.txt", "r") as file:
		rows=[]
		for row in file:
			rows.append(row)

	#initialise empty A and B arrays
	A_matrix=np.zeros(shape=(len(PARAMS), len(rows)))
	B_matrix=np.zeros(shape=(len(PARAMS), len(PARAMS), len(rows)))
       
        #EFT_param_vector=np.array([vals[i][1] for i in range(len(vals))])

	EFT_param_list=[]
	EFT_param_names=[]
	
	#create a list of the EFT parameter names, EFT_param_names
	for i in PARAMS.keys():
		EFT_param_list.append(i.split('_'))
		EFT_param_names.append(i.split('_')[0])

	"""
	
	"""
	index1= EFT_param_names.index('cWWMinuscB')
	index2= EFT_param_names.index('cWWPluscB')

        vals[index1][1]= 0.5*(vals[index1][1]+vals[index2][1])
	vals[index2][1]= 0.5*(vals[index2][1]-vals[index1][1])
 	EFT_param_vector=np.array([vals[i][1] for i in range(len(vals))])
        
        EFT_param_names[index1]='cWW'
        EFT_param_names[index2]='cB'

	#print(EFT_param_names)	

	scaling_functions_decomposed = []
	for j in rows:
		split_row = splitRow(j)
		del split_row[0] #removes ggH:1 bit
		scaling_functions_decomposed.append(split_row)
		print(split_row)
	
	#print(scaling_functions_decomposed)

	for y in range(len(scaling_functions_decomposed)): #for each bin
		for term in scaling_functions_decomposed[y]: #for each term in equation for this bin
			#print(term)
			coeff, params = decomposeTerm(term)
			#print(term, params)
			if len(params) == 1:
				param_index = EFT_param_names.index(params[0])
				if len(EFT_param_list[param_index]) == 2: #if there is a scaling bit
					scaling = EFT_param_list[param_index][1] #what to scale the terms by "x02"
					scaling = int(scaling.strip("x"))
				else:
					scaling = 0	 
				A_matrix[param_index][y] = coeff/10**scaling
			elif len(params) == 2:
				param_index1 = EFT_param_names.index(params[0])
				if len(EFT_param_list[param_index1]) == 2: #if there is a scaling bit
					scaling1 = EFT_param_list[param_index1][1] #what to scale the terms by "x02"
                                	scaling1 = int(scaling1.strip("x"))
				else:
					scaling1 = 0
			
				param_index2 = EFT_param_names.index(params[1])
				if len(EFT_param_list[param_index2]) == 2: #if there is a scaling bit
					scaling2 = EFT_param_list[param_index2][1] #what to scale the terms by "x02"
                               		scaling2 = int(scaling2.strip("x"))
				else:
					scaling2 = 0

				B_matrix[param_index1][param_index2][y] = coeff/10**(scaling1+scaling2)
				B_matrix[param_index2][param_index1][y] = coeff/10**(scaling1+scaling2) #make it symmetric
			else:
				print("Something went wrong")	

	return A_matrix, B_matrix
				
A, B = new_get_x()		
						



	




	


	
