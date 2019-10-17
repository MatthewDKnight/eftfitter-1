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
	#find where the first term starts by finding lowercase letter
	term = term.replace("*", "") #remove *

	for i in range(len(term)):
		if ord(term[i])>97 and ord(term[i])<122: #use ascii
			first_half = term[:i] #contains coefficient
			second_half = term[i:] #contains term(s)
			break
	first_half = first_half.replace(" ", "") #remove spaces
	coefficient = float(first_half) 
	
	terms = second_half.split(" ", 1)
	
	if terms[1] == "": #if only one term
		del terms[1]
	else:
		terms[1] = terms[1].replace(" ", "")
	return coefficient, terms 


def new_get_x(vals=np.zeros(shape=(39, 2)).tolist(), MINDEX=0):
	with open("equations/stage1_1_xs.txt", "r") as file:
		rows=[]
		for row in file:
			rows.append(row)
	A_matrix=np.zeros(shape=(len(PARAMS), len(rows)))
	B_matrix=np.zeros(shape=(len(PARAMS), len(PARAMS), len(rows)))
        #EFT_param_vector=np.array([vals[i][1] for i in range(len(vals))])

	EFT_param_list=[]
	EFT_param_names=[]

	for i in PARAMS.keys():
		EFT_param_list.append(i.split('_'))
		EFT_param_names.append(i.split('_')[0])
	index1= EFT_param_names.index('cWWMinuscB')
	index2= EFT_param_names.index('cWWPluscB')

        vals[index1][1]= 0.5*(vals[index1][1]+vals[index2][1])
	vals[index2][1]= 0.5*(vals[index2][1]-vals[index1][1])
 	EFT_param_vector=np.array([vals[i][1] for i in range(len(vals))])
        
        EFT_param_names[index1]='cWW'
        EFT_param_names[index2]='cB'

	print(EFT_param_names)
	
	"""
	#print (EFT_param_names)
	scaling_functions_decomposed=[]
	for j in rows:
		row_decomposed=[]
		totally_split_row=list(j)
		unsplit_terms=[]
		for k in range(len(totally_split_row)):
			if totally_split_row[k]=='+':
				row_decomposed.append("".join(unsplit_terms))
				unsplit_terms=['+']
               		elif totally_split_row[k]=='-':
				row_decomposed.append("".join(unsplit_terms))
				unsplit_terms=['-']
			else:
				unsplit_terms.append(totally_split_row[k])
		scaling_functions_decomposed.append(row_decomposed)
	"""
	scaling_functions_decomposed = []
	for j in rows:
		split_row = splitRow(j)
		del split_row[0] #removes ggH:1 bit
		scaling_functions_decomposed.append(split_row)
		print(j)
		print(split_row)

	"""	
	for y in range(len(scaling_functions_decomposed)):
		for term in scaling_functions_decomposed[y]:
			print('term; '+term)
			for z in range(len(EFT_param_names)):
				if ' '+EFT_param_names[z]+' ' in term:
					#print('EFT_param; '+EFT_param_names[z])
					if term.count('*') == 1:
						#print('EFT_param; '+EFT_param_names[z])
					#if not any(x in term.replace(EFT_param_names[z], '', 1) for x in EFT_param_names):
						if '+ ' in term:
							#print('1Dnew_entry; '+term.replace(' * ' +EFT_param_names[z], '', 1).replace('+ ', '', 1))
							new_entry=float(term.replace( ' * '+EFT_param_names[z], '', 1).replace('+ ', '', 1))
						elif '- ' in term:
							#print('1Dnew_entry; '+term.replace(' * '+EFT_param_names[z], '', 1).replace('- ', '', 1))
							new_entry=-float(term.replace(' * '+EFT_param_names[z], '', 1).replace('- ', '', 1))
						else:
							#print('1Dnew_entry; '+term.replace(' * '+EFT_param_names[z], '', 1))
							new_entry=float(term.replace(' * '+EFT_param_names[z], '', 1)) 
						A_matrix[z][y]=new_entry
					elif term.count('*') == 2:
						EFT_param_1=EFT_param_names[z]
						#print('EFT_param_1; '+EFT_param_1)
						term=term.replace(' * '+EFT_param_names[z], '', 1)
						for a in range(len(EFT_param_names)):
							if ' '+EFT_param_names[a]+' ' in term:
								EFT_param_2=EFT_param_names[a]
								#print('EFT_param_2; '+EFT_param_2)
								if '+ ' in term:
									#print('2Dnew_entry; '+term.replace(' * '+EFT_param_2, '', 1).replace('+ ', ''))
									new_entry=float(term.replace(' * '+EFT_param_2, '', 1).replace('+ ', ''))
								elif '- ' in term:
									#print('2Dnew_entry; '+term.replace(' * '+EFT_param_2, '', 1).replace('- ', ''))
									new_entry=-float(term.replace(' * '+EFT_param_2, '', 1).replace('- ', ''))
								else:
									#print('2Dnew_entry; '+term.replace(' * '+EFT_param_2, '', 1))
									new_entry=float(term.replace(' * '+EFT_param_2, '', 1))
								B_matrix[z][a][y]=new_entry
	"""					
	
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
						



	




	


	
