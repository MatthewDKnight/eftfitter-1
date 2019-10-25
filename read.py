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

	#add exception for cG and cA
	exceptions = ["cG", "cA"]
	for param in params:
		if param in exceptions:
			coefficient = coefficient * (4*np.pi)**2

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
        cWW and cB do not appear independently in the list of parameters
	but they do in the equations. 	
	"""
	index1= EFT_param_names.index('cWWMinuscB')
	index2= EFT_param_names.index('cWWPluscB')

        vals[index1][1]= 0.5*(vals[index1][1]+vals[index2][1])
	vals[index2][1]= 0.5*(vals[index2][1]-vals[index1][1])
 	EFT_param_vector=np.array([vals[i][1] for i in range(len(vals))])
        
        EFT_param_names[index1]='cWW'
        EFT_param_names[index2]='cB'

	#print(EFT_param_names)

	#split every row into the individual terms like "-1 * cG "
	scaling_functions_decomposed = []
	for j in rows:
		split_row = splitRow(j)
		print(split_row[0])
		del split_row[0] #removes ggH:1 bit
		scaling_functions_decomposed.append(split_row)
		print(split_row)
	
	#print(scaling_functions_decomposed)                                   

	#Makes the A and B matrices
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

def initialiseScalingFunctions(PARAMS=PARAMS, filenames=["stage0_xs.txt","decay.txt"]):
	"""
	PARAMS is the dictionary of EFT parameters from
	the EFT config file.
	This function will read in the equation files and 
	initialise the functions dictionary with the A and B
	matrices for each bin.
	"""
	functions = {}

	#read in text file, row by r
	rows = []
	for filename in filenames:
	        with open("equations/%s"%filename, "r") as file:
        		for row in file:
        			rows.append(row)
                                                                               
        EFT_param_list=[] #contains scaling bit plus name
        EFT_param_names=[] #list of names
        
        #create a list of the EFT parameter names and their scaling factors
        for i in PARAMS.keys():
		#try to split key into name and scaling part
		split_key = i.split("_")
		name = split_key[0]
		if len(split_key) == 2: #if there was a scaling factor
			scaling_factor = split_key[1]
			scaling_factor = int(scaling_factor.strip("x"))
		else: #if no scaling factor
			scaling_factor = 0	
		EFT_param_list.append([name, scaling_factor])
        	EFT_param_names.append(name)
                                                                               
        """
        cWW and cB do not appear independently in the list of parameters
        but they do in the equations. 	
        """
        index1= EFT_param_names.index('cWWMinuscB')
        index2= EFT_param_names.index('cWWPluscB')
        
        EFT_param_names[index1]='cWW'
        EFT_param_names[index2]='cB'
        for row in rows:
		#initialise empty A and B matrices
		A = np.zeros(len(PARAMS))
		B = np.zeros((len(PARAMS), len(PARAMS)))

		terms = splitRow(row)
		bin_name = terms[0].split(":")[0] #removes the ":1 " at end of string
		del terms[0] #remove "ggH:1" type bit from list of terms
	
		#for term in equation	
		for term in terms:                                                                                    		
                        coeff, params = decomposeTerm(term)
		                                                                                                                
                        if len(params) == 1: #if its an A term
                        	param_index = EFT_param_names.index(params[0])
                        	scaling = EFT_param_list[param_index][1] #what to scale the terms by "x02"
                        	A[param_index] = coeff/10**scaling
                        elif len(params) == 2: #if its a B term
                        	param_index1 = EFT_param_names.index(params[0])
                       		scaling1 = EFT_param_list[param_index1][1] #what to scale the terms by "x02"
                        
                        	param_index2 = EFT_param_names.index(params[1])
                       		scaling2 = EFT_param_list[param_index2][1] #what to scale the terms by "x02"
                                                                                                                            
                        	B[param_index1][param_index2] = coeff/10**(scaling1+scaling2)
                        	B[param_index2][param_index1] = coeff/10**(scaling1+scaling2) #make it symmetric
                        else:
                        	print("Something went wrong")
		functions[bin_name] = [A, B]
	return functions

if __name__=="__main__":
	A, B = new_get_x()

	filenames = ["stage0_xs.txt", "decay.txt"]
	functions = initialiseScalingFunctions(PARAMS, filenames)
                                                               
