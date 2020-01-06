#with open("equations/stage0_xs.txt", "r") as file:
#	rows = []
#	for row in file:
#		rows.append(row)


#from __future__ import print_function
#import tensorflow as tf
#import sys

import numpy as np
import parameters_config_EFT
from parameters_config_EFT import PARAMS
import string
import re
import sys

#functions, name_ordering = initialiseScalingFunctions() #dictionary for scaling functions

#decay_names = ["hmm", "hzg", "hzz", "hbb", "hww", "htt", "hgg", "hgluglu", "hcc", "tot"]
#self.EFT = config.PARAMS

EXCLUDE_CROSS_TERMS = False

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

def initialiseScalingFunctions(PARAMS=PARAMS, filenames=["stage0_xs.txt","decay.txt"], POIs=PARAMS.keys()):
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
		#sys.path.append('~/eftfitter-1')
	        with open("/home/hep/jd2616/eftfitter-1/equations/%s"%filename, "r") as file:
        		for row in file:
        			rows.append(row)
                                                                               
        EFT_param_list=[] #contains scaling bit plus name
        EFT_param_names=[] #list of names
       
	name_ordering = []
 
        #create a list of the EFT parameter names and their scaling factors
        for i in PARAMS.keys():
		name_ordering.append(i)
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
        #index1= EFT_param_names.index('cWWMinuscB')
        #index2= EFT_param_names.index('cWWPluscB')
        
        #EFT_param_names[index1]='cWW'
	#name_ordering[index1] = 'cWW'
	#EFT_param_list[index1] = ["cWW", 0]
        #EFT_param_names[index2]='cB'
	#name_ordering[index2] = 'cB'
	#EFT_param_list[index2] = ["cB", 0]

	#print(len(EFT_param_names))
	#print(len(PARAMS))

	new_EFT_param_names = []
	new_EFT_param_list = []
	new_name_ordering = []

	#print(name_ordering)

	for j in range(len(PARAMS)):
		if EFT_param_names[j] in POIs:
			new_EFT_param_names.append(EFT_param_names[j])
			new_name_ordering.append(name_ordering[j])
		if EFT_param_list[j][0] in POIs:
			new_EFT_param_list.append(EFT_param_list[j])

	if 'cWWMinuscB' in new_EFT_param_names:
		index1= new_EFT_param_names.index('cWWMinuscB')
		new_EFT_param_names[index1]='cWW'
	        new_name_ordering[index1] = 'cWW'
		new_EFT_param_list[index1] = ["cWW", 0]


	if 'cWWPluscB' in new_EFT_param_names:
        	index2= EFT_param_names.index('cWWPluscB')
        	EFT_param_names[index2]='cB'
		name_ordering[index2] = 'cB'
		EFT_param_list[index2] = ["cB", 0]



	#print(new_name_ordering)
	#print(EFT_param_list)
	#print(name_ordering)

        for row in rows:
		#initialise empty A and B matrices
		A = np.zeros(len(POIs))
		B = np.zeros((len(POIs), len(POIs)))

		terms = splitRow(row)
		bin_name = terms[0].split(":")[0] #removes the ":1 " at end of string
		del terms[0] #remove "ggH:1" type bit from list of terms
	
		#for term in equation	
		for term in terms:                                                                                    		
                        coeff, params = decomposeTerm(term)
		                                                                                                                
                        if len(params) == 1: #if its an A term
				if params[0] in new_EFT_param_names:
                        		param_index = new_EFT_param_names.index(params[0])
                        		scaling = new_EFT_param_list[param_index][1] #what to scale the terms by "x02"
                        		#A[param_index] = coeff/10**scaling
					A[param_index]=coeff
                        elif len(params) == 2: #if its a B term
				if params[0] in new_EFT_param_names and params[1] in new_EFT_param_names:
					param_index1 = new_EFT_param_names.index(params[0])
					scaling1 = new_EFT_param_list[param_index1][1] #what to scale the terms by "x02"
				
					param_index2 = new_EFT_param_names.index(params[1])
					scaling2 = new_EFT_param_list[param_index2][1] #what to scale the terms by "x02"
					
					if EXCLUDE_CROSS_TERMS:
						if param_index1 == param_index2:
							#B[param_index1][param_index2] = coeff/10**(scaling1+scaling2)
							B[param_index1][param_index2] = coeff
					else:
						#B[param_index1][param_index2] = coeff/10**(scaling1+scaling2)
						B[param_index1][param_index2] = coeff
					B[param_index2][param_index1]= coeff	
					#B[param_index2][param_index1] = coeff/10**(scaling1+scaling2) #make it symmetric
                        else:
                        	print("Something went wrong")
		functions[bin_name] = [A, B]
		scaling_list=[i[1] for i in new_EFT_param_list]
	return functions, new_name_ordering, scaling_list






if __name__=="__main__":
	filenames = ["stage0_xs.txt", "decay.txt"]
	functions, new_name_ordering = initialiseScalingFunctions(PARAMS, filenames)
                                                               
