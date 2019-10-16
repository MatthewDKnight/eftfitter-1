#with open("equations/stage0_xs.txt", "r") as file:
#	rows = []
#	for row in file:
#		rows.append(row)

import numpy as np
import parameters_config_EFT
from parameters_config_EFT import PARAMS
import string
import re


def new_get_x(vals=np.zeros(shape=(39, 2)).tolist(), MINDEX=0):
	with open("equations/stage0_xs.txt", "r") as file:
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
	
	for y in range(len(scaling_functions_decomposed)):
		for term in scaling_functions_decomposed[y]:
			print('term; '+term)
			for z in range(len(EFT_param_names)):
				if ' '+EFT_param_names[z]+' ' in term:
					#print('EFT_param; '+EFT_param_names[z])
					if term.count('*') == 1:
						print('EFT_param; '+EFT_param_names[z])
					#if not any(x in term.replace(EFT_param_names[z], '', 1) for x in EFT_param_names):
						if '+ ' in term:
							print('1Dnew_entry; '+term.replace(' * ' +EFT_param_names[z], '', 1).replace('+ ', '', 1))
							new_entry=float(term.replace( ' * '+EFT_param_names[z], '', 1).replace('+ ', '', 1))
						elif '- ' in term:
							print('1Dnew_entry; '+term.replace(' * '+EFT_param_names[z], '', 1).replace('- ', '', 1))
							new_entry=-float(term.replace(' * '+EFT_param_names[z], '', 1).replace('- ', '', 1))
						else:
							print('1Dnew_entry; '+term.replace(' * '+EFT_param_names[z], '', 1))
							new_entry=float(term.replace(' * '+EFT_param_names[z], '', 1)) 
						A_matrix[z][y]=new_entry
					elif term.count('*') == 2:
						EFT_param_1=EFT_param_names[z]
						print('EFT_param_1; '+EFT_param_1)
						term=term.replace(' * '+EFT_param_names[z], '', 1)
						for a in range(len(EFT_param_names)):
							if ' '+EFT_param_names[a]+' ' in term:
								EFT_param_2=EFT_param_names[a]
								print('EFT_param_2; '+EFT_param_2)
								if '+ ' in term:
									print('2Dnew_entry; '+term.replace(' * '+EFT_param_2, '', 1).replace('+ ', ''))
									new_entry=float(term.replace(' * '+EFT_param_2, '', 1).replace('+ ', ''))
								elif '- ' in term:
									print('2Dnew_entry; '+term.replace(' * '+EFT_param_2, '', 1).replace('- ', ''))
									new_entry=-float(term.replace(' * '+EFT_param_2, '', 1).replace('- ', ''))
								else:
									print('2Dnew_entry; '+term.replace(' * '+EFT_param_2, '', 1))
									new_entry=float(term.replace(' * '+EFT_param_2, '', 1))
								B_matrix[z][a][y]=new_entry
						
	

	print([A_matrix[i][0] for i in range(len(A_matrix))])
						
						



	




	


	
