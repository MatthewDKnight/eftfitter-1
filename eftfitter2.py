#Simply python fitting for STXS->EFT interpretation :
from scipy.optimize import minimize
from scipy import linalg
import array,numpy,sys
import numpy as np

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import matplotlib.cm as cm

#import ROOT as r

from read import initialiseScalingFunctions

VERB=True #print debugging messages
NSCANPOINTS=60 #determines granularity of the scan

import matplotlib.pyplot as plt
#import HIG_17_031 as data
import Higgs_PT as data
import parameters_config_EFT as config

class eft_fitter:
    def __init__(self, config):   # for now lets just play, user interface later
        self.data_sets = [] #list of data sets
        self.functions, self.name_ordering = initialiseScalingFunctions() #dictionary for scaling functions

        self.decay_names = ["hmm", "hzg", "hzz", "hbb", "hww", "htt", "hgg", "hgluglu", "hcc", "tot"]
        self.EFT = config.PARAMS
        self.POIs = config.MYPARAMS	

    def processDataSet(self,data_set):
        """
        data_set - data file

        Adds model to the fitter.
        Normalises the weights.
        Converts correlation dictionary (from model file) into a matrix and finds
        the error matrix.
        """

        #Normalise the weights for the bins
        for x in data_set.X.items(): #for each bin in data set
            names = x[1][0] #the list of bin names and corresponding weights
            if len(names) == 0:
                x[1][0] = names = [[1.,x[0]]]
            else:
               tsc=0 #counter for all weights
               for name in names:
                    weight = float(name[0])
                    tsc+=weight
               for i in range(len(names)): #for every weight+bin combinations
                   x[1][0][i][0] /= tsc #renormalise

        # covert the correlation dict into a flattened array
        ccorr = []
        #for every combination of bin
        for x in data_set.X.items():
          for y in data_set.X.items():
                if (x[0],y[0]) in data_set.correlation.keys():
                    rho = data_set.correlation[(x[0],y[0])]
                elif (y[0],x[0]) in data_set.correlation.keys():
                    rho = data_set.correlation[(y[0],x[0])]
                else:
                  if x[0]==y[0]:
                      rho=1.
                  else:
                      rho = 0.
                  print("WARNING - Assuming correlation in %s (no info given) rho(%s,%s) = %g"%(data_set,x[0],y[0],rho))
                ccorr.append(rho)

        data_set.correlation = array.array('d',ccorr)

        # symmetrize the errors
        error_vector = []
        for x in data_set.X.items():
         error_vector.append(x[1][3])

        # do some squarification and inverting
        data_set.nbins = len(data_set.X.items())
        v = data_set.correlation
        data_set.square_correlation = [v[i:i+data_set.nbins] for i in range(0,len(v),data_set.nbins)]
        data_set.variance = error_vector
        data_set.square_covariance = [ [ data_set.square_correlation[i][j] * (data_set.variance[i]*data_set.variance[j])\
                      for i in range(data_set.nbins)]\
                      for j in range(data_set.nbins)]

        data_set.err_mat = numpy.array(data_set.square_covariance)
        data_set.err_mat = linalg.inv(data_set.err_mat)
	
	#print(data_set.err_mat)

        self.data_sets.append(data_set)

    def getMeasurements(self,data_set_no):
      """
      Gets measured values of the bins for a particular data set.
      """
      X = self.data_sets[data_set_no].X
      return [x[1][2] for x in X.items()]

    def getPredictions(self, data_set_no, include_names=False):
        """
        Find STXS predictions given the EFT parameters from the scaling functions.
        """

	EFT_copy = self.EFT.copy()
	EFT_copy["cWW"] = [[0,0],0,0]
	EFT_copy["cB"] = [[0,0],0,0]
	cWWMinuscB = float(EFT_copy["cWWMinuscB_x02"][1]) / 100
	cWWPluscB = float(EFT_copy["cWWPluscB_x03"][1]) / 1000
	cWW = 0.5 * (cWWMinuscB + cWWPluscB)
	cB = 0.5 * (cWWPluscB - cWWMinuscB)
	EFT_copy["cWW"][1] = cWW
	EFT_copy["cB"][1] = cB
	
	EFT_vector = []
	for name in self.name_ordering:
		EFT_vector.append(EFT_copy[name][1])	

        # note that some models may use a different naming for the SCALING function string (eg stage1 vs 1.1)
        dataset = self.data_sets[data_set_no]

        for x in dataset.X.items():
            names = x[1][0]
            tsc = 0.
            for name in names:
                weight = float(name[0])
                name = name[1]
                if "R_BR" in name: # in this case, we have a ratio of ratios model, expect parameter BR_hxx_BR_hyy - THIS IS VERY SPECIFIC TO SOME MODELS (eg STXS combination)!
                    Bxx = name.split("BR_")[1]
                    Byy = name.split("BR_")[2]
                    A_vector_Bxx=self.functions[Bxx][0]
                    B_matrix_Bxx=self.functions[Bxx][1]
                    A_vector_Byy=self.functions[Byy][0]
                    B_matrix_Byy=self.functions[Byy][1]
                    nom = 1+np.dot(A_vector_Bxx, EFT_vector)+np.dot(np.transpose(EFT_vector), np.dot(B_matrix_Bxx, EFT_vector))
                    dnom = 1+np.dot(A_vector_Byy, EFT_vector)+np.dot(np.transpose(EFT_vector), np.dot(B_matrix_Byy, EFT_vector))
                    sc = nom/dnom
                else:
                    split_name=name.split('_')
                    if split_name[-1] in self.decay_names:
                        prod_name='_'.join(split_name[:-1])
                        decay_name=split_name[-1]
                        A_vector_prod=self.functions[prod_name][0]
                        B_matrix_prod=self.functions[prod_name][1]
                        A_vector_decay=self.functions[decay_name][0]
                        B_matrix_decay=self.functions[decay_name][1]
                        A_vector_decay_total=self.functions["tot"][0]
                        B_matrix_decay_total=self.functions["tot"][1]
                        prod = 1+np.dot(A_vector_prod, EFT_vector)+np.dot(np.transpose(EFT_vector), np.dot(B_matrix_prod, EFT_vector))
                        decay =(1+np.dot(A_vector_decay, EFT_vector)+np.dot(np.transpose(EFT_vector), np.dot(B_matrix_decay, EFT_vector)))
                        decay_total= (1+np.dot(A_vector_decay_total, EFT_vector)+np.dot(np.transpose(EFT_vector), np.dot(B_matrix_decay_total, EFT_vector)))
                        sc = prod*decay/decay_total

                    elif name in self.functions.keys():
                        A_vector=self.functions[name][0]
                        B_matrix=self.functions[name][1]
                        sc= 1+np.dot(A_vector, EFT_vector)+np.dot(np.transpose(EFT_vector), np.dot(B_matrix, EFT_vector))

                    else:
                        print("Could not extract production/decay names")

                tsc+=weight*sc
            dataset.X[x[0]][1]=tsc
        if include_names:
            #return [(x[0]+"_"+dataset.decay,x[1][1]) for x in dataset.X.items()]
            return [(x[0],x[1][1]) for x in dataset.X.items()]
        else:
            return [x[1][1] for x in dataset.X.items()]

    def calculate_x(self):
      for i in range(len(self.data_sets)):
          self.getPredictions(i)

    def reset(self):
      # finds predictions given the resettable EFT parameter values
      for e in self.EFT.keys():
          self.EFT[e][1] = self.EFT[e][2]
      self.calculate_x()

    def neg_log_likelihood(self,x,*args):
      eft_keys = args[0]

      for i in range(len(eft_keys)):
          self.EFT[eft_keys[i]][1]=x[i]

      constr=0

      for i, M in enumerate(self.data_sets):
        predictions = self.getPredictions(i)
        measurements = self.getMeasurements(i)

        xarr  = numpy.array([xx-xx0 for xx,xx0 in zip(predictions,measurements)])
        xarrT = xarr.T

        constr += 0.5*(xarrT.dot(M.err_mat.dot(xarr)))

      return constr

    def minimizer(self,rv={},constrained=False,params_list=[]):  # params_list is now list of POI
        if constrained:
            for i in range(len(params_list)):
                 self.EFT[params_list[i]][1]=rv[params_list[i]]

	POI_dict = dict(E for E in filter(lambda x: x[0] in fitter.POIs, self.EFT.items()))
	POI_dict = dict(E for E in filter(lambda x: x[0] not in params_list, POI_dict.items())) #takes out scanning param

	#init_CFG = [[e[0],float(e[1][0][1])] for e in POI_dict.items()]

        init_CFG = [[e[0],float(e[1][1])] for e in POI_dict.items()]
        eft_keys = [i[0] for i in init_CFG]
        init = [i[1] for i in init_CFG]

        bounds = [(self.EFT[v][0][0],self.EFT[v][0][1]) for v in eft_keys]
        xbest = minimize(self.neg_log_likelihood,init,eft_keys,bounds=bounds)
	results = [[eft_keys[i],xbest.x[i]] for i in range(len(eft_keys))]
	
	chi2 = 2 * xbest.fun	

	return results, chi2

    def scan_LH(self,param, R,do_profile=True):
        """
        param - EFT parameter to make scan for
        R - list of sampling points
        Make a 1D scan for a particular EFT parameter.
        """

        # make a 1D scan of a particular EFT parameter, choose whether to profile remaining parameters or leave at 0

        self.reset() #reset all predictions

        if do_profile:
            print("Profiling %s in range [%g,%g]"%(param,R[0],R[-1]))
        else:
            print("Scanning %s in range [%g,%g]"%(param,R[0],R[-1]))

        # in the scan, keep track of the scaling functions ...
        floated_POI_values = [] #list of floated EFT params
        scaling_functions = []
        LH = []
        minll = 9999
        for r in R:
            #self.reset() #reset all predictions
            if do_profile:
                res = self.minimizer(rv={param:r},constrained=True,params_list=[param])
                if res[1] < minll:
                    minll = res[1]
                LH.append(res[1])
                floated_POI_values.append(res[0])
            else:
                res = 2*self.neg_log_likelihood([r],[param])
                if res< minll:
                    minll = res
                LH.append(res)
             # now for every process, get the value of it
            self.reset() #reset all predictions
            mus = []
            for data_set_no in range(len(self.data_sets)):
                self.EFT[param][1] = r
                predictions = self.getPredictions(data_set_no, True)
		#predictions = self.get_x(data_set_no, True)
                for prediction in predictions:
                    mus.append(prediction)
            scaling_functions.append(mus)

        LH = [lh-minll for lh in LH]
        return LH,floated_POI_values,scaling_functions

    def scan(self,param):

        pv = self.EFT[param] #get min, max, current and nominal value for this EFT param
        np = NSCANPOINTS #number of points to sample when scanning
        R = numpy.linspace(pv[0][0],pv[0][1],np) #array of points between min and max

        C_RES_PROF  = self.scan_LH(param,R,1)
        C_RES_FIXED = self.scan_LH(param,R,0)

        C_prof  = C_RES_PROF[0]
        C_fixed = C_RES_FIXED[0]

        profiled_POIs = C_RES_PROF[1]
        scaling_functions = C_RES_FIXED[2]

        fig, ax1 = plt.subplots()
        ax1.plot(R,C_prof,color='black',linewidth=3,linestyle='-',label="Profiled")
        ax1.plot(R,C_fixed,color='black',linewidth=3,linestyle='--',label="Scan")

        ax1.set_ylabel("$\Delta\chi^{2}$",fontsize=20)
        ax1.set_xlabel("%s"%param,fontsize=20)
        plt.ylim(0,10)

        if len(profiled_POIs[0]):
            ax2 = ax1.twinx()
            poilabels = []
            for P in [p[0] for p in profiled_POIs[0]]:
                poilabels.append(P)
            for i,P in enumerate(poilabels):
                vals = [p[i][1] for p in profiled_POIs]
                ax2.plot(R,vals, label=P)

            ax2.set_ylabel("Profiled EFT coeff.")
            ax2.legend(fontsize=9,loc='upper right')

        ax1.axvline(0., linestyle='--', color='k')
        ax1.axhline(1., linestyle='--', color='r') # horizontal lines
        ax1.axhline(4., linestyle='--', color='r') # horizontal lines
        ax1.legend(fontsize=9,loc='upper left')

        plt.savefig("eftfitter2_Graphs/%s.pdf"%(param))
        plt.savefig("eftfitter2_Graphs/%s.png"%(param))

        plt.clf()
        plt.cla()
        plt.close()

        fig, ax1 = plt.subplots()
        styles = [".","v","+","o","*","s","x","p"]
        if len(scaling_functions):
            scalers = []
            for P in [ p[0] for p in scaling_functions[0] ]:
                scalers.append(P)
            for i,P in enumerate(scalers):
                vals = [p[i][1] for p in scaling_functions]
                ax1.plot(R,vals,styles[i//7], label=P)

            ax1.set_ylabel("$\mu$",fontsize=20)
            ax1.set_xlabel("%s"%param,fontsize=20)
            box = ax1.get_position()
            ax1.set_position([box.x0*0.8, box.y0, box.width * 0.7, box.height])

            # Put a legend to the right of the current axis
            ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=6,ncol=2)
            plt.savefig("eftfitter2_Graphs/stxs_scaling_vs_%s.pdf"%(param))
            plt.savefig("eftfitter2_Graphs/stxs_scaling_vs_%s.png"%(param))

            plt.clf()
            plt.cla()
            plt.close()

if __name__=="__main__":
    fitter = eft_fitter(config)
    fitter.processDataSet(data)

    for e in config.MYPARAMS:
        fitter.scan(e)
