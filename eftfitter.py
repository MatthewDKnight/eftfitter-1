# Simply python fitting for STXS->EFT interpretation : 
from scipy.optimize import minimize
from scipy import linalg 
import array,numpy,sys
from matplotlib import pyplot as plt
import matplotlib.cm as cm 

import ROOT as r 

VERB=False


class eft_fitter:
  def __init__(self, EFT_PARAMETERS):   # for now lets just play, user interface later
    
    self.EFT_PARAMETERS = EFT_PARAMETERS 

    fw = r.TFile.Open("result.root")    
    self.w = fw.Get("w")    
  
    self.EFT = {
    "clW_x02"           :[[-5 , 5    ]		,0,0]    # current value and nominal / resettable value 
    ,"cHu_x02"          :[[-1.1 , 1.1]		,0,0]
    ,"c2W"              :[[-5 , 5  	]	,0,0]
    ,"cl"               :[[-5 , 5  	]	,0,0]
    ,"cdG_x02"          :[[-5 , 5  	]	,0,0]
    ,"cH_x01"           :[[-1.4 , 1.94]		,0,0]
    ,"cpHL_x02"         :[[-5 , 5]		,0,0]
    ,"c2B"              :[[-5 , 5 ] 		,0,0]
    ,"cG_x04"           :[[-50. , 40.]		,0,0]
    ,"tcA_x04"          :[[-12 , 12 	]	,0,0]
    ,"cT_x03"           :[[-4.3 , 3.3 ] 	,0,0]
    ,"tc3W_x01"         :[[-1.8 , 1.8  ]	,0,0]
    ,"cWWPluscB_x03"    :[[-3.3 , 1.8  ]	,0,0]
    ,"cpHQ_x03"         :[[-4.4 , 4.4  ]	,0,0]
    ,"cHud_x02"         :[[-5 , 5  	]	,0,0]
    ,"cHe_x03"          :[[-1.8 , 0.25] 	,0,0]
    ,"cA_x04"           :[[-10000 , 80000 ]		,0,0]
    ,"cWWMinuscB_x03"   :[[-35 , 150  	]	,0,0]
    ,"tcHB_x01"         :[[-2.4 , 2.4]		,0,0]
    ,"cHQ_x03"          :[[-1.9 , 6.9 ] 	,0,0]
    ,"c3W_x02"          :[[-8.3 , 4.5  ]	,0,0]
    ,"cuB_x02"          :[[-5 , 5 	]	,0,0]
    ,"c2G_x04"          :[[-1.6 , 1.6 ] 	,0,0]
    ,"cu_x02"           :[[-40. , 20. ]		,0,0]
    ,"cHB_x02"          :[[-4.5 , 7.5  ]	,0,0]
    ,"c3G_x04"          :[[-1.6 , 1.6  ]	,0,0]
    ,"cdW_x02"          :[[-5 , 5  	]	,0,0]
    ,"cHW_x02"          :[[-20. , 30.]		,0,0]
    ,"c6"               :[[-5 , 5  	]	,0,0]
    ,"tcHW_x02"         :[[-6 , 6  	]	,0,0]
    ,"tcG_x04"          :[[-1.2 , 1.2]		,0,0]
    ,"cHL_x02"          :[[-5 , 5  	]	,0,0]
    ,"cdB_x02"          :[[-5 , 5  	]	,0,0]
    ,"cuW_x02"          :[[-5 , 5  	]	,0,0]
    ,"cHd_x02"          :[[-4.2 , 0.44] 	,0,0]
    ,"cd_x02"           :[[-19.8 , 8.8 ]	,0,0]
    ,"clB_x02"          :[[-5 , 5  	]	,0,0]
    ,"cuG_x02"          :[[-5 , 5	]	,0,0]
    ,"tc3G_x04"         :[[-1.6 , 1.6]		,0,0]
    }
    
    self.MODELS = []

  def processModel(self,model,decay):

   # make weights sum to 1
   for x in model.X.items():
    names = x[1][0] 
    if len(names) == 0: 
       x[1][0] = names = [[1.,x[0]]]
    elif len(names) == 0: 
       x[1][0][0] = 1. 
    else: 
     tsc=0
     for name in names: 
	  weight = float(name[0])
	  tsc+=weight
     for i in range(len(names)): 
	  x[1][0][i][0] /= tsc # renormalise


   # covert the correlation dict into a matrix 
   ccorr = []
   for x in model.X.items():
    for y in model.X.items(): 
      if (x[0],y[0]) in model.correlation.keys() :
         rho = model.correlation[(x[0],y[0])]
      else: rho = model.correlation[(y[0],x[0])] 
      #print x[0],y[0],rho 
      ccorr.append(rho)

   model.correlation = array.array('d',ccorr)
   model.decay = decay 

   # symmetrize the errors 
   error_vector = []  
   for x in model.X.items(): error_vector.append(x[1][3])

   # do some squarification and inverting 
   model.nbins = len(model.X.items())
   v = model.correlation
   model.square_correlation = [v[i:i+model.nbins] for i in range(0,len(v),model.nbins)]
   model.variance = error_vector
   model.square_covariance = [ [ model.square_correlation[i][j] * (model.variance[i]*model.variance[j])\
				    for i in range(model.nbins)]\
				    for j in range(model.nbins)]

   model.err_mat = numpy.array(model.square_covariance)
   model.err_mat = linalg.inv(model.err_mat)
   # finally lets make all of the names extend with _decay 
   self.MODELS.append(model)


  def print_EFT(self):
    print " ---- EFT Parameter Values ---- " 
    for eft in self.EFT.items():
      print "%s = %g"%(eft[0],eft[1][1])
    print " ------------------------------ " 

  def print_X(self):
    print " ---- STXS Parameter Values ---- " 
    for i,M in enumerate(self.MODELS): 
      X = M.X
      print "Data set %d"%i
      for x in X.items():
       print "%s = %g +/- %g (measured), %g (predicted at EFT values) "%(x[0]+"_"+M.decay,x[1][2],x[1][3],x[1][1])
    print " ------------------------------ " 

  def get_x0(self,MINDEX):
   
    X = self.MODELS[MINDEX].X
    return [x[1][2] for x in X.items()]

  ############## Dummy Function to test Gaussian Constraint! #####################
  """
  def calculate_x(self,vals): 
    #print " I'm being asked to set the following " , vals
    for v in vals:
      self.EFT[v[0]][1] = v[1]
      self.w.var(v[0]).setVal(v[1])      

    for x in self.X.items():
      tsc=0.
      if x[0]=='ggH_0J': tsc=self.w.var('cG_x04').getVal()
      elif x[0]=='ggH_1J_low': tsc=self.w.var('cHW_x02').getVal()
      elif x[0]=='ggH_1J_med': tsc=self.w.var('cWWMinuscB_x03').getVal()
      self.X[x[0]][2]=tsc+1.
    return [x[1][2] for x in self.X.items()]
  """
  
  def get_x(self,vals,MINDEX,include_names=False): 
    #if VERB: print " Setting following (to recalculate STXS bins) --> " , vals
    for v in vals:
      self.EFT[v[0]][1] = v[1]
      self.w.var(v[0]).setVal(v[1])      

    model = self.MODELS[MINDEX] 
    for x in model.X.items():
      names = x[1][0] 
      if not len(names) : 
      	names = [[1.,x[0]]]
      tsc = 0.
      #print "New list -> ", x[0],names 
      for name in names: 
        weight = float(name[0])
	name = name[1]
	if "BR" in name: # in this case, we have a ratio of ratios model, expect parameter BR_hxx_BR_hyy 
	  Bxx = name.split("BR_")[1] 
	  Byy = name.split("BR_")[2] 
	  nom  = self.w.function("scaling_%s"%(Bxx)).getVal(r.RooArgSet())
	  dnom = self.w.function("scaling_%s"%(Byy)).getVal(r.RooArgSet())
	  sc = nom/dnom
	else: sc = self.w.function("stxs1toeft_scaling_%s_%s_13TeV"%(name,model.decay)).getVal(r.RooArgSet())
	tsc+=weight*sc 
      model.X[x[0]][1]=tsc
    #if VERB: 
    #  self.print_EFT()
    #  self.print_X()
    if include_names: return [(x[0]+"_"+model.decay,x[1][1]) for x in model.X.items()]
    else : return [x[1][1] for x in model.X.items()]
  
  def calculate_x(self,vals): 
    if VERB: print " Setting following (to recalculate STXS bins) --> " , vals
    for i in range(len(self.MODELS)): self.get_x(vals,i)
    if VERB:
     self.print_EFT()
     self.print_X()

  def neg_log_likelihood(self,ECFG,*args):
    #print " my current EFT ", E
    args= args[0]
    #print ECFG
    E = [ [i,e] for i,e in zip(args['eft_keys'],ECFG)]

    constr=0
    for i, M in enumerate(self.MODELS):
      x  = self.get_x(E,i)
      x0 = self.get_x0(i)

      xarr  = numpy.array([xx-xx0 for xx,xx0 in zip(x,x0)])
      xarrT = xarr.T

      constr += 0.5*(xarrT.dot(M.err_mat.dot(xarr)))
    
    return constr
 
  def minimizer(self,rv=0,constrained=False,params_list=[]):  # params_list is now list of POI

   if constrained:
     self.EFT[params_list[0]][1]=rv
     E=[[e[0],float(e[1][1])] for e in self.EFT.items()]
     self.calculate_x(E)

   self.EFT_safe = self.EFT.copy()

   #if constrained and params_list: 
   self.EFT = dict(E for E in filter(lambda x: x[0] not in params_list, self.EFT.items()))

   init_CFG = [[e[0],float(e[1][1])] for e in self.EFT.items()]
   eft_keys = {"eft_keys":[i[0] for i in init_CFG]}
   init = [i[1] for i in init_CFG]

   if VERB: print "My EFT parameter llist is --> ", params_list
   results = []
   if params_list: results = [params_list[0],rv]

   if len(init) :
        bounds = [(self.EFT[v][0][0],self.EFT[v][0][1]) for v in eft_keys['eft_keys']]
   	xbest = minimize(self.neg_log_likelihood,init,eft_keys,bounds=bounds)
	results = [[e[0],i] for e,i in zip(self.EFT.items(),xbest.x)]
   
   self.EFT = self.EFT_safe.copy()
   
   if VERB: 
     print "Finished a minimization ... " 
     print " ... Results = ", results
   
   self.calculate_x(results)
   
   if VERB: 
    self.print_EFT()
    self.print_X()

   return results, 2*self.neg_log_likelihood([r[1] for r in results],eft_keys)

  def prep(self):
    # 1. Set all the things to 0 
    #print self.calculate_x([0 for i in self.EFT.items()])
    # 2. remove the useless parameters from the list (user asks for only some of them anyway)

    fi = r.TFile("inputs_converted.root","RECREATE")
    for i,M in enumerate(self.MODELS): 
     d = fi.mkdir("DataSet_%d"%i)
     # make a weird ROOT file of the results why not ?
     hcorr = r.TH2F("h2corr","Correlations",M.nbins,0,M.nbins,M.nbins,0,M.nbins)
     hcov  = r.TH2F("h2cov","Covariance",M.nbins,0,M.nbins,M.nbins,0,M.nbins)
     hgr   = r.TH1F("hgr","Fitted values",M.nbins,0,M.nbins)
     hgr.SetMarkerStyle(20); hgr.SetMarkerSize(1.0); hgr.SetLineWidth(3)
     for i,x in enumerate(M.X.items()): 
      hgr.SetBinContent(i+1,x[1][1])
      hgr.SetBinError(i+1,M.variance[i])
      hgr.GetXaxis().SetBinLabel(i+1,x[0])
      hcorr.GetXaxis().SetBinLabel(i+1,x[0])
      hcov.GetXaxis().SetBinLabel(i+1,x[0])
      for j,y in enumerate(M.X.items()):
        hcorr.GetYaxis().SetBinLabel(j+1,y[0])
        hcorr.SetBinContent(i+1,j+1,M.square_correlation[i][j])
        hcov.GetYaxis().SetBinLabel(j+1,y[0])
        hcov.SetBinContent(i+1,j+1,M.square_covariance[i][j])
     hgr.GetYaxis().SetTitle("#mu #pm #sigma") 
     d.cd()
     hcorr.Write()
     hcov.Write()
     hgr.Write() 

    self.EFT = dict(E for E in filter(lambda x: x[0] in self.EFT_PARAMETERS, self.EFT.items()))

    for v in self.EFT.keys(): 
      xmin,xmax = self.EFT[v][0][0],self.EFT[v][0][1]
      self.w.var(v).setMin(xmin)
      self.w.var(v).setMax(xmax)

    if VERB: 
     print "------- Setup state ------>"   
     self.print_X()
     self.print_EFT()
     print "-------------------------->"   

   
  def global_fit(self): 
    best_fit,nll2 = self.minimizer()
    self.calculate_x(best_fit)      # Note that this also sets the values of the EFT vector to the ones from the fit!
    for e in best_fit : self.EFT[e[0]][2]=e[1]
    self.print_EFT()

  def reset(self):
    # resets EFT parameters to 0 
    self.calculate_x([[e,self.EFT[e][2]] for e in self.EFT.keys()])

  def scan2d(self, px, py): # set do_profile off here!
    # make a 2D scan of a likelihood, don't profile other things !    
    self.reset()
    np = 50 
    pxx = self.EFT[px]
    pyy = self.EFT[py]

    xx = numpy.linspace(pxx[0][0],pxx[0][1],np)
    yy = numpy.linspace(pyy[0][0],pyy[0][1],np)

    C = []
    minll = 99999999
    for i in range(np):
      cc = []
      for j in range(np):
        nll2 = 2*self.neg_log_likelihood([xx[i],yy[j]],{'eft_keys':[px,py]})
	if nll2<minll: minll = nll2
        cc.append(nll2)
	#print " -> ", xx[i],yy[j], nll2
      C.append(cc)
    C = numpy.array(C) 
    #print " ----> ? ", 0, 0, 2*self.neg_log_likelihood([0,0],{'eft_keys':[px,py]})
    # always start and end with a reset in any scan
    for c in range(len(C)): 
     for i in range(len(cc)): C[c][i] -= minll
    self.reset()


    a2D = plt.subplot(111)
    conts = plt.contour(yy,xx,C,levels=[2.3,5.99], colors='b')  # the way I constructed C, y, is the faster variable (think like a matrix)
    plt.clabel(conts, fontsize=9, inline=1)
    plt.contourf(yy,xx,C,levels=numpy.arange(0,6,0.2),cmap=cm.gray)  # the way I constructed C, y, is the faster variable (think like a matrix)
    plt.colorbar()
    a2D.set_ylabel(px)
    a2D.set_xlabel(py)
    a2D.axhline(0., linestyle='--', color='k') # horizontal lines
    a2D.axvline(0., linestyle='--', color='k') # vertical lines

    plt.savefig("scan_2d_%s_%s.pdf"%(px,py));
    plt.savefig("scan_2d_%s_%s.png"%(px,py));
    
    plt.clf()
    plt.cla()
    plt.close()
        
    
  
  def scan_LH(self,param, R,do_profile=True): 
    # make a 1D scan of a particular EFT parameter, choose whether to profile remaining parameters or leave at 0 
    self.reset()
    if param not in self.EFT.keys() :
     print "No EFT parameter named: %s"%param
     exit()

    pv = self.EFT[param]
    print "Scanning %s in range [%g,%g]"%(param,pv[0][0],pv[0][1])

    #if not do_profile: params_list = self.EFT.keys()
    #else: params_list = [param]
    #if do_profile : print "List of profiled EFT params ? -> ", filter(lambda x: x not in params_list, self.EFT.keys())
    #else: print "Fixing all other EFT parameters in scan -> " ,filter(lambda x: x not in [param], self.EFT.keys())
    
    # in the scan, keep track of the scaling functions ... 
    scalers = []
    proc_scalers = []
    C = []
    minll = 9999
    for r in R : 
      if do_profile : 
        res = self.minimizer(rv=r,constrained=True,params_list=[param])
	if res[1] < minll : minll = res[1]
        C.append(res[1])
	scalers.append(res[0])
      else: 
        res = 2*self.neg_log_likelihood([r],{'eft_keys':[param]}) 
	if res< minll: minll = res
      	C.append(res)
      # now for every process, get the value of it 
      pscaler = []
      for MINDEX in range(len(self.MODELS)): 
        items = self.get_x([[param,r]],MINDEX,True)
	for item in items: pscaler.append(item)
      proc_scalers.append(pscaler)
      if VERB : 
        self.calculate_x([ [e[0],e[1][1]] for e in self.EFT.items() ])
	self.print_EFT()
        self.print_X()
    C = [c-minll for c in C]
    self.reset()
    return C,scalers,proc_scalers

  def scan(self,param): 
    

    """
      for x in self.X.items(): 
	scaler.append(x[1][2]) 
      scalers.append(scaler)
    """
    pv = self.EFT[param]
    np = 40
    R = numpy.linspace(pv[0][0],pv[0][1],np)

    C_RES_PROF  = self.scan_LH(param,R,1)
    C_RES_FIXED = self.scan_LH(param,R,0)

    C_prof  = C_RES_PROF[0]
    C_fixed = C_RES_FIXED[0]

    profiled_POIs = C_RES_PROF[1]
    scaling_functions = C_RES_FIXED[2]
    
    fig, ax1 = plt.subplots()
    ax1.plot(R,C_prof,color='black',linewidth=3,linestyle='-',label="Profiled")
    ax1.plot(R,C_fixed,color='black',linewidth=3,linestyle='--',label="Scan")

    ax1.set_ylabel("$\chi^{2}$",fontsize=20)
    ax1.set_xlabel("%s"%param,fontsize=20)
    #plt.show()
    if param in ["cG_x04","cHW_x02"]: 
      plt.ylim(0,10)
    
    if len(profiled_POIs[0]):
      ax2 = ax1.twinx()
      poilabels = []
      for P in [ p[0] for p in profiled_POIs[0] ]: poilabels.append(P)
      for i,P in enumerate(poilabels):
        vals = [p[i][1] for p in profiled_POIs]
        ax2.plot(R,vals, label=P)
        
        
      
      #for i,x in enumerate(self.X.items()): ax2.plot(R,[scalers[j][i] for j in range(len(R))], label=x[0])
      ax2.set_ylabel("Profiled EFT coeff.")
      ax2.legend(fontsize=9,loc=0)

    ax1.axvline(0., linestyle='--', color='k') # horizontal lines
    ax1.legend(fontsize=9,loc=1)
    plt.savefig("%s.pdf"%(param))
    plt.savefig("%s.png"%(param))
     
    plt.clf()
    plt.cla()
    plt.close()

    fig, ax1 = plt.subplots()
    styles = [".","v","+","o","*"]
    if len(scaling_functions): 
      scalers = []
      for P in [ p[0] for p in scaling_functions[0] ]: scalers.append(P)
      for i,P in enumerate(scalers):
        vals = [p[i][1] for p in scaling_functions]
        ax1.plot(R,vals,styles[i//7], label=P)

      ax1.set_ylabel("$\mu$",fontsize=20)
      ax1.set_xlabel("%s"%param,fontsize=20)
      box = ax1.get_position()
      ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])

      # Put a legend to the right of the current axis
      ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=7)
      plt.savefig("stxs_scaling_vs_%s.pdf"%(param))
      plt.savefig("stxs_scaling_vs_%s.png"%(param))
      
      plt.clf()
      plt.cla()
      plt.close()


################### Import datasets ##############

import ATLAS36 as model 


################# Pick EFT parameters to care about and make the fitter
EFT_PARAMETERS = ["cG_x04","cA_x04","cu_x02","cHW_x02","cWWMinuscB_x03"] 
#EFT_PARAMETERS = ["cG_x04","cHW_x02","cWWMinuscB_x03"] 
fitter = eft_fitter(EFT_PARAMETERS)

############### CHOOSE YOUR DATA SETS TO INCLUDE, no correlations between them ##############
fitter.processModel(model,"hzz")
#fitter.processModel(model_stxs1_hgg,"hgg")
#fitter.processModel(model_stxs0_hgg,"hgg")
#fitter.processModel(model_stxs_h4l,"hzz")
#############################################################################################

fitter.prep()
# Uncomment to set the other parameters in the model to their best fits in the fixed scans !
"""
fitter.global_fit()
fitter.reset()  # < - Now the nominal is at the best fits!
"""
# -------------------------------------------------------
#fitter.scan("cu_x02",0)
#fitter.scan("cG_x04")
#sys.exit()
#fitter.scan("cWWMinuscB_x03")
for e in EFT_PARAMETERS: fitter.scan(e)
fitter.scan2d("cG_x04","cHW_x02")
fitter.scan2d("cWWMinuscB_x03","cHW_x02")
fitter.scan2d("cWWMinuscB_x03","cG_x04")
fitter.scan2d("cu_x02","cG_x04")
fitter.scan2d("cA_x04","cG_x04")
fitter.scan2d("cA_x04","cHW_x02")
fitter.scan2d("cA_x04","cWWMinuscB_x03")
fitter.scan2d("cHW_x02","cu_x02")

