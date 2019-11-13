# Stage-0 Combination from CMS 2016 (36/fb) dataset - HIG-17-031 

# Best-fit values and (symmetrized) uncertainties from : http://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-17-031/CMS-HIG-17-031_Table_005.pdf

# Interpret each XS as decay to ZZ --> will be scaled by BR(H->ZZ) as in measurement from CMS 

X = {
    'GG2H_PTH_0_20'                  : [[], 1,1.25,0.28]
   ,'GG2H_PTH_20_45'                 : [[], 1,0.67,0.37]
   ,'GG2H_PTH_45_80'                 : [[], 1,1.35,0.35]
   ,'GG2H_PTH_80_120'                : [[], 1,0.80,0.46]
   ,'GG2H_PTH_120_200'               : [[], 1,1.15,0.46]
   ,'GG2H_PTH_GT200'               : [[], 1,0.70,0.64]
}

# Correlation values from: http://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-17-031/CMS-HIG-17-031_Figure-aux_002.pdf
correlation = {
('GG2H_PTH_0_20','GG2H_PTH_0_20')    : 1.        
,('GG2H_PTH_0_20','GG2H_PTH_20_45')  : -0.598     
,('GG2H_PTH_0_20','GG2H_PTH_45_80')  : 0.246
,('GG2H_PTH_0_20','GG2H_PTH_80_120') : 0.022 
,('GG2H_PTH_0_20','GG2H_PTH_120_200'): 0.011 
,('GG2H_PTH_0_20','GG2H_PTH_GT200'): 0.014
      
,('GG2H_PTH_20_45','GG2H_PTH_20_45')  : 1.    
,('GG2H_PTH_20_45','GG2H_PTH_45_80')  : -0.368
,('GG2H_PTH_20_45','GG2H_PTH_80_120') : 0.104
,('GG2H_PTH_20_45','GG2H_PTH_120_200'): 0.028 
,('GG2H_PTH_20_45','GG2H_PTH_GT200'): 0.003
   
,('GG2H_PTH_45_80','GG2H_PTH_45_80')  : 1.
,('GG2H_PTH_45_80','GG2H_PTH_80_120') : -0.211
,('GG2H_PTH_45_80','GG2H_PTH_120_200'): 0.082
,('GG2H_PTH_45_80','GG2H_PTH_GT200'): 0.034

,('GG2H_PTH_80_120','GG2H_PTH_80_120') : 1.
,('GG2H_PTH_80_120','GG2H_PTH_120_200'): -0.117
,('GG2H_PTH_80_120','GG2H_PTH_GT200'): 0.045

,('GG2H_PTH_120_200','GG2H_PTH_120_200'): 1
,('GG2H_PTH_120_200','GG2H_PTH_GT200'): -0.025

,('GG2H_PTH_GT200','GG2H_PTH_GT200'): 1
}
