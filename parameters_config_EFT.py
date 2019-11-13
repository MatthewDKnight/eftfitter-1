
PARAMS = {
    "clW_x02"           :[[-5 , 5    ]		,0,0]    # current value and nominal / resettable value 
    ,"cHu_x02"          :[[-1.1 , 1.1]		,0,0]
    ,"c2W"              :[[-5 , 5  	]	,0,0]
    ,"cl"               :[[-5 , 5  	]	,0,0]
    ,"cdG_x02"          :[[-5 , 5  	]	,0,0]
    ,"cH_x01"           :[[-1.4 , 1.94]		,0,0]
    ,"cpHL_x02"         :[[-5 , 5]		,0,0]
    ,"c2B"              :[[-5 , 5 ] 		,0,0]
    ,"cG_x05"     	:[[-10. ,10.]		,0,0]
    ,"cA_x04"     	:[[-10 , 10 ]		,0,0]
    ,"cWWMinuscB_x02"   :[[-15 , 15  	]	,0,0]
    ,"cu_x01"           :[[-20. , 10. ]		,0,0]
    ,"cuW_x02"          :[[-20. , 30.  	]	,0,0]
    ,"tcA_x04"          :[[-12 , 12 	]	,0,0]
    ,"cT_x03"           :[[-4.3 , 3.3 ] 	,0,0]
    ,"tc3W_x01"         :[[-1.8 , 1.8  ]	,0,0]
    ,"cWWPluscB_x03"    :[[-3.3 , 1.8  ]	,0,0]
    ,"cpHQ_x03"         :[[-4.4 , 4.4  ]	,0,0]
    ,"cHud_x02"         :[[-5 , 5  	]	,0,0]
    ,"cHe_x03"          :[[-1.8 , 0.25] 	,0,0]
    ,"tcHB_x01"         :[[-2.4 , 2.4]		,0,0]
    ,"cHQ_x03"          :[[-1.9 , 6.9 ] 	,0,0]
    ,"c3W_x02"          :[[-8.3 , 4.5  ]	,0,0]
    ,"cuB_x02"          :[[-5 , 5 	]	,0,0]
    ,"c2G_x04"          :[[-1.6 , 1.6 ] 	,0,0]
    ,"cHB_x02"          :[[-4.5 , 7.5  ]	,0,0]
    ,"c3G_x04"          :[[-1.6 , 1.6  ]	,0,0]
    ,"cdW_x02"          :[[-5 , 5  	]	,0,0]
    ,"cHW_x02"          :[[-12. , 16.]		,0,0]
    ,"c6"               :[[-5 , 5  	]	,0,0]
    ,"tcHW_x02"         :[[-6 , 6  	]	,0,0]
    ,"tcG_x04"          :[[-1.2 , 1.2]		,0,0]
    ,"cHL_x02"          :[[-5 , 5  	]	,0,0]
    ,"cdB_x02"          :[[-5 , 5  	]	,0,0]
    ,"cHd_x02"          :[[-4.2 , 0.44] 	,0,0]
    ,"cd_x01"           :[[-20. , 10. ]		,0,0]
    ,"clB_x02"          :[[-5 , 5  	]	,0,0]
    ,"cuG_x02"          :[[-5 , 5	]	,0,0]
    ,"tc3G_x04"         :[[-1.6 , 1.6]		,0,0]
    }

"""
SCALING_FUNC_STR = "stxstoeft_scaling"
COMBINE_WS = "summer2019/result.root"  # <- output from text2workspace with scaling functions 
"""
COMBINE_WS="result.root"
SCALING_FUNC_STR ="scaling"
#MYPARAMS = ["cG_x05","cA_x04","cu_x01","cHW_x02","cWWMinuscB_x02","cd_x01","cl"]  
MYPARAMS = ["cG_x05","tcG_x04","c2G_x04","c3G_x04"] 
#MYPARAMS = ["cG_x05","tcG_x04"]
