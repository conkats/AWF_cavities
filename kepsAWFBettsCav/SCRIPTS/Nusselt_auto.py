# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:09:06 2020

@author: Constantinos katsamis
requirement: python3, pandas, numpy lib
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#allows to pass arguments when calling the script
import argparse
import sys

#i.e python3 scriptName.py firstArgument secondArgument 
#Provide the case and experimental data path to plot Nu profile
runDirectory = sys.argv[1]
expDirecotry = sys.argv[2] 
#eg python3 Nusselt_auto.py ../RESU/BETTS_KEAWF_i0_/ hot_wall_heatflux_engauge_digitizer_app_loRa.csv

pathsim =['../RESU/'+runDirectory+'/']

fig = plt.figure(figsize = (10.0,8.0))
fig.subplots_adjust(top=0.88, hspace = 0.35)
axes1 = fig.add_subplot(1, 1, 1)

modelList = ['$k-{\epsilon}$ AWF SGDH']           
markerColorList=['r*-']

for IndexPath in range(len(pathsim)):

    qexpdataRead = pd.read_csv(expDirecotry,sep=',')
    q_betts = qexpdataRead.iloc[0:,0].values
    y_nondim_betts = qexpdataRead.iloc[0:,1].values
    nu_betts = (q_betts*0.076)/(19.6*0.0253)
    axes1.plot(y_nondim_betts, nu_betts, 'ko', label='Exp. Betts and Bokhari',ms=9,markeredgewidth=3)

    for modelIndex in range(len(modelList)):
        simRead = pd.read_table(pathsim[IndexPath]+'Nusselt.dat',comment='#', sep='\s+')
        y = simRead.iloc[0:,0].values
        y_non_dim = y/2.18 
        nusselt=(simRead.iloc[0:,1].values)
        nu_mean=np.mean([nusselt])
        
        axes1.plot(y_non_dim, nusselt, markerColorList[modelIndex], label=modelList[modelIndex],linewidth=3,ms=5)

axes1.xaxis.set_tick_params(labelsize=24)
axes1.yaxis.set_tick_params(labelsize=24)
axes1.legend(loc='upper right', shadow=True, fontsize = 20)
axes1.set_xlabel('$y/h$',fontsize =20)
axes1.set_ylabel('$Nu_{hot}$',fontsize =20)
fig.suptitle(r'$Rayleigh~0.86x10^6$',  fontsize =22)
plt.savefig('../RESU/'+ runDirectory +'/Nusselt.png',dpi=200,bbox_inches="tight")
