import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.nddata.utils import block_reduce
import os , sys, glob
from common_functions import simple_symmetric_map, correct_map, clean_data, gaussian, Smooth, shift_image
from scipy.special import erf
import subprocess
from scipy.optimize import curve_fit
from subprocess import check_output
import time
from scipy.interpolate import interp1d
from astropy.table import Table , Column ,vstack,hstack
 
eta=2*0.6744897501
kms_to_pcmyr=1.022

cmdargs = sys.argv


L=float(cmdargs[-1])                        # Image size in parsecs
name_file=cmdargs[-2]

vort_map=np.load('../Maps/'+name_file+'_vort.npy')

N=len(vort_map)                       # Size image
DX=L/N                                      # Spatial resolution in pc
DA=DX**2                                    # Area of pixel
p_min=4.0/L
p_max=1.0/(4.0*DX)
n1_min=0.801#1.333-0.8
n1_max=1.901#1.333+0.8
n2_min=2.0153#2.000-0.8
n2_max=2.9801#2.000+0.8

sigma_t=1.0*kms_to_pcmyr*2*DX               # total error in circulation

"""
Creating the resolution array which are integers factors
"""

Scales=np.exp(np.linspace(np.log(2),np.log(0.125*N),50))
Scales=np.insert(Scales,0,1)


resolutions=np.array(Scales*DX)


Resolution=''
for re in resolutions:
    Resolution+=' '+str(re)

print('Writing parameters')
parameters=' '+str(n1_min)+' '+str(n1_max)+' '+str(n2_min)+' '+str(n2_max)+' '+str(p_min*2)+' '+str(p_max*0.5)
linea=name_file+' '+str(L)+' '+str(N) +' '+ parameters + Resolution
variables=open('Files_model_2/'+name_file+'_psg.txt','w')
variables.write(linea)
variables.close()

