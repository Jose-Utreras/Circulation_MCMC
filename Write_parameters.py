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

R1=cmdargs[-2]
R2=cmdargs[-1]

Rmax=float(cmdargs[-1])
Rmin=float(cmdargs[-2])
L=float(cmdargs[-3])                        # Image size in parsecs
name_file=cmdargs[-4]
snapshot=name_file
print(1)
directory='/data6/jutreras/Circulation-square/Bin_Tables/'+name_file[:-4]+'/'

rbin_string='-%05d' %int(Rmin)+ '-%05d' %int(Rmax)
rbin_prefix='_%05d' %int(Rmin)+ '_%05d' %int(Rmax)

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
dv_min=-5.0
dv_max=5.0

sigma_t=1.0*kms_to_pcmyr*2*DX               # total error in circulation

"""
Creating the resolution array which are integers factors
"""
print(2)
Scales=np.exp(np.linspace(np.log(2),np.log(0.125*N),50))
Scales=np.insert(Scales,0,1)

lista=glob.glob(directory+'*'+name_file+'*'+rbin_string)
lista=sorted(lista)

#print('lista')
#print(lista)

resolutions=np.array(Scales*DX)

Number_of_resolutions=len(Scales)


stds=np.zeros_like(resolutions,dtype=float)
serr=np.zeros_like(resolutions,dtype=float)

ratio=[]

res_temp=np.exp(np.linspace(np.log(0.75*Scales.min()),np.log(1.5*Scales.max()),1000))

aux_res=np.insert(Scales,0,Scales[0]/2.0)
aux_res=np.insert(aux_res,len(aux_res),Scales[-1]*2)

name_file+=rbin_prefix
print(3)
for ij in range(Number_of_resolutions):
    R		  = Scales[ij]
    c_sig=sigma_t/R**1.5
    tabla	  = Table.read(lista[ij],path='data')
    vorticity = tabla['vort']/R**2
    omega     = tabla['av_vort']/R**2

    #print(ij, "cuociente de las sumas", vorticity.sum()/omega.sum())
    #print(ij, "suma del cuociente", (vorticity/omega).sum())
    mapa=vorticity-omega
    mapa-=np.mean(mapa)
    #print(lista[ij],np.mean(mapa))
    try:
        c_len=len(mapa)
        pu,pl=np.percentile(mapa,[75,25])
        c_width=16*(pu-pl)*c_len**(-0.333)
        c_bins=1.0/c_width

        c_max=mapa.max()+0.5*c_width
        c_min=mapa.min()-0.5*c_width
        c_bins=max(int((c_max-c_min)*c_bins),4)
        c_width=1.0*(c_max-c_min)/(1.0*c_bins)

    except:
        c_len=len(mapa)
        c_max=mapa[0]+0.5*c_width
        c_min=mapa[0]-0.5*c_width
        c_bins=4
        c_width=1.0*(c_max-c_min)/(1.0*c_bins)
    if c_len<2:
        c_max=mapa[0]+2*c_sig
        c_min=mapa[0]-2*c_sig
        c_bins=4
        c_width=1.0*(c_max-c_min)/(1.0*c_bins)

    np.savetxt('Files_model_2/'+name_file+'_circ_vector.txt',mapa)
    parameters = name_file +' '+str(c_min)+' '+str(c_max)+' '+str(c_bins)+' '
    parameters+= str(c_width)+' '+str(c_sig)+' '+str(c_len)
    subprocess.call(['./histogram_m_2 '+parameters], stdin = sys.stdin, shell=True)
            
    y,errors=np.loadtxt('Files_model_2/'+name_file+'_hist_data.txt',unpack=True)
    x=np.linspace(c_min,c_max,c_bins+1)
    x=0.5*(x[1:]+x[:-1])
    errors=np.sqrt(y+errors**2)
    errors[errors<1.0]=1.0

    std_lim=np.average((x-np.average(x,weights=y))**2,weights=y)**0.5
    a_min,a_max=np.average(y),y.max()*2
    xaux=np.insert(x,0,2*x[0]-x[1])
    yaux=np.insert(y,0,0)
    ysum=np.cumsum(yaux)/yaux.sum()
    fsum=interp1d(ysum,xaux)
    mean_lim=max(np.abs(fsum(0.6)),np.abs(fsum(0.4)))
    mu_min,mu_max=-mean_lim,mean_lim
    sig_lim=0.5*(fsum(0.84)-fsum(0.16))
    sig_min,sig_max=min(sig_lim,std_lim)/2.0,1.5*max(sig_lim,std_lim)

   
    try:
        popt, pcov   = curve_fit(gaussian,x,y,sigma=errors,bounds=([a_min,mu_min,sig_min],[a_max,mu_max,sig_max]))
        popt1,pcov1  = curve_fit(gaussian,x,y,sigma=errors,method='trf',bounds=([a_min,mu_min,sig_min],[a_max,mu_max,sig_max]))
        chi_0        = (y-gaussian(x,*popt))/errors
        chi_1        = (y-gaussian(x,*popt1))/errors
        chi_0        = chi_0**2
        chi_1        = chi_1**2
        chi_0=chi_0.sum()
        chi_1=chi_1.sum()
        if chi_0>chi_1:
            popt=popt1
            pcov=pcov1
    except:
        pass
    try:
        stds[ij]=np.abs(popt[-1])
        serr[ij]=np.sqrt(pcov.diagonal()[-1])
    except:
        stds[ij]=c_sig
        serr[ij]=c_sig


    outlier=False
    if ij>2:
        aux_ratio = np.array(ratio)
        mu        = np.mean(aux_ratio)
        dev       = np.std(aux_ratio)
        QQ        = 0.5*(c_max-c_min)/stds[ij]
        if (QQ<mu-5*dev)|(mu+5*dev<QQ):
            outlier=True
        del aux_ratio,QQ

    ratio.append(0.5*(c_max-c_min)/stds[ij])
    if outlier:
        print('HAY UN OUTLIER!!!!')
        stds[ij]=0.5*(c_max-c_min)/mu
        serr[ij]=dev

#print(Scales)
print(4)
#print(serr)
if (serr[-1]>100.0)|(np.isnan(serr[-1])):
    serr[-1]=stds[-1]
infinitos=(serr>100.0)|(serr==np.inf)|(np.isnan(serr))
finitos=~infinitos
function=interp1d(Scales[finitos],serr[finitos])
serr[infinitos]=function(Scales[infinitos])

aux_stds=np.insert(stds,0,stds[0])
aux_stds=np.insert(aux_stds,len(aux_stds),stds[-1])

fun=interp1d(aux_res,aux_stds)
stds_temp=fun(res_temp)

smooth_stds=Smooth(np.log(res_temp),stds_temp,n=40,b=1)

fun_smooth=interp1d(res_temp,smooth_stds)

diff=np.abs(stds_temp-smooth_stds)
dev=np.std(diff[diff>0])
diff=np.abs(stds-fun_smooth(Scales))
diff=np.sqrt(dev**2+diff**2)
stds=fun_smooth(Scales)
stot=np.sqrt(serr**2+diff**2)

#print(stds)
#print(stot)

#stot=np.sqrt(stot**2+(0.1*stds)**2)
stot=np.maximum(stot,0.05*stds)

#print(stot)

frac_stot=stot/stds
right_region=np.where(Scales>Scales[stds==stds.max()])
#print('right region', right_region)
best=np.where(frac_stot==frac_stot[right_region].min())[0][0]

stds_string='stds'
stot_string='stot'
print(5)
for ele1,ele2 in zip(stds,stot):
    stds_string+=','+str(ele1)
    stot_string+=','+str(ele2)
stds_string+='\n'
stot_string+='\n'

#print(stds_string)
#print(stot_string)

Resolution=''
sigma_array=''
error_array=''
for sa,re,se in zip(stds,resolutions,stot):
    sigma_array+=' '+str(sa)
    Resolution+=' '+str(re)
    error_array+=' '+str(se)

print('Writing parameters')
parameters=' '+str(n1_min)+' '+str(n1_max)+' '+str(n2_min)+' '+str(n2_max)+' '+str(p_min*2)+' '+str(p_max*0.5)+' '+str(dv_min)+' '+str(dv_max)
linea=name_file+' '+snapshot+' '+str(L)+' '+str(N) +' '+ parameters + Resolution + sigma_array + error_array
variables=open('Files_model_2/'+name_file+'_reduce_parameters_m_2.txt','w')
variables.write(linea)
variables.close()

