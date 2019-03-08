import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.nddata.utils import block_reduce
import os , sys, glob
from common_functions import colorplot
from scipy.special import erf
import subprocess
from scipy.optimize import curve_fit
from subprocess import check_output
import time
from scipy.interpolate import interp1d
from astropy.table import Table , Column ,vstack,hstack

print('a')

def sigmas_tab(name,n1,n2,pc,vo,N):

    psg_file='Files_model_2/'+name+'_psg.txt'
    print(psg_file)
    Parameters=np.loadtxt(psg_file,dtype=str)

    n1_min=float(Parameters[3])
    n1_max=float(Parameters[4])
    n2_min=float(Parameters[5])
    n2_max=float(Parameters[6])
    pc_min=float(Parameters[7])
    pc_max=float(Parameters[8])

    del Parameters

    aux= 1.0*(N-1)*1.0*(n1-n1_min)/(n1_max-n1_min)
    l1 = int(aux)
    u1 = l1+1

    number_file = l1/3;

    if l1%3==0: i1 = 0
    if l1%3==1: i1 = N*N
    if l1%3==2: i1 = 2*N*N

    aux = (N-1)*1.0*(n2-n2_min)/(n2_max-n2_min)

    l2 = int(aux)
    u2 = l2+1

    i2 = l2*N
    aux = (N-1)*1.0*(np.log(1.0*pc/pc_min)/np.log(1.0*pc_max/pc_min)) +1

    l3 = int(aux)
    i3 = l3

    fp = open('Files_model_2/'+name+'_sampler%03d.txt' %number_file)

    iline_1  = i1+i2+i3-1
    iline_2  = i1+i2+i3
    iline_3  = i1+i2+i3+N-1
    iline_4  = i1+i2+i3+N

    for i, line in enumerate(fp):
        if i == iline_1:
            sline_1 = line
        if i == iline_2:
            sline_2 = line
        if i == iline_3:
            sline_3 = line
        if i == iline_4:
            sline_4 = line
        if i>iline_4:
            break

    number_file = int((l1+1)/3);

    if (l1+1)%3==0: i1 = 0
    if (l1+1)%3==1: i1 = N*N
    if (l1+1)%3==2: i1 = 2*N*N

    aux = (N-1)*1.0*(n2-n2_min)/(n2_max-n2_min)

    l2 = int(aux)
    u2 = l2+1

    i2 = l2*N
    aux = (N-1)*1.0*(np.log(1.0*pc/pc_min)/np.log(1.0*pc_max/pc_min)) +1

    l3 = int(aux)
    i3 = l3

    fp = open('Files_model_2/'+name+'_sampler%03d.txt' %number_file)

    iline_5  = i1+i2+i3-1
    iline_6  = i1+i2+i3
    iline_7  = i1+i2+i3+N-1
    iline_8  = i1+i2+i3+N

    for i, line in enumerate(fp):
        if i == iline_5:
            sline_5 = line
        if i == iline_6:
            sline_6 = line
        if i == iline_7:
            sline_7 = line
        if i == iline_8:
            sline_8 = line
        if i>iline_8:
            break

    l1=sline_1.split('\t')
    l2=sline_2.split('\t')
    l3=sline_3.split('\t')
    l4=sline_4.split('\t')
    l5=sline_5.split('\t')
    l6=sline_6.split('\t')
    l7=sline_7.split('\t')
    l8=sline_8.split('\t')

    l1=np.array(l1,dtype=float)
    l2=np.array(l2,dtype=float)
    l3=np.array(l3,dtype=float)
    l4=np.array(l4,dtype=float)
    l5=np.array(l5,dtype=float)
    l6=np.array(l6,dtype=float)
    l7=np.array(l7,dtype=float)
    l8=np.array(l8,dtype=float)

    v1=1.0/np.abs((n1-l1[0])*(n2-l1[1])*(pc-l1[2]))
    v2=1.0/np.abs((n1-l2[0])*(n2-l2[1])*(pc-l2[2]))
    v3=1.0/np.abs((n1-l3[0])*(n2-l3[1])*(pc-l3[2]))
    v4=1.0/np.abs((n1-l4[0])*(n2-l4[1])*(pc-l4[2]))
    v5=1.0/np.abs((n1-l5[0])*(n2-l5[1])*(pc-l5[2]))
    v6=1.0/np.abs((n1-l6[0])*(n2-l6[1])*(pc-l6[2]))
    v7=1.0/np.abs((n1-l7[0])*(n2-l7[1])*(pc-l7[2]))
    v8=1.0/np.abs((n1-l8[0])*(n2-l8[1])*(pc-l8[2]))

    print(n1,l1[0],l2[0],l3[0],l4[0],l5[0],l6[0],l7[0],l8[0])
    print(n2,l1[1],l2[1],l3[1],l4[1],l5[1],l6[1],l7[1],l8[1])
    print(pc,l1[2],l2[2],l3[2],l4[2],l5[2],l6[2],l7[2],l8[2])
    vtot=v1+v2+v3+v4+v5+v6+v7+v8
    v1/=vtot
    v2/=vtot
    v3/=vtot
    v4/=vtot
    v5/=vtot
    v6/=vtot
    v7/=vtot
    v8/=vtot

    sigma=v1*l1[3:]+v2*l2[3:]+v3*l3[3:]+v4*l4[3:]+v5*l5[3:]+v6*l6[3:]+v7*l7[3:]+v8*l8[3:]
    return vo*sigma


eta=2*0.6744897501
kms_to_pcmyr=1.04589549342

cmdargs   = sys.argv
L         = float(cmdargs[-1])                        # Image size in parsecs
name_file = cmdargs[-2]

directory ='/data6/jutreras/Circulation-square/Bin_Tables/'+name_file[:-4]+'/'
vort_map  = np.load('../Maps/'+name_file+'_vort.npy')


N     = len(vort_map)                       # Size image
DX    = L/N                                      # Spatial resolution in pc
DA    = DX**2                                    # Area of pixel
pmin=4.0/L
pmax=(N-2.0)/(4.0*L)

del vort_map


lista = glob.glob(directory+'*'+name_file+'*')

rmax_array  = []
rmin_array  = []

for li in lista:
    elements=li.split('-')

    rmax_array.append(elements[-1])
    rmin_array.append(elements[-2])

rmax_set=sorted(set(rmax_array))

print(rmax_set)

for k, rmax in enumerate(rmax_set):
    f1         = plt.figure(figsize=(10,8))
    ax1        = f1.add_subplot(111)
    new_lista  = sorted(glob.glob(directory+'*'+name_file+'*-'+rmax))

    rmin       = new_lista[0].split('-')[-2]
    scales_str = []
    scales_flt = []
    resol      = ''
    for nl in new_lista:
        scales_str.append(nl.split('-')[-3])
        scales_flt.append(float(nl.split('-')[-3]))
        resol += nl.split('-')[-3] +' '
    walker_file = 'Files_model_2/'+name_file+'_'+rmin+'_'+rmax+'_mcmc'
    print(walker_file)
    print(os.path.isfile(walker_file))

    if os.path.isfile(walker_file):
        parametros = Table.read(walker_file,path='data')

        n_arr       = parametros['n1']
        m_arr       = parametros['n2']
        pc_arr      = parametros['pc']
        vo_arr      = parametros['vo']

        ##### exponent range #####

        N16, N50, N84 = np.percentile(n_arr, [16,50,84])

        ##### second exponent ####

        M16, M50, M84 = np.percentile(m_arr, [16,50,84])

        ##########################
        ##### kscale  range  #####

        K16, K50, K84 = np.percentile(pc_arr, [16,50,84])

        ##########################
        ##### velocity range #####

        V16, V50, V84 = np.percentile(vo_arr, [16,50,84])

        ##########################
        del parametros

        mean_walker = sigmas_tab(name_file,N50,M50,K50,V50,60)
        ##mean_walker*=DX**2

        wn_1  = sigmas_tab(name_file,N16,M50,K50,V50,60)
        wn_2  = sigmas_tab(name_file,N84,M50,K50,V50,60)
        wm_1  = sigmas_tab(name_file,N50,M16,K50,V50,60)
        wm_2  = sigmas_tab(name_file,N50,M84,K50,V50,60)
        wk_1  = sigmas_tab(name_file,N50,M50,K16,V50,60)
        wk_2  = sigmas_tab(name_file,N50,M50,K84,V50,60)
        wv_1  = sigmas_tab(name_file,N50,M50,K50,V16,60)
        wv_2  = sigmas_tab(name_file,N50,M50,K50,V84,60)

        std_walker=np.sqrt(((wv_2-wv_1)/2)**2 + ((wm_2-wm_1)/2)**2 + ((wn_2-wn_1)/2)**2 + ((wk_2-wk_1)/2)**2)

        #break
        n_scales        = len(new_lista)
        res             = np.zeros(n_scales,dtype=float)

        p10            = np.zeros(n_scales,dtype=float)
        p20            = np.zeros(n_scales,dtype=float)
        p30            = np.zeros(n_scales,dtype=float)
        p40            = np.zeros(n_scales,dtype=float)
        p50            = np.zeros(n_scales,dtype=float)
        p60            = np.zeros(n_scales,dtype=float)
        p70            = np.zeros(n_scales,dtype=float)
        p80            = np.zeros(n_scales,dtype=float)
        p90            = np.zeros(n_scales,dtype=float)

        q10            = np.zeros(n_scales,dtype=float)
        q20            = np.zeros(n_scales,dtype=float)
        q30            = np.zeros(n_scales,dtype=float)
        q40            = np.zeros(n_scales,dtype=float)
        q50            = np.zeros(n_scales,dtype=float)
        q60            = np.zeros(n_scales,dtype=float)
        q70            = np.zeros(n_scales,dtype=float)
        q80            = np.zeros(n_scales,dtype=float)
        q90            = np.zeros(n_scales,dtype=float)

        qu10            = np.zeros(n_scales,dtype=float)
        qu20            = np.zeros(n_scales,dtype=float)
        qu30            = np.zeros(n_scales,dtype=float)
        qu40            = np.zeros(n_scales,dtype=float)
        qu50            = np.zeros(n_scales,dtype=float)
        qu60            = np.zeros(n_scales,dtype=float)
        qu70            = np.zeros(n_scales,dtype=float)
        qu80            = np.zeros(n_scales,dtype=float)
        qu90            = np.zeros(n_scales,dtype=float)

        ql10            = np.zeros(n_scales,dtype=float)
        ql20            = np.zeros(n_scales,dtype=float)
        ql30            = np.zeros(n_scales,dtype=float)
        ql40            = np.zeros(n_scales,dtype=float)
        ql50            = np.zeros(n_scales,dtype=float)
        ql60            = np.zeros(n_scales,dtype=float)
        ql70            = np.zeros(n_scales,dtype=float)
        ql80            = np.zeros(n_scales,dtype=float)
        ql90            = np.zeros(n_scales,dtype=float)

        radius  = 0.5*(float(new_lista[0].split('-')[-1])+float(new_lista[0].split('-')[-2]))

        for ij, nl in enumerate(new_lista):
            res[ij] = float(nl.split('-')[-3])
            R       = res[ij]/res[0]
            tabla   = Table.read(nl,path='data')
            X       = tabla['vort']/R**2/DA
            Y       = tabla['av_vort']/R**2/DA

            Z       = Y + np.random.normal(0,mean_walker[ij],size=len(Y))/DA
            Zu      = Y + np.random.normal(0,mean_walker[ij]+std_walker[ij],size=len(Y))/DA

            if mean_walker[ij]-std_walker[ij] <=0:
                Zl  = Y
            else:
                Zl  = Y + np.random.normal(0,mean_walker[ij]-std_walker[ij],size=len(Y))/DA

            p10[ij], p20[ij], p30[ij], p40[ij], p50[ij], p60[ij], p70[ij], p80[ij], p90[ij] = np.percentile(X,[10,20,30,40,50,60,70,80,90])
            q10[ij], q20[ij], q30[ij], q40[ij], q50[ij], q60[ij], q70[ij], q80[ij], q90[ij] = np.percentile(Z,[10,20,30,40,50,60,70,80,90])

            qu10[ij], qu20[ij], qu30[ij], qu40[ij], qu50[ij], qu60[ij], qu70[ij], qu80[ij], qu90[ij] = np.percentile(Zu,[10,20,30,40,50,60,70,80,90])
            ql10[ij], ql20[ij], ql30[ij], ql40[ij], ql50[ij], ql60[ij], ql70[ij], ql80[ij], ql90[ij] = np.percentile(Zl,[10,20,30,40,50,60,70,80,90])

        psum=p10+p20+p30+p40+p50+p60+p70+p80+p90
        qsum=q10+q20+q30+q40+q50+q60+q70+q80+q90

        delta=(psum-qsum)/9.0

        ax1.plot(res, p10, color=colorplot(11,1), label='simulation')
        ax1.plot(res, p20, color=colorplot(11,2), label='')
        ax1.plot(res, p30, color=colorplot(11,3), label='')
        ax1.plot(res, p40, color=colorplot(11,4), label='')
        ax1.plot(res, p50, color=colorplot(11,5), label='')
        ax1.plot(res, p60, color=colorplot(11,6), label='')
        ax1.plot(res, p70, color=colorplot(11,7), label='')
        ax1.plot(res, p80, color=colorplot(11,8), label='')
        ax1.plot(res, p90, color=colorplot(11,9), label='')

        ax1.plot(res, q10+delta, color=colorplot(11,1), label='model',linestyle='-.')
        ax1.fill_between(res, ql10+delta, qu10+delta, color=colorplot(11,1), alpha=0.2)

        ax1.plot(res, q20+delta, color=colorplot(11,2), label='',linestyle='-.')
        ax1.fill_between(res, ql20+delta, qu20+delta, color=colorplot(11,2), alpha=0.2)

        ax1.plot(res, q30+delta, color=colorplot(11,3), label='',linestyle='-.')
        ax1.fill_between(res, ql30+delta, qu30+delta, color=colorplot(11,3), alpha=0.2)

        ax1.plot(res, q40+delta, color=colorplot(11,4), label='',linestyle='-.')
        ax1.fill_between(res, ql40+delta, qu40+delta, color=colorplot(11,4), alpha=0.2)

        ax1.plot(res, q50+delta, color=colorplot(11,5), label='',linestyle='-.')
        ax1.fill_between(res, ql50+delta, qu50+delta, color=colorplot(11,5), alpha=0.2)

        ax1.plot(res, q60+delta, color=colorplot(11,6), label='',linestyle='-.')
        ax1.fill_between(res, ql60+delta, qu60+delta, color=colorplot(11,6), alpha=0.2)

        ax1.plot(res, q70+delta, color=colorplot(11,7), label='',linestyle='-.')
        ax1.fill_between(res, ql70+delta, qu70+delta, color=colorplot(11,7), alpha=0.2)

        ax1.plot(res, q80+delta, color=colorplot(11,8), label='',linestyle='-.')
        ax1.fill_between(res, ql80+delta, qu80+delta, color=colorplot(11,8), alpha=0.2)

        ax1.plot(res, q90+delta, color=colorplot(11,9), label='',linestyle='-.')
        ax1.fill_between(res, ql90+delta, qu90+delta, color=colorplot(11,9), alpha=0.2)

        ax1.set_ylabel(r'$\Gamma/ \Delta^2 \ {\rm 1/Myr}$',size=20)
        ax1.set_xlabel(r'Scale [pc]',size=20)
        ax1.set_xscale('log')
        ax1.legend(loc='upper right', prop={'size':15}, ncol = 2)
        f1.savefig('Plots/'+name_file+'_%05d_percentiles.png' %int(radius))
        plt.close(f1)

        os.system('rclone copy Plots/'+name_file+'_%05d_percentiles.png uchile:Double_Power/Percentiles/'  %int(radius))

        tabla=Table()
        tabla['resolution']   = Column(np.array(res)       , description='resolution in pc')

        tabla['p10']          = Column(np.array(p10))
        tabla['p20']          = Column(np.array(p20))
        tabla['p30']          = Column(np.array(p30))
        tabla['p40']          = Column(np.array(p40))
        tabla['p50']          = Column(np.array(p50))
        tabla['p60']          = Column(np.array(p60))
        tabla['p70']          = Column(np.array(p70))
        tabla['p80']          = Column(np.array(p80))
        tabla['p90']          = Column(np.array(p90))

        tabla['q10']          = Column(np.array(q10+delta))
        tabla['q20']          = Column(np.array(q20+delta))
        tabla['q30']          = Column(np.array(q30+delta))
        tabla['q40']          = Column(np.array(q40+delta))
        tabla['q50']          = Column(np.array(q50+delta))
        tabla['q60']          = Column(np.array(q60+delta))
        tabla['q70']          = Column(np.array(q70+delta))
        tabla['q80']          = Column(np.array(q80+delta))
        tabla['q90']          = Column(np.array(q90+delta))

        tabla['ql10']          = Column(np.array(ql10+delta))
        tabla['ql20']          = Column(np.array(ql20+delta))
        tabla['ql30']          = Column(np.array(ql30+delta))
        tabla['ql40']          = Column(np.array(ql40+delta))
        tabla['ql50']          = Column(np.array(ql50+delta))
        tabla['ql60']          = Column(np.array(ql60+delta))
        tabla['ql70']          = Column(np.array(ql70+delta))
        tabla['ql80']          = Column(np.array(ql80+delta))
        tabla['ql90']          = Column(np.array(ql90+delta))

        tabla['qu10']          = Column(np.array(qu10+delta))
        tabla['qu20']          = Column(np.array(qu20+delta))
        tabla['qu30']          = Column(np.array(qu30+delta))
        tabla['qu40']          = Column(np.array(qu40+delta))
        tabla['qu50']          = Column(np.array(qu50+delta))
        tabla['qu60']          = Column(np.array(qu60+delta))
        tabla['qu70']          = Column(np.array(qu70+delta))
        tabla['qu80']          = Column(np.array(qu80+delta))
        tabla['qu90']          = Column(np.array(qu90+delta))

        tabla.write('Tables/'+name_file+'_%05d_percentiles'%int(radius) ,path='data',format='hdf5',overwrite=True)
        os.system('rclone copy Tables/'+name_file+'_%05d_percentiles uchile:Double_Power/Percentiles/'  %int(radius))
