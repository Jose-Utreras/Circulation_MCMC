import numpy as np
import os , sys, glob
from astropy.table import Table , Column ,vstack,hstack
import h5py

cmdargs = sys.argv
name_file=cmdargs[-1]


os.system("cat Files_model_2/"+name_file+"_third[0-9]* > Files_model_2/"+name_file+"_mcmc.txt")

n1,n2,pc,vo,prob = np.loadtxt("Files_model_2/"+name_file+"_mcmc.txt",unpack=True)

tabla=Table()
tabla['n1'] = Column(np.array(n1))
tabla['n2'] = Column(np.array(n2))
tabla['pc'] = Column(np.array(pc))
tabla['vo'] = Column(np.array(vo))
tabla['prob'] = Column(np.array(prob))
tabla.write('Files_model_2/'+name_file+'_mcmc',path='data',format='hdf5',overwrite=True)


