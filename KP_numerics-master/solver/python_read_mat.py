import scipy.io as sp
import numpy as np

param = sp.loadmat('Example/parameters.mat');
print(param['Ny'])

for n in range(0,np.size(param['t'])):
    d = sp.loadmat('Example/'+'{:05d}'.format(n)+'.mat');
    