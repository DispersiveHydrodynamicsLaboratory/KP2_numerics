import sys
import os
import glob
import pdb
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt


def main():
    main_dir = '/home/magbags/Downloads/Research/Numerics/KP2/';
    data_dir = main_dir + 'KP3/';
    print('Directory '+main_dir+'exists?');
    fnames = sorted(glob.glob(data_dir+"0*.mat"));
    # remove IC
    fnames = fnames[1:len(fnames)];
    # print(fnames);
    # pdb.set_trace()
    
    # Load numerical parameters
    param = sio.loadmat(data_dir + 'parameters.mat');
    # print(param['Ny']); # Print parameter for debugging
    
    # Reassign parameters needed
    Lx = np.asscalar(param['Lx']);
    Nx = np.asscalar(param['Nx']);
    Ly = np.asscalar(param['Ly']);
    Ny = np.asscalar(param['Ny']);
    
    # Initialize domain
    # ISSUE: negatives don't work the way you think!
    x = (2*Lx/Nx)*np.linspace(Nx*-1/2,Nx*1/2-1,Nx);
    y = (2*Ly/Ny)*np.linspace(Ny*-1/2,Ny*1/2-1,Ny);
    [X,Y] = np.meshgrid(x,y);
    
    for n in fnames: #range(0,np.size(param['t'])):
        d = sio.loadmat(n);
        print(d['tnow']);
        # pdb.set_trace()
        # Plot
        plt.figure(1)
        plt.contourf(X,Y,d['u'])
        plt.colorbar()
        plt.show()
        # plt.savefig('out.png')
        
if __name__=='__main__':
    sys.exit(main())
    