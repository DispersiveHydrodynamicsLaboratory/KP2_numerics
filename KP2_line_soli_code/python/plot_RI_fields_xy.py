# Plot possible Riemann Invariant combinations
# Given m

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Initialize Da = (a_- - a_+) and Dq = q_l - q_r
Da      = np.linspace(-5,5,101)
Dq      = Da
[DA,DQ] = np.meshgrid(Da,Dq)
RWposs  = 1*(np.abs(DA)<=DQ)+1*(np.abs(DA)<=-DQ)
RWposs[50,50] = 1


# Plot middle speed options
plt.figure(1)
plt.contourf(DA,DQ,RWposs)
plt.colorbar()
plt.show()
plt.savefig('out.png')
