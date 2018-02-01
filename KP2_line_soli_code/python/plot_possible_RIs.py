# Plot possible Riemann Invariant combinations
# Given m

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# m, angle over which jump occurs, chosen
m = -1

# Initialize R1 and R2
r1 = np.linspace(-5,5,100)
r2 = r1
(R1r,R2l) = np.meshgrid(r1,r2)

# Middle Speed options
ms1 = int( (R1r+R2l>=0) & (R2l>=-1/m) & (R1r<=2/m+R2l) )
ms2 = int( (R1r+R2l<=0) & (R2l<=-1/m) & (R1r>=2/m+R2l) )
ms  = (ms1*ms2)

# Plot middle speed options
plt.figure(1)
plt.plot(R1r,R2l,ms)
plt.show()
plt.savefig('out.png')
   
