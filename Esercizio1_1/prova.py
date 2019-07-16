import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math


y,erry,x = np.loadtxt("risultati1.dat",unpack=True)
y1,erry1,x1 = np.loadtxt("risultati2.dat",unpack=True)


ax1 = plt.subplot(211)
ax1.plot(x,y,linewidth= 1, label='risultati1', lw=3)
ax1.errorbar(x,y,erry)
plt.xlabel('$block$')
plt.ylabel('$<r>$')

ax2 = plt.subplot(212)
ax2.plot(x1,y1,linewidth= 1, label='risultati2',lw=20)
ax2.errorbar(x1,y1,erry1)
plt.xlabel('$block$')
plt.ylabel('$<(r-1/2)^2>$')


plt.show()
