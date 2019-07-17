import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

y,erry,x = np.loadtxt("ave_temp.out",unpack=True)
#y1,erry1,x1 = np.loadtxt("S_directput.dat",unpack=True)
#y2,erry2,x2 = np.loadtxt("S_discretizedcall.dat",unpack=True)
#y3,erry3,x3 = np.loadtxt("S_discretizedput.dat",unpack=True)

spettro=plt.figure()
ax=plt.axes()
ax.plot(x,y,linewidth= 1, label='risultati', lw=3)
ax.errorbar(x,y,erry)
plt.xlabel('$block$')
plt.ylabel('$<\pi>$')
plt.xlim([0,104])
plt.ylim([-0.7,2.5])
plt.title('Temperatura ')


#ax2 = plt.subplot(222)
#ax2.plot(x1,y1,linewidth= 1, label='risultati2',lw=30)
#ax2.errorbar(x1,y1,erry1)
#plt.xlabel('$block$')
#plt.ylabel('$<P>$')
#plt.xlim([0,110])
#plt.ylim([5.0,6.0])
#plt.title('European call-option prices, $P[S(0),0]$ (direct) ')

#ax3 = plt.subplot(223)
#ax3.plot(x2,y2,linewidth= 1, label='risultati1', lw=2)
#ax3.errorbar(x2,y2,erry2)
#plt.xlabel('$block$')
#plt.ylabel('$<C>$')
#plt.xlim([0,110])
#plt.ylim([14.4,15.75])
#plt.title('European call-option prices, $C[S(0),0]$ (discretized) ')


#ax4 = plt.subplot(224)
#ax4.plot(x3,y3,linewidth= 1, label='risultati2',lw=30)
#ax4.errorbar(x3,y3,erry3)
#plt.xlabel('$block$')
#plt.ylabel('$<P>$')
#plt.xlim([0,110])
#plt.ylim([5.0,6.0])
#plt.title('European call-option prices, $P[S(0),0]$ (discretized) ')

#plt.subplots_adjust(hspace = 0.4)




plt.show()
