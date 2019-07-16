import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.pyplot import *
from matplotlib import rc

#bs1,vs1,be1,ve1,be1,ve1= np.loadtxt("isto_1.dat", delimiter=' ',unpack=True)
#bs2,vs2,be2,ve2,be2,ve2= np.loadtxt("isto_2.dat", delimiter=' ',unpack=True)
#bs10,vs10,be10,ve10,be10,ve10 = np.loadtxt("isto_10.dat", delimiter=' ',unpack=True)
#bs100,vs100,be100,ve100,be100,ve100 = np.loadtxt("isto_100.dat", delimiter=' ',unpack=True)

plt.rc('text', usetex=False)
plt.rc('font', family='arial')

br1,vr1= np.loadtxt("isto_1.dat", usecols=(0,1), delimiter=' ', unpack='true')
br2,vr2= np.loadtxt("isto_2.dat", usecols=(0,1), delimiter=' ', unpack='true')
br10,vr10= np.loadtxt("isto_10.dat", usecols=(0,1), delimiter=' ', unpack='true')
br100,vr100= np.loadtxt("isto_100.dat", usecols=(0,1), delimiter=' ', unpack='true')

be1,ve1= np.loadtxt("isto_1.dat", usecols=(2,3), delimiter=' ', unpack='true')
be2,ve2= np.loadtxt("isto_2.dat", usecols=(2,3), delimiter=' ', unpack='true')
be10,ve10= np.loadtxt("isto_10.dat", usecols=(2,3), delimiter=' ', unpack='true')
be100,ve100= np.loadtxt("isto_100.dat", usecols=(2,3), delimiter=' ', unpack='true')

bl1,vl1= np.loadtxt("isto_1.dat", usecols=(4,5), delimiter=' ', unpack='true')
bl2,vl2= np.loadtxt("isto_2.dat", usecols=(4,5), delimiter=' ', unpack='true')
bl10,vl10= np.loadtxt("isto_10.dat", usecols=(4,5), delimiter=' ', unpack='true')
bl100,vl100= np.loadtxt("isto_100.dat", usecols=(4,5), delimiter=' ', unpack='true')



plt.figure(1)

# linear
plt.subplot(221)
plt.bar(br1,vr1)
plt.title('random n=1')
plt.grid(True)


# log
plt.subplot(222)
plt.bar(br2,vr2)
plt.title('random n=2')
plt.grid(True)


# symmetric log
plt.subplot(223)
plt.bar(br10,vr10)
plt.title('random n=10')
plt.grid(True)

# logit
plt.subplot(224)
plt.bar(br100,vr100)
plt.title('random n=100')
plt.grid(True)
# Format the minor tick labels of the y-axis into empty strings with
# `NullFormatter`, to avoid cumbering the axis with too many labels.
plt.gca().yaxis.set_minor_formatter(NullFormatter())
# Adjust the subplot layout, because the logit one may take more space
# than usual, due to y-tick labels like "1 - 10^{-3}"
plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)

###################################################esponenziale#######################################

plt.figure(2)

# linear
plt.subplot(221)
plt.bar(be1,ve1)
plt.title('Esponenziale n=1')
plt.grid(True)


# log
plt.subplot(222)
plt.bar(be2,ve2)
plt.title('Esponenziale n=2')
plt.grid(True)


# symmetric log
plt.subplot(223)
plt.bar(be10,ve10)
plt.title('Esponenziale n=10')
plt.grid(True)

# logit
plt.subplot(224)
plt.bar(be100,ve100)
plt.title('Esponenziale n=100')
plt.grid(True)
# Format the minor tick labels of the y-axis into empty strings with
# `NullFormatter`, to avoid cumbering the axis with too many labels.
plt.gca().yaxis.set_minor_formatter(NullFormatter())
# Adjust the subplot layout, because the logit one may take more space
# than usual, due to y-tick labels like "1 - 10^{-3}"
plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)


##########################################################lorenziana########################################

plt.figure(3)

# linear
plt.subplot(221)
plt.bar(bl1,vl1)
plt.title('Lorenziana n=1')
plt.grid(True)


# log
plt.subplot(222)
plt.bar(bl2,vl2)
plt.title('Lorenziana n=2')
plt.grid(True)


# symmetric log
plt.subplot(223)
plt.bar(bl10,vl10)
plt.title('Lorenziana n=10')
plt.grid(True)

# logit
plt.subplot(224)
plt.bar(bl100,vl100)
plt.title('Lorenziana n=100')
plt.grid(True)
# Format the minor tick labels of the y-axis into empty strings with
# `NullFormatter`, to avoid cumbering the axis with too many labels.
plt.gca().yaxis.set_minor_formatter(NullFormatter())
# Adjust the subplot layout, because the logit one may take more space
# than usual, due to y-tick labels like "1 - 10^{-3}"
plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,
                    wspace=0.35)




plt.show()


