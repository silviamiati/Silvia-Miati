import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
import subprocess
from shutil import *
from glob import glob


y,x = np.loadtxt("config.final",unpack=True)
#T=np.arange(10,200,100,dtype=int)


with open('input.dat', 'r') as file:
  data = file.readlines()

for i in np.arange(15,0.01,-0.1):
  data[0] = str(i) + "\n" 
  with open('input.dat', 'w') as file:
    file.writelines(data)
  cmd = "./ese10"
  value = subprocess.call(cmd, shell = True) 

'''

def plotTSP(path, points, num_iters=1):

    """
    path: List of lists with the different orders in which the nodes are visited
    points: coordinates for the different nodes
    num_iters: number of paths that are in the path list
    
    """

    # Unpack the primary TSP path and transform it into a list of ordered 
    # coordinates

    x = []; y = []
    for i in paths[0]:
       x.append(points[i][0])
       y.append(points[i][1])
    
    plt.plot(x, y, 'co')

    # Set a scale for the arrow heads (there should be a reasonable default for this, WTF?)
    a_scale = float(max(x))/float(30)

    # Draw the older paths, if provided
    if num_iters > 1:

        for i in range(1, num_iters):

            # Transform the old paths into a list of coordinates
            xi = []; yi = [];
            for j in paths[i]:
                xi.append(points[j][0])
                yi.append(points[j][1])

            plt.arrow(xi[-1], yi[-1], (xi[0] - xi[-1]), (yi[0] - yi[-1]), 
                    head_width = a_scale, color = 'r', 
                    length_includes_head = True, ls = 'dashed',
                    width = 0.001/float(num_iters))
            for i in range(0, len(x) - 1):
                plt.arrow(xi[i], yi[i], (xi[i+1] - xi[i]), (yi[i+1] - yi[i]),
                        head_width = a_scale, color = 'r', length_includes_head = True,
                        ls = 'dashed', width = 0.001/float(num_iters))

    # Draw the primary path for the TSP problem
    plt.arrow(x[-1], y[-1], (x[0] - x[-1]), (y[0] - y[-1]), head_width = a_scale, 
            color ='g', length_includes_head=True)
    for i in range(0,len(x)-1):
        plt.arrow(x[i], y[i], (x[i+1] - x[i]), (y[i+1] - y[i]), head_width = a_scale,
                color = 'g', length_includes_head = True)
		#plt.annotate('Start', (x[-1], y[-1]))
    #Set axis too slitghtly larger than the set of x and y
    plt.xlim(min(x)*1.1, max(x)*1.1)
    plt.ylim(min(y)*1.1, max(y)*1.1)
			
    plt.show()



if __name__ == '__main__':
    # Run an example
    
    # Create a randomn list of coordinates, pack them into a list

    points = []
    for i in range(0, len(x)):
        points.append((x[i], y[i]))

    # Create two paths, teh second with two values swapped to simulate a 2-OPT
    # Local Search operation
    path=np.arange(30)
    # Pack the paths into a list
    paths = [path]
    
    # Run the function
    plotTSP(path, points, 1)












y,x = np.loadtxt("rw.dat",unpack=True)
y1,x1 = np.loadtxt("rwlibero.dat",unpack=True)


ax1 = plt.subplot(211)
ax1.plot(x,y,linewidth= 1, label='risultati1', lw=2)
#ax1.errorbar(x,y,erry)
plt.xlabel('$passi$')
plt.ylabel('$<|{r}_{N}|>_{RW}$')
plt.title('Random Walk reticolo')

ax2 = plt.subplot(212)
ax2.plot(x1,y1,linewidth= 1, label='risultati2',lw=30)
#ax2.errorbar(x1,y1,erry1)
plt.xlabel('$passi$')
plt.ylabel('$<|{r}_{N}|>_{RW}$')
plt.title('Random Walk nel continuo')

plt.subplots_adjust(hspace = 0.5)

plt.show()'''
