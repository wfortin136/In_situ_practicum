import numpy as np
import pylab as pl

data = np.loadtxt('run_results.txt')

pl.title('Runtime vs Nodes')

pl.plot(data[:,2], data[:,14], 'ro')
pl.xlabel('Number of Nodes')
pl.ylabel('Runtime (sec)')

pl.yscale('log')

pl.show()

