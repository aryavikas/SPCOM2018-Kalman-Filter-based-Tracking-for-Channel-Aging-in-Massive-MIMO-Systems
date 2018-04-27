from pylab import *
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
import numpy as np
x = [20,30,40,50,60,70]

datal1=[12.9557,13.8327,13.9096,14.0589,14.2832,14.4532]
datal3=[13.0853,14.0242,14.3613,14.4707,14.6653,14.7782]
datal9=[15.025,16.6993,17.2367,18.204,18.2042,18.2996]
datal15=[25.3013,28.0978,30.5899,33.0768,32.7496,33.0706]



plt.hold(True)
line1, = plt.plot(x,datal1, 'r^-', label="$m = 1$", lw=3, ms=9.0)
line2, = plt.plot(x,datal3, 'mo-', label="$m = 3$", lw=3, ms=9.0)
line3, = plt.plot(x,datal9, 'g*-', label="$m = 9$", lw=3, ms=12.0)
line4, = plt.plot(x,datal15, 'bv-', label="$m = 15$", lw=3, ms=9.0)
first_legend = plt.legend(handles=[line1,line2,line3,line4], loc='best')

plt.xlabel('Number of antennas')
plt.ylabel('Average Sum Rate [bps/Hz]')

plt.grid()
plt.show()




