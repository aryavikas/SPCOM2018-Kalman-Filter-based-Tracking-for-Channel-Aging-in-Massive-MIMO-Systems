#for l=6
from pylab import *
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
import numpy as np
x = [20,30,40,50,60,70]

datal1=[51.1389,55.241,58.0812,59.6638,61.5756,62.7718]
datal2=[44.0442,47.7625,49.8974,52.8254,54.1335,55.3828]
datal3=[32.466,34.1956,35.8139,37.0752,37.9049,38.7998]
datal4=[23.4549,25.2685,26.3975,26.7502,27.807,28.7189]


plt.hold(True)
line1, = plt.plot(x,datal1, 'r^-', label="$f_DT_s: 0.005-0.05$", lw=3, ms=9.0)
line2, = plt.plot(x,datal2, 'mo-', label="$f_DT_s: 0.05-0.1$", lw=3, ms=9.0)
line3, = plt.plot(x,datal3, 'b*-', label="$f_DT_s: 0.1-0.2$", lw=3, ms=12.0)
line4, = plt.plot(x,datal4, 'gv-', label="$f_DT_s: 0.5-0.6$", lw=3, ms=9.0)

first_legend = plt.legend(handles=[line1,line2,line3,line4], loc=9)

plt.xlabel('Number of antennas')
plt.ylabel('Average Sum Rate [bps/Hz]')

plt.grid()
plt.show()






