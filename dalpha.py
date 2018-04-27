#for l=6
from pylab import *
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
import numpy as np
x = [20,30,40,50,60,70]

datal1=[65.7623,69.0477,69.554,70.8325,71.0676,71.3576]
datal2=[50.3904,52.5126,54.1197,54.3704,56.0471,56.6293]
datal3=[20.813,21.591,21.888,22.334,22.634,22.99]
datal4=[9.627,10.8828,11.7528,11.7014,11.8876,12.6821]

plt.hold(True)
line1, = plt.plot(x,datal1, 'r^-', label="$f_DT_s: 0.005-0.05$", lw=3, ms=9.0)
line2, = plt.plot(x,datal2, 'mv-', label="$f_DT_s: 0.05-0.1$", lw=3, ms=9.0)
line3, = plt.plot(x,datal3, 'go-', label="$f_DT_s: 0.1-0.2$", lw=3, ms=9.0)
line4, = plt.plot(x,datal4, 'b*-', label="$f_DT_s: 0.5-0.6$", lw=3, ms=12.0)

first_legend = plt.legend(handles=[line1,line2,line3,line4], loc=7,bbox_to_anchor=(1.0,0.45))
#first_legend = plt.legend(handles=[line1,line2,line3,line4], loc=7)
plt.xlabel('Number of antennas')
plt.ylabel('Average Sum Rate [bps/Hz]')

plt.grid()
plt.show()






