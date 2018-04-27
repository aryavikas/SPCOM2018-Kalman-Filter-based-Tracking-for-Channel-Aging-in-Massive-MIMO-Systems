from pylab import *
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
import numpy as np
x = [20,30,40,50,60,70]

datal1=[41.957,44.4117,46.2795,47.5331,48.9205,49.8835]
datal3=[46.5947,49.7684,52.4763,54.3212,56.52,57.8057]
datal9=[54.3871,58.247,60.7366,62.7584,64.4172,65.298]
datal15=[58.2962,62.4194,64.2804,66.765,68.9646,69.028]


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




