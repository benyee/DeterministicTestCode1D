from data_set import data_set
from numpy import array,cosh,sinh,loadtxt
from matplotlib import pyplot as plt
import math

data = loadtxt("../interSoln.txt");
data_len = 11;

for i in (0,50,100,200,400,599):
    plt.plot(data[i*data_len:(i+1)*data_len,0],data[i*data_len:(i+1)*data_len,1],label = "Iteration {0}".format(i))

plt.legend(loc='best')
plt.show()