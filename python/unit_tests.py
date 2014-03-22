from data_read import data_read
from numpy import array
import math

X_DATA = []
PHI_DATA = []

X = 20
dx = (4.,0.01,)
eps = (0.1,1.,)
alpha = ("DD","SC","new")
    
for j in range(0,len(dx)):
    for m in range(0,len(alpha)):
        for k in range(0,len(eps)):
            path = "../OutputFiles/Run11/output_11_"+alpha[m]+"_"+str(j)+"_"+str(j)+"_0_.txt"
            
            #Calculate where we should start looking:
            num_lines = sum(1 for line in open(path,'r'))
            num_xvals = math.ceil(X/dx[j]*2)
            startline = max(1,num_lines - num_xvals*2-30)
            
            data_read(path,int(startline))