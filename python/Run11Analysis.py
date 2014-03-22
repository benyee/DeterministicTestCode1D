from data_read import data_read
from numpy import array

X = 20
dx = (4,2,1,0.5,0.25,0.01,)
eps = (0.03,0.1,0.25,0.5,1,)
alpha = ("DD","SC","FV")

for j in range(0,len(dx)):
    for m in range(0,len(alpha)):
        for k in range(0,len(eps)):
            print dx[j],alpha[m],eps[k]