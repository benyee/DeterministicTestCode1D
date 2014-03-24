from data_set import data_set
from numpy import array,cosh,sinh
import math

def diff_sol(x):
    sqrt3 = math.sqrt(3)
    return (1 - 2*cosh(sqrt3*x)/(2*cosh(20*sqrt3)+4*sinh(sqrt3*20)/sqrt3))

X = 20
eps = (0.03,0.1,0.25,0.5,1,)
dx = (4,2,1,0.5,0.25,0.01,)
alpha = ("new","SC","DD")

data_list = data_set(X,dx,alpha,eps,path_base="../OutputFiles/Run12/output_12_",dx_eps_scale = 1)
data_list.plot("../OutputFiles/Run12/",show = 0,fun=diff_sol,funlabel="Diff. Soln.")
