from data_set import data_set
from numpy import array,cosh,sinh
import math

def diff_sol(x):
    return 3.0 * ( 400 - x * x ) / 2 + 40

X = 20
eps = (0.03,0.1,0.25,0.5,1,)
dx = (4,2,1,0.5,0.25,0.01,)
alpha = ("new","SC")

data_list = data_set(X,dx,alpha,eps,path_base="../OutputFiles/Run14/output_soln_14_",dx_eps_scale = 1)
data_list.plot("../OutputFiles/Run14/",show = 0,fun=diff_sol,funlabel="Diff. Soln.")
