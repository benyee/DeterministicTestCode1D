from data_set import data_set

X = 20
dx = (4.,0.01,)
eps = (0.1,1.,)
alpha = ("DD","SC","new")

data_list = data_set(X,dx,alpha,eps)
data_list.plot('../figures/unit_test')
print data_list