from data_set import data_set

X = 20
dx = (4,2,1,0.5,0.25,0.01,)
eps = (0.03,0.1,0.25,0.5,1,)
alpha = ("DD","SC","new")

data_list = data_set(X,dx,alpha,eps)
data_list.plot('./figures/',show = 0)
