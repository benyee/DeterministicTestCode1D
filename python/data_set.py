from data_object import data_object
from numpy import array
from matplotlib import pyplot

class data_set:
    def __init__(self,X,dx,alpha,eps):
        self.data_list = []
        for j in range(0,len(dx)):
            for m in range(0,len(alpha)):
                for k in range(0,len(eps)):
                    path = "../OutputFiles/Run11/output_11_"+alpha[m]+"_"+str(j)+"_"+str(j)+"_0_.txt"
                    
                    print "Reading in "+path
                    temp = data_object(path,X,dx[j],alpha[m],eps[k])
                    temp.data_read()
                    
                    self.data_list.append(array(temp))

    def plot(self,path = "../figures/",mode = 0):
        pyplot.plot(self.data_list[0].x,self.data_list[0].phi)
        pyplot.show()
        return 0

    def __str__(self):
        out_str = '\n'
        for obj in self.data_list:
            out_str = out_str+str(obj)+'\n';
        return out_str