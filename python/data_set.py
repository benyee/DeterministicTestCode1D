from data_object import data_object
from matplotlib import pyplot
from numpy import array, where

class data_set:
    symbols = ['s','*','o','x','d','1','p','2','3','4']
    linestyles = ['-',':','--','-.']
    colors = 'bgrkmc'
    
    #dx and eps inputs should be tuples
    def __init__(self,X,dx,alpha,eps = -1):
        self.data_list = []
        self.X = X
        self.dx = array(dx)
        self.alpha = alpha
        self.eps = eps
        for j in range(0,len(dx)):
            for m in range(0,len(alpha)):
                for k in range(0,len(eps)):
                    path = "../OutputFiles/Run11/output_11_"+alpha[m]+"_"+str(j)+"_"+str(j)+"_0_.txt"
                    
                    print "Reading in "+path
                    temp = data_object(path,X,dx[j],alpha[m],eps[k])
                    temp.data_read()
                    
                    self.data_list.append(temp)

    def plot(self,path = "../figures/",color = 0,show = 1,mode = 0):
        if self.eps == -1:
            return
        
        #Set of handles and labels for legend:
#        handles_list1 = [[]]*len(self.dx)
#        labels_list1 = [[]]*len(self.dx)
#        handles_list2 = [[]]*len(self.dx)
#        labels_list2 = [[]]*len(self.dx)

        for data_obj in self.data_list:
            #Figure out which plot it should go on:
            fignum = where(data_obj.dx==self.dx)
            fignum = fignum[0][0]+1
            
            #Figure out how to annotate the plot:
            plotmod = ''
            eps_num = 0
            if color:
                eps_num = where(data_obj.eps == array(self.eps))
            else:
                eps_num = where(data_obj.eps==array(self.eps))
            eps_num = eps_num[0][0]
            plotmod = plotmod+data_set.symbols[eps_num]
            
            linesty_num = where(data_obj.alpha==array(self.alpha))
            linesty_num = linesty_num[0][0]
            plotmod = plotmod+data_set.linestyles[linesty_num]
            
            #Plot in a way to get the legend entries:
            pyplot.figure(fignum)
            handles1 = []
            labels1 = []
            handles2 = []
            labels2 = []
            if linesty_num == 0:
                temp = pyplot.plot(data_obj.x,data_obj.phi,plotmod)
                handles1.append(temp)
                labels1.append("$\epsilon = "+str(self.eps[eps_num])+"$")
            elif eps_num == 0:
                temp = pyplot.plot(data_obj.x,data_obj.phi,plotmod)
                handles2.append(temp)
                labels2.append(self.alpha[linesty_num])
            else:
                pyplot.plot(data_obj.x,data_obj.phi,plotmod)
                    

        for i in range(0,len(self.dx)):
            pyplot.figure(i+1)
            pyplot.title("$\Delta x = "+str(self.dx[i])+"$")
            pyplot.xlabel("$x$")
            pyplot.ylabel("Scalar Flux")
            #handles, legends = pyplot.get_legend_handles_labels()
            pyplot.legend(loc=0)
        
        pyplot.show()
        return 0

    def __str__(self):
        out_str = '\n'
        for obj in self.data_list:
            out_str = out_str+str(obj)+'\n';
        return out_str