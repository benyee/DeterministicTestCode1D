from data_object import data_object
from matplotlib import pyplot
from numpy import array, where

class data_set:
    #Plot symbols:
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
                    path = "../OutputFiles/Run11/output_11_"+alpha[m]+"_"+str(k)+"_"+str(j)+"_0_.txt"
                    
                    print "Reading in "+path
                    temp = data_object(path,X,dx[j],alpha[m],eps[k])
                    temp.data_read()
                    
                    self.data_list.append(temp)

    def plot(self,path = "./figures/",color = 0,show = 1,mode = 0):
        if self.eps == -1:
            return
        
        #Set of handles and labels for legend:
        handles1 = [['']*len(self.eps)]*len(self.dx)
        labels1 = [['']*len(self.eps)]*len(self.dx)
        handles2 = [['']*len(self.alpha)]*len(self.dx)
        labels2 = [['']*len(self.alpha)]*len(self.dx)
        
        for data_obj in self.data_list:
            #Figure out which plot it should go on:
            fignum = where(data_obj.dx==self.dx)
            fignum = fignum[0][0]
            
            #Figure out how to annotate the plot:
            plotmod = ''
            eps_num = 0
            if color:
                eps_num = where(data_obj.eps == array(self.eps))
            else:
                eps_num = where(data_obj.eps==array(self.eps))
                plotmod = plotmod+'k'
            eps_num = eps_num[0][0]
            plotmod = plotmod+data_set.symbols[eps_num]
            
            linesty_num = where(data_obj.alpha==array(self.alpha))
            linesty_num = linesty_num[0][0]
            plotmod = plotmod+data_set.linestyles[linesty_num]
            
            #Plot in a way to get the legend entries right:
            pyplot.figure(fignum+1)
            if linesty_num == 0:
                labels1[fignum][eps_num] = "$\epsilon = "+str(self.eps[eps_num])+"$"
                handles1[fignum][eps_num], = pyplot.plot(data_obj.x,data_obj.phi,plotmod)
                if eps_num == 0:
                    handles2[fignum][eps_num], =pyplot.plot(data_obj.x[0],data_obj.phi[0],data_set.linestyles[linesty_num])
                    labels2[fignum][eps_num] =self.alpha[0]
            else:
                pyplot.plot(data_obj.x,data_obj.phi,plotmod)
                if eps_num == 0:
                    handles2[fignum][linesty_num], = pyplot.plot(data_obj.x[0],data_obj.phi[0],data_set.linestyles[linesty_num])
                    labels2[fignum][linesty_num]=self.alpha[linesty_num]
                        
        #Format figures and save:
        for i in range(0,len(self.dx)):
            pyplot.figure(i+1)
            pyplot.title("$\Delta x = "+str(self.dx[i])+"$")
            pyplot.xlabel("x")
            pyplot.ylabel("Scalar Flux")
            pyplot.legend(handles1[i]+handles2[i],labels1[i]+labels2[i],loc=0)
            pyplot.savefig(path+"dx_"+str(i)+".png")
        
        if show:
            pyplot.show()
        return 0

    def __str__(self):
        out_str = '\n'
        for obj in self.data_list:
            out_str = out_str+str(obj)+'\n';
        return out_str