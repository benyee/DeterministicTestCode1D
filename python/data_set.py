from data_object import data_object
from matplotlib import pyplot
from numpy import array, where

class data_set:
    #Plot symbols:
    symbols = ['s','*','o','x','d','1','p','2','3','4']
    linestyles = ['-',':','--','-.']
    colors = 'bgrkmc'
    
    #dx and eps inputs should be tuples
    #dx_eps_scale says whether or not dx should be scaled by epsilon
    def __init__(self,X,dx,alpha,eps = -1,path_base="../OutputFiles/Run11/output_11_",dx_eps_scale = 0):
        self.data_list = []
        self.X = X
        self.dx = array(dx)
        self.alpha = alpha
        self.eps = eps
        self.dx_eps_scale = dx_eps_scale
        for j in range(0,len(dx)):
            for m in range(0,len(alpha)):
                for k in range(0,len(eps)):
                    path = path_base+alpha[m]+"_"+str(k)+"_"+str(j)+"_0_.txt"
                    
                    print "Reading in "+path
                    if dx_eps_scale and eps != -1:
                        temp = data_object(path,X,dx[j]*eps[k],alpha[m],eps[k])
                        temp.data_read()
                        
                        self.data_list.append(temp)
                    else:
                        temp = data_object(path,X,dx[j],alpha[m],eps[k])
                        temp.data_read()
                        
                        self.data_list.append(temp)

    def plot(self,path = "./figures/",color = 0,show = 1,mode = 0,fun="none",funlabel=''):
        if self.eps == -1:
            return
        
        #Set of handles and labels for legend:
        handles1 = [[0 for x in range(0,len(self.eps))] for x in range(0,len(self.dx))]
        labels1 = [[0 for x in range(0,len(self.eps))] for x in range(0,len(self.dx))]
        handles2 = [[0 for x in range(0,len(self.alpha))] for x in range(0,len(self.dx))]
        labels2 = [[0 for x in range(0,len(self.alpha))] for x in range(0,len(self.dx))]
        
        for data_obj in self.data_list:
            #Figure out which plot it should go on:
            fignum = 0
            if self.dx_eps_scale:
                fignum = where(data_obj.dx==self.dx*data_obj.eps)
                fignum = fignum[0][0]
            else:
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
                handles1[fignum][eps_num], = pyplot.plot(data_obj.x,data_obj.phi,plotmod,markevery=max(len(data_obj.x)/10,1))
                if eps_num == 0:
                    handles2[fignum][eps_num], =pyplot.plot(data_obj.x,data_obj.phi,data_set.linestyles[linesty_num]+(1-color)*'k',markevery=max(len(data_obj.x)/10,1))
                    labels2[fignum][eps_num] =self.alpha[0]
                    if fun!="none":
                        funhand, = pyplot.plot(data_obj.x,fun(data_obj.x),"k-.p",markevery=max(len(data_obj.x)/20,1))
                        if funlabel != '':
                            handles2[fignum].append(funhand)
                            labels2[fignum].append(funlabel)
            else:
                pyplot.plot(data_obj.x,data_obj.phi,plotmod,markevery=max(len(data_obj.x)/10,1))
                if eps_num == 0:
                    handles2[fignum][linesty_num], = pyplot.plot(data_obj.x[0],data_obj.phi[0],data_set.linestyles[linesty_num]+(1-color)*'k')
                    labels2[fignum][linesty_num]=self.alpha[linesty_num]
            
        #Format figures and save:
        for i in range(0,len(self.dx)):
            pyplot.figure(i+1)
            
            #Annotate:
            if self.dx_eps_scale:
                pyplot.title("$\Delta x = "+str(self.dx[i])+"$ mfp")
            else:
                pyplot.title("$\Delta x = "+str(self.dx[i])+"$")
            pyplot.xlabel("x")
            pyplot.ylabel("Scalar Flux")
            pyplot.legend(handles1[i]+handles2[i],labels1[i]+labels2[i],loc=0,ncol=2)
            
            #Rescale y-axis on plot if necessary because it gets weird for some reason:
            gca = pyplot.gca()
            curr_ylim = gca.get_ylim()
            curr_data = gca.get_lines()
            curr_data = curr_data[0].get_ydata()
            curr_data = float(max(curr_data))
            if (curr_ylim[1]-curr_ylim[0])/(curr_data-curr_ylim[0]) < 1.1:
                gca.set_ylim((curr_ylim[0],(curr_data - curr_ylim[0])*0.201+curr_ylim[1]))
            
            #Save:
            pyplot.savefig(path+"dx_"+str(i)+".png")
        
        if show:
            pyplot.show()
        return 0

    def __str__(self):
        out_str = '\n'
        for obj in self.data_list:
            out_str = out_str+str(obj)+'\n';
        return out_str