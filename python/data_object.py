import math
from numpy import array

class data_object:
    def __init__(self,path,X,dx,alpha,eps = -1):
        self.path = path
        self.X = X
        self.dx = dx
        self.alpha = alpha
        self.eps = eps
    
    def data_read(self):
        #Calculate where we should start looking:
        num_lines = sum(1 for line in open(self.path,'r'))
        num_xvals = math.ceil(self.X/self.dx*2)
        startline = max(1,num_lines - num_xvals*2-30)
        
        print "Reading "+self.path
        x = []
        phi = []
        with open(self.path,'r') as f:
            foundData = 0;
            foundData2 = 0;
            
            i = 0
            for line in f:# and f.readline().contains('<'):
                #Find out where to start:
                i = i+1
                if(i < startline):
                    continue
        
                        
                if foundData2:
                    try:
                        temp = line.split();
                        x.append(temp[0]);
                        phi.append(temp[1]);
                    except:
                        print "Reached end of data"
                        break;
                        
                if foundData:
                    foundData2 = 1
                
                if line[0]=='<':
                    foundData = 1

        self.x = array(x);
        self.phi = array(phi);
        return 1

    def __str__(self):
        return "<----START---->\ndx ="+str(self.dx)+"\nalpha ="+self.alpha+"\neps ="+str(self.eps)+"\nx-vector = "+str(self.x)+", \n\nphi = "+str(self.phi)+"\n<----END---->"
