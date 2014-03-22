
def data_read(path,startline = 1):
    print "Reading "+path
    x = []
    phi = []
    with open(path,'r') as f:
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

    return x,phi