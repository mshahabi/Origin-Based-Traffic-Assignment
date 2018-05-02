import pandas as pd
import numpy as np
from scipy import sparse
def getDemand(filename):
    
    """This function returns the demand for each O-D pair in the TAP problem
    
    INPUTS:
    filename: String which the represents the name of the file from which the demand is read
    
    OUTPUTS:
    d: N-by-K Scipy sparse array where d[:,k] is the demand vector for OD pair k and where each row corresponds to a node
    """
    
    #Here we open the file
    f  = open(filename, "r")
             
        #Parse through the file line by line
    origin_set = set()
    destination_set = set()
    for line in f.readlines():
    #Get the number of O_D pairs and creates a matrix for these O-D pairs
        if "<NUMBER OF ZONES>" in line: 
              n_OD = int(line.split(">")[1])
              D = np.zeros([n_OD,n_OD])
              
              OD_pair = []
            #Check number of Origin
        elif "Origin" in line:
              index_i=int(line.split("\t")[1])
              origin_set.update([int(index_i)])
            #Fills row for Origin index_i
        else:
              index_j =line.replace(";",":").split(":")
              for k in range(0,len(index_j)-2,2):
                    j = int(index_j[k].replace(" ",""))-1
                    destination_set.update([int(index_j[k].replace(" ",""))])
                    D[index_i-1,j]=float(index_j[k+1].replace(" ",""))
                    if (index_i != j+1):
                        OD_pair.append((index_i,j+1))
    
    qrs = np.zeros([n_OD,n_OD])
    for o in origin_set:
        for d in destination_set: 
            if o==d:
                qrs[o-1,o-1] = D[o-1,:].sum()
            else:
                qrs[o-1,d-1] = -D[o-1,d-1]
    
    return D, OD_pair, origin_set, destination_set, qrs
####test 
#od_pair_data  = curdir + "\\Data\\SiouxFalls_OD.txt"
#od_pair_demand, od_pair_list, Origin_set, Destination_set, qrs = getDemand(od_pair_data)