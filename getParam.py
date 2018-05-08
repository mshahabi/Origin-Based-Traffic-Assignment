import pandas as pd
import numpy as np
from scipy import sparse
def getParam(filename):
    """This function returns the parameters for the TAP problem
    
    INPUTS:
    filename: (relative) Path of the file from which the parameters for the problem are read.
               This is a string
    
    OUTPUTS:
    param: Dictionary with the following
        c: list of link capacities, in the same order as A
        t_0: list of link free flow times
        alpha: list of link "Power" (exponent in BPR function)
        beta: list of link "B" (constant in BPR function) 
    A: list of tuples A[k]=(i,j) where i (j) is the head (tail) of link k
    M: Sparse node incidence matrix M
    N: number of nodes
    L = number of links
    """
    #Reading the parameters from the file
    Param = pd.read_csv(filename,  sep = '\t',usecols=['Capacity (veh/h)',	"Length (ft)",	"Free Flow Time (min)",	"B"])
    
    #Renaming the columns of the file according to our conventions
    Param_Columns = Param.columns
    print(Param_Columns)
#    Param.columns = ['C','t_0','beta','alpha']
     
    LinkList = pd.read_csv(filename,sep = '\t', usecols=['Tail', 'Head'])
    From = LinkList['Tail']  
    To   = LinkList['Head']
    NumberOfNodes_Tail = set(LinkList['Tail'])
    NumberOfNodes_Head = set(LinkList['Head']) 
    NumberOfNodes_Tail.update(NumberOfNodes_Head)
    
    
    LinkList = list(zip(LinkList['Tail'], LinkList['Head']))
    
    
    Param.index =   LinkList 
    NumberofLinks = Param.shape[0]
    return Param, LinkList, NumberofLinks, NumberOfNodes_Tail, From, To
### SETTING UP the DATA PATH and Directories
curdir = os.getcwd()
data_path = curdir + "\\Data\\SiouxFalls_net.txt"
partition_network_path = curdir + "\\Data\\SFOutfile.csv"
node_pos_path = curdir + "\\Data\\SiouxFalls_Nodes_Pos.txt"
od_pair_data  = curdir + "\\Data\\SiouxFalls_OD.txt"
    
Param, Links, LinkNum, NodeSet,From, To = getParam(data_path)