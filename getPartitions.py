import pandas as pd
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition
from getParam import getParam
from getDemand import getDemand
import networkx as nx
import matplotlib.pyplot as plt  
import os
### SETTING UP the DATA PATH and Directories
curdir = os.getcwd()
data_path = curdir + "\\Data\\SiouxFalls_net.txt"
partition_network_path = curdir + "\\Data\\SFOutfile.csv"
node_pos_path = curdir + "\\Data\\SiouxFalls_Nodes_Pos.txt"
od_pair_data  = curdir + "\\Data\\SiouxFalls_OD.txt"

class Partitions:
    
    def __init__(self,data_path,num_cluster,max_cardinality,node_position_data):     
        self.Param, self.Links, self.LinkNum, self.NodeSet, self.From, self.To = getParam(data_path)
        self.od_pair_demand, self.od_pair_list, self.Origin_set, self.Destination_set, self.qrs = getDemand(od_pair_data)
        self.num_cluster = range(1,num_cluster+1)
        self.max_cardinality = max_cardinality
        self.node_position_data = node_position_data
        self.Create_Model()   
        
    def Create_Model(self):
        self.model = ConcreteModel()
        self.model.nodes   = Set(initialize=self.NodeSet)
        self.model.arcs    = Set(initialize=self.Links)
        self.model.num_cluster = Set(initialize=self.num_cluster)
        self.model.max_cardinality = Param(initialize=self.max_cardinality)
        self.model.x       = Var(self.model.nodes*self.model.num_cluster, within=Binary)
        self.model.n_k     = Var(self.model.num_cluster, within=NonNegativeReals) 
          
        def obj_rule(model):
            first_term  = sum(sum(self.model.x[i,k]*self.model.x[j,k] for i, j in self.model.arcs) for k in self.model.num_cluster) 
            return first_term
                    
        self.model.obj     = Objective(rule=obj_rule,sense=maximize)
                    
        def assignment_rule(model,n):
            return sum(self.model.x[n,k] for k in self.model.num_cluster) == 1   
                        
        self.model.assignment = Constraint(self.model.nodes,rule=assignment_rule)
        
        def cluster_rule(model,k):
            return sum(self.model.x[n,k] for n in self.model.nodes) >= self.model.n_k[k]   
                        
        self.model.cluster = Constraint(self.model.num_cluster,rule=cluster_rule)
                    
        def cardinality_rule(model,k):
            return self.model.max_cardinality <= self.model.n_k[k]  

        self.model.cardinality = Constraint(self.model.num_cluster, rule=cardinality_rule)
        
    def solve(self):
        opt = SolverFactory("baron")
        self.results = opt.solve(self.model)
        return self.results 
    
    def extract_partitions(self):
        LinkToCluster = {}
        ClusterToLink = {}
        ClusterToBrokenLinks = {}
        ClusterToBrokenNodes = {}
        ClusterToNode = {}
        ClusterOrigins = {}
        
        for i in self.num_cluster:ClusterToNode[i] = set()
        for i in self.num_cluster:ClusterOrigins[i] = set()
        for i in self.num_cluster:ClusterToBrokenNodes[i] = set()
        for i in self.num_cluster:ClusterToBrokenLinks[i] = []

        for i in self.num_cluster:ClusterToLink[i] = []
        for p,s in self.model.arcs:
            for k in self.model.num_cluster:
               if self.model.x[p,k].value>0 and self.model.x[s,k].value>0:
                   LinkToCluster[p,s] = k
                   ClusterToLink[k].append((p,s))
                   ClusterToNode[k].update([p,s])
                   if p in self.Origin_set:
                       ClusterOrigins[k].update([p])
               elif self.model.x[p,k].value==1 and self.model.x[s,k].value==0:                 
                   ClusterToBrokenLinks[k].append((p,s))
                   if p in self.Origin_set:
                       ClusterOrigins[k].update([p])
                   ClusterToBrokenNodes[k].update([p])
               elif self.model.x[p,k].value==0 and self.model.x[s,k].value==1:
                   ClusterToBrokenLinks[k].append((p,s))  
                   ClusterToBrokenNodes[k].update([s])
                   
                   
        return LinkToCluster,ClusterToLink,ClusterToBrokenLinks,ClusterToNode,ClusterOrigins,ClusterToBrokenNodes
      

    def node_positions(self,filename):
        pos = pd.read_csv(filename,sep = '\t', usecols=['Node', 'X', 'Y']) 
        return pos        
    
    def draw_original_graph(self):
        node_pos = self.node_positions(self.node_position_data)
        G = nx.DiGraph()
        G.add_edges_from(self.Links)
        
        for i in node_pos.Node:
            G.add_node(i,pos=(node_pos.X[i-1],node_pos.Y[i-1]))
        pos=nx.get_node_attributes(G,'pos') 
   
        plt.figure(1,figsize=(12,12)) 
        nx.draw_networkx_nodes(G, pos, node_size = 5)
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos,  edge_color='r')
        plt.show()   

    def draw_partitions(self):
        node_pos = self.node_positions(self.node_position_data)
        LinkToCluster,ClusterToLink,ClusterToBrokenLinks,ClusterToNode,ClusterOrigins,ClusterToBrokenNodes = self.extract_partitions()
        G = nx.DiGraph()
        for i in node_pos.Node:
            G.add_node(i,pos=(node_pos.X[i-1], node_pos.Y[i-1]))
        pos=nx.get_node_attributes(G,'pos') 
        G.add_edges_from(LinkToCluster)

        plt.figure(1,figsize=(12,12)) 
        nx.draw_networkx_nodes(G, pos, node_size = 5)
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos,  edge_color='r', arrows=True)
        plt.show()

#
#
#netClustring= Partitions(data_path,2,12,node_pos_path)
#clusters = netClustring.solve()
#netClustring.draw_partitions() 
#netClustring.draw_original_graph()  
##a=Partitions(data_path,2,3,node_pos_path)                
##a.solve()
#b,d,d2,d3,d4,d5=netClustring.extract_partitions()
##a.draw_partitions()











#    def getPartitions(filename):
#        
#       Param = pd.read_csv(filename, usecols=['i', 'j', 'k'])    
#       num_of_partitions = max(set(Param['k']))
#       
#       partitionToLinks = {}
#       for i in range(1,num_of_partitions+1):
#           partitionToLinks[i] = []
#       
#       for i in range(len(Param)):
#           partitionToLinks[Param['k'][i]].append((Param['i'][i],Param['j'][i]))
#        
#       return partitionToLinks
#    
#    def node_position(filename):
#        
#        pos = pd.read_csv(filename,sep = '\t', usecols=['Node', 'X', 'Y']) 
#        
#        return pos
#        
