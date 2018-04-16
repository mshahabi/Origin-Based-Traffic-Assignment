# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 19:43:54 2017
min x1 + 2x2
s.t. 3x1 + 4x2 ≥ 1
2x1 + 5x2 ≥ 2
x1,x2 ≥ 0
@author: mehrdad
"""
""" Install the pyomo package """
#conda install -c conda-forge pyomo 
#install IPOPT a linear and mixed integer solver
#conda install -q -y --channel cachemeorg ipopt_bin
#install glpk a linear and mixed integer solver
#conda install -c conda-forge glpk 
import os
curdir = os.getcwd()
data_path = curdir + "\\Data\\SiouxFalls_net.txt"
partition_network_path = curdir + "\\Data\\SFOutfile.csv"
node_pos_path = curdir + "\\Data\\SiouxFalls_Nodes_Pos.txt"

import networkx as nx
import matplotlib.pyplot as plt    
from getParam import getParam
from getPartitions import getPartitions, node_position

#def getNetworkParameters(data_path):
Param,Links,LinkNum, NodeSet,From,To = getParam(data_path)

Links_Partitions_List = getPartitions(partition_network_path)

node_pos = node_position(node_pos_path )

    
G = nx.Graph()
G.add_edges_from(Links)
for i in node_pos.Node:
    G.add_node(i,pos=(node_pos.X[i-1],node_pos.Y[i-1]))

pos=nx.get_node_attributes(G,'pos') 
   
plt.figure(1,figsize=(12,12)) 
nx.draw_networkx_nodes(G, pos, node_size = 5)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos,  edge_color='r')
plt.show()


G = nx.Graph()
for i in node_pos.Node:
    G.add_node(i,pos=(node_pos.X[i-1],node_pos.Y[i-1]))
    
pos=nx.get_node_attributes(G,'pos') 
G.add_edges_from(Links_Partitions_List[2])

plt.figure(1,figsize=(12,12)) 
nx.draw_networkx_nodes(G, pos, node_size = 5)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos,  edge_color='r')
plt.show()

G = nx.Graph()
for i in node_pos.Node:
    G.add_node(i,pos=(node_pos.X[i-1],node_pos.Y[i-1]))
    
pos=nx.get_node_attributes(G,'pos') 
G.add_edges_from(Links_Partitions_List[1])

plt.figure(1,figsize=(12,12)) 
nx.draw_networkx_nodes(G, pos, node_size = 5)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos,  edge_color='r')
plt.show()