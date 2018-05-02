# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 18:16:53 2018

@author: mshahabi
"""
import pandas as pd
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition

import networkx as nx
import matplotlib.pyplot as plt  
import os
import time

from getPartitions import Partitions
from getParam import getParam
from getDemand import getDemand

import matplotlib.pyplot as plt
import pickle 

curdir = os.getcwd()
data_path = curdir + "\\Data\\SiouxFalls_net.txt"
partition_network_path = curdir + "\\Data\\SFOutfile.csv"
node_pos_path = curdir + "\\Data\\SiouxFalls_Nodes_Pos.txt"
od_pair_data  = curdir + "\\Data\\SiouxFalls_OD.txt"

clust_num = 4
node_per_clust = 50

netClustring= Partitions(data_path,clust_num,node_per_clust,node_pos_path)
netClustring.solve()

from origin_based_assignment_1 import Origin_Based_TAP


class cluster_assignement(Origin_Based_TAP):
    def __init__(self,data_path,od_pair_data,cluster_number,total_clusters,x_iter,lambdaa_iter,etta):
        super(cluster_assignement,self).__init__(data_path,od_pair_data)  
        
        total_clusters = [i for i in range(1,total_clusters+1)]         
        self.lkToClust,self.clustToLk,self.clustToBrknLk,self.clustToNd,self.clustOrgn,self.clustToBrknNd = netClustring.extract_partitions()
        print("Cluster Info", self.lkToClust,self.clustToLk,self.clustToBrknLk,self.clustToNd,self.clustOrgn,self.clustToBrknNd)
        self.NodeSet = self.clustToNd[cluster_number]
        self.Links   = self.clustToLk[cluster_number]
        self.Origin_clust = self.clustOrgn[cluster_number]
        self.Broken_Nodes = self.clustToBrknNd[cluster_number]
        self.Broken_Links = self.clustToBrknLk[cluster_number]
        self.Broken_Origin = set.intersection(self.Origin_set,self.Broken_Nodes)
        self.total_clusters = total_clusters
        self.cluster_number = cluster_number
        self.lambdaa = lambdaa_iter
        self.x_bar  = x_iter 
        self.create_model(x_iter,lambdaa_iter,etta)
        self.create_obj()
        self.flow_con()
        self.solve()
#        print("hhhh",self.Broken_Nodes,self.clustToNd) 
        
    def create_model(self,x_iter,lambdaa_iter,etta):
        self.model = ConcreteModel()
        self.model.nodes         = Set(initialize=self.NodeSet)
        self.model.arcs          = Set(initialize=self.Links)
        self.model.origins       = Set(initialize=self.Origin_set)
        self.model.destinations  = Set(initialize=self.Destination_set)
        self.model.borkenLinks   = Set(initialize=self.Broken_Links )
        self.model.borkenOrigins = Set(initialize=self.Broken_Origin )
        self.model.borkenNodes   = Set(initialize=self.Broken_Nodes )
        self.model.clusters      = Set(initialize=self.total_clusters)
        self.model.clusterOrigin = Set(initialize=self.Origin_clust)
        self.model.x             = Var(self.model.arcs*self.model.origins, within=NonNegativeReals)
        self.model.x_brk         = Var(self.model.borkenLinks*self.model.origins*self.model.clusters, within=NonNegativeReals)
        self.model.etta          = Param(initialize=etta)
  
    def create_obj(self):
        
        def obj_rule(model):
            lagrange = 0
            cluster_free_flow  = sum(self.Param['Free Flow Time (min)'][e]*(sum(self.model.x[e,o] for o in self.model.origins)) for e in self.model.arcs)
            cluster_bpr        = sum(self.Param['Free Flow Time (min)'][e]*self.bc[e]*(sum(self.model.x[e,o] for o in self.model.origins))**5/5 for e in self.model.arcs)
            broken_free_flow   = sum(sum(self.Param['Free Flow Time (min)'][e]*(sum(self.model.x_brk[e,o,c] for o in self.model.origins)) for e in self.model.borkenLinks) for c in self.model.clusters if c==self.cluster_number)
            broken_bpr         = sum(sum(self.Param['Free Flow Time (min)'][e]*self.bc[e]*(sum(self.model.x_brk[e,o,c] for o in self.model.origins))**5/5 for e in self.model.borkenLinks) for c in self.model.clusters if c==self.cluster_number)
            for c in self.model.clusters:
                if c==self.cluster_number:
                    for e in self.model.borkenLinks:
                        for o in self.model.origins:
                            if e[0] in self.NodeSet:
                                   lagrange           = lagrange + self.lambdaa[e,o,c]*(self.model.x_brk[e,o,c]- self.x_bar[e,o])  
                            if e[1] in self.NodeSet:
                                   lagrange           = lagrange + self.lambdaa[e,o,c]*(-self.model.x_brk[e,o,c]+ self.x_bar[e,o])
            augmented          = sum(sum(sum(self.model.etta*(self.model.x_brk[e,o,c]- self.x_bar[e,o])**2 for e in self.model.borkenLinks) for o in self.model.origins)for c in self.model.clusters if c==self.cluster_number)
            return cluster_free_flow + cluster_bpr + broken_free_flow + broken_bpr + lagrange + augmented 
        self.model.obj     = Objective(rule=obj_rule,sense=minimize) 
        self.model.obj.pprint()
#       
        
    def flow_con(self):
        def flow_balance_rule(model,n,o):
            pred = [];succ=[];pred_brk = [];succ_brk = [];  
            if n not in self.model.borkenNodes:
                for arc in self.model.arcs:
                    if arc[1] == n:
                        pred.append(arc[0])
                    if arc[0] == n:
                        succ.append(arc[1])
                return sum(self.model.x[n,p,o] for p in succ) - sum(self.model.x[s,n,o] for s in pred) == self.qrs[o-1,n-1]        
            if n in self.model.borkenNodes:
                 for arc in self.model.arcs:
                    if arc[1] == n:
                        pred.append(arc[0])
                    if arc[0] == n:
                        succ.append(arc[1])
                 for arc_brk in self.model.borkenLinks :
                    if arc_brk[1] == n:
                        pred_brk.append(arc_brk[0])
                      
                    if arc_brk[0] == n:
                        succ_brk.append(arc_brk[1])
#                print("aaa",succ_brk,pred_brk)         
                 return sum(self.model.x[n,p,o] for p in succ) + sum(sum(self.model.x_brk[n,p,o,c] for p in succ_brk)for c in self.model.clusters if c==self.cluster_number)-sum(self.model.x[s,n,o] for s in pred)-sum(sum(self.model.x_brk[s,n,o,c] for s in pred_brk)for c in self.model.clusters if c==self.cluster_number) == self.qrs[o-1,n-1]    
            
        self.model.flowbal = Constraint(self.model.nodes,self.model.origins, rule=flow_balance_rule)     
#        self.model.flowbal.pprint()
#        
    def solve(self):
        start_time = time.time()
        opt = SolverFactory("IPOPT")
#        opt.options["nlpgap"] = 0
        results = opt.solve(self.model)
#        print(results)
        print("solving each partition took --- %s seconds ---" % (time.time() - start_time))    
        x_bar = {}
        x_clust = {}
        test_flow_1 = {}
        for e in  self.model.borkenLinks:test_flow_1[e] = 0 
        for e in  self.model.borkenLinks:
            for o in self.model.origins:
                for c in self.model.clusters:
                    if c==self.cluster_number:
                        x_bar[e,o,c]=self.model.x_brk[e,o,c].value
                        print("x_bar",e,o,c,x_bar[e,o,c])     
                        test_flow_1[e] = self.model.x_brk[e,o,c].value + test_flow_1[e]      
        for e in  self.model.arcs:
            for o in self.model.origins: 
                x_clust[e,o]=self.model.x[e,o].value
#                print("x clust",e,o,x_clust[e,o])        
                
        self.x_c   = x_clust  
        self.x_bar = x_bar
        self.test_flow_1 = test_flow_1
        
class Pricing:
    def __init__(self,total_cluster_number,x_bar_,lambdaa_,etta_):
        self.x_bar_ =x_bar_
        self.lambdaa_ = lambdaa_
        self.total_cluster_number = total_cluster_number
        self.total_clusters = set([i for i in range(1, self.total_cluster_number+1)])
        self.get_optimal_brk_x(x_bar_,lambdaa_,etta_)
                         
        
    def get_optimal_brk_x(self,x_bar_,lambdaa_,etta_):
        x_brk_bar = {}
        
        lambdaa_bar = {}
        origin_set_pricing = set() 
        link_set_pricing = []
        link_set = []
        ClustToBrknLnks = {}
        x_clust = {}
        total_system_travel_time = 0
        node_clust = {}
        for i in self.total_clusters:ClustToBrknLnks[i] = set()
        for i in self.total_clusters:node_clust[i] = set()
        for clust in range(1,self.total_cluster_number+1):
            cluster_solution = cluster_assignement(data_path,od_pair_data,clust,self.total_cluster_number,self.x_bar_,self.lambdaa_,etta_)
            total_system_travel_time = total_system_travel_time + cluster_solution.model.obj.expr()
            x_brk_bar.update(cluster_solution.x_bar)
            x_clust.update(cluster_solution.x_c)
            lambdaa_bar.update(cluster_solution.lambdaa)
            etta_bar = cluster_solution.model.etta
            origin_set_pricing.update(cluster_solution.model.origins)
            link_set_pricing.append(cluster_solution.Broken_Links)
            ClustToBrknLnks[clust] = cluster_solution.Broken_Links
            node_clust[clust]  = cluster_solution.NodeSet
        for items in link_set_pricing:
            for item in items:
                if item not in link_set:
                    link_set.append(item)
        self.x_brk_bar = x_brk_bar
#        print(self.x_brk_bar)
        self.lambdaa_bar = lambdaa_bar
        self.etta_bar = etta_bar
        self.origin_set_pricing = origin_set_pricing
        self.link_set = link_set
        self.ClustToBrknLnks = ClustToBrknLnks
        self.cluster_solution = cluster_solution
        self.x_clust = x_clust
        self.total_system_travel_time = total_system_travel_time
        self.Node = node_clust
        
    def create_pricing_model(self):
        self.model = ConcreteModel()
        self.model.arcs          = Set(initialize=self.link_set)
        self.model.origins       = Set(initialize=self.origin_set_pricing)
        self.model.x_p           = Var(self.model.arcs*self.model.origins, within=NonNegativeReals)
  
    def create_obj(self):
        lagrange = {}
        augmented_lagrange = {}
        for c in self.total_clusters: lagrange[c] = 0;augmented_lagrange[c] = 0 
        def obj_rule(model):            
            for c in self.total_clusters: 
                for e in self.ClustToBrknLnks[c]: 
                    for o in self.model.origins: 
                        if e[0] in self.Node[c]:
                            print("here",self.Node[c])
                            lagrange[c]      =  lagrange[c] + self.lambdaa_bar[e,o,c]*(self.x_brk_bar[e,o,c]-self.model.x_p[e,o])
                        else:
                            print("there",self.Node[c])
                            lagrange[c]      =  lagrange[c] + self.lambdaa_bar[e,o,c]*(-self.x_brk_bar[e,o,c]+self.model.x_p[e,o])
            
                        augmented_lagrange[c]  = augmented_lagrange[c] + self.etta_bar*(self.x_brk_bar[e,o,c]-self.model.x_p[e,o])**2
                        obj_value = sum(lagrange[c]+augmented_lagrange[c] for c in self.total_clusters)
            return  obj_value   
        self.model.obj     = Objective(rule=obj_rule,sense=minimize) 
        self.model.obj.pprint()      
        
    def solve_pricing_problem(self):
        start_time = time.time()
        self.create_pricing_model()
        self.create_obj()
        opt = SolverFactory("IPOPT")
#        opt.options["nlpgap"] = 0
        print("solving prcing problem took --- %s seconds ---" % (time.time() - start_time)) 
        self.results = opt.solve(self.model)
        link_flow ={}
        for c in self.total_clusters:
            for e in self.ClustToBrknLnks[c]:
                for o in self.model.origins:
                    link_flow[e,o]=self.model.x_p[e,o].value
#                    print(e,o,link_flow[e,o])
        self.link_flow = link_flow
        self.primal_residule  = self.model.obj.expr()
        
    def updating_lagrangian_price(self):
        self.solve_pricing_problem()
        flow_results = self.link_flow

        for c in self.total_clusters:
            for e in self.ClustToBrknLnks[c]:
                for o in self.model.origins:
                    if e[0] in self.Node[c]:
                        self.lambdaa_bar[e,o,c] = self.lambdaa_bar[e,o,c]+ self.etta_bar*(self.x_brk_bar[e,o,c]-flow_results[e,o])
                    if e[1] in self.Node[c]:
                        self.lambdaa_bar[e,o,c] = self.lambdaa_bar[e,o,c]+ self.etta_bar*(-self.x_brk_bar[e,o,c]+flow_results[e,o])

lc,cL,cbl,cn,co,cbn = netClustring.extract_partitions()
os = set()
with open('flow.pickle', 'rb') as handle:
    b = pickle.load(handle)
for c in range(1,clust_num+1): os.update(co[c])
la={}; x={}
for c in range(1,clust_num+1):
      for ee in cbl[c]:
            for oo in os:
                 la[ee,oo,c] = 10
                 x[ee,oo] = b[ee,oo]

total_link_flow = {}
  
bl_set = set()      
for c in range(1,clust_num+1):  
     bl_set.update(cbl[c])
     

                     
                 
param, L, l2, l3, l4,l5 = getParam(data_path)
A, B, C, D, E = getDemand(od_pair_data)     
bc =param['B']/param['Capacity (veh/h)']**4 
 
tstt = []  
xx = []   
etta = 0.01   
tstt_ = 10000000000
for iteration in range(1,1000):
    
#    if iteration>1:
#        etta = 0.00755
#        0.008-3*iteration/100000
    print("etta is " , etta)     
    ADMM_Instance=Pricing(clust_num,x,la,etta)
    ADMM_Instance.updating_lagrangian_price()
    
    la   = ADMM_Instance.lambdaa_bar
    x   = ADMM_Instance.link_flow
    y   = ADMM_Instance.x_clust
    print("iteration %s has been passed"%iteration)
    print("primal residule is ", ADMM_Instance.primal_residule)
              
    for c in range(1,clust_num+1): 
        for  l in cbl[c]:
           total_link_flow[l] = 0

    for c in range(1,clust_num+1):         
         for  l in cL[c]:
           total_link_flow[l] = 0 
    
    
    for oo in C:
         for bl in  bl_set:
            total_link_flow[bl] = total_link_flow[bl] + x[bl,oo]
    
    for cc in range(1,clust_num+1): 
        for oo in C:       
            for ll in  cL[cc]:
                total_link_flow[ll] = total_link_flow[ll] + y[ll,oo]
 
           
    fft = sum(param['Free Flow Time (min)'][ee] for ee in L)
    bpr = sum(param['Free Flow Time (min)'][ee]*bc[ee]*(total_link_flow[ee])**4 for ee in L)
    tstt.append(fft+bpr)
    xx.append(iteration)
    plt.scatter(xx,tstt)
    plt.show()
    print("Total System Travel Time is ", fft+ bpr)
    if iteration < 50:
        etta = 0.8
    elif iteration > 50 and  iteration < 100:
        etta = 0.8
    elif iteration >= 100 and  iteration < 120:    
        etta = 0.8
    else:    
        etta = 0.7
    tstt_ =  fft+ bpr 
fft = sum(param['Free Flow Time (min)'][ee]*total_link_flow[ee] for ee in L)
bpr = sum(param['Free Flow Time (min)'][ee]*bc[ee]*(total_link_flow[ee])**5 for ee in L)    
print(fft+bpr)
for ee in L:
 print(ee,total_link_flow[ee],param['Free Flow Time (min)'][ee]*total_link_flow[ee]+10*param['Free Flow Time (min)'][ee]*bc[ee]*(total_link_flow[ee])**5)    
#a1.origin_set_pricing
#a1.x_brk_bar[(5,9),1,1]
#a1.lambdaa_bar[5,9,1,1]
#a1.ClustToBrknLnks 
#b=a1.link_flow
#a=cluster_assignement(data_path,od_pair_data,2,3)
#a.model.x_brk.display()
#a.model.obj.display()
#a.model.pprint()
##a.model.lambdaa.display()
##a.create_constraint()
##a.arcs
#
#a.model.origins.display()
#a.model.borkenOrigins.display()
#a.Broken_Links
#a.Broken_Origin
#a.create_model()
#a.create_obj()
#a.flow_con()
#a.model.flowbal[1,2]
#d=a.solve()
#a.model.lambdaa.display()
#a.model.clusters
#a.clustOrigin
#a.Destination_set