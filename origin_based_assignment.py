
import os

from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition
      
from scipy import sparse

from getParam import getParam
from getDemand import getDemand
### SETTING UP the DATA PATH and Directories
curdir = os.getcwd()
data_path = curdir + "\\Data\\SiouxFalls_net.txt"
partition_network_path = curdir + "\\Data\\SFOutfile.csv"
node_pos_path = curdir + "\\Data\\SiouxFalls_Nodes_Pos.txt"
od_pair_data  = curdir + "\\Data\\SiouxFalls_OD.txt"

class OriginBasedTAP():
    
    def __init__(self, data_path,od_pair_data):
        self.Param, self.Links, self.LinkNum, self.NodeSet, self.From, self.To = getParam(data_path)
        self.od_pair_demand, self.od_pair_list, self.Origin_set, self.Destination_set, self.qrs = getDemand(od_pair_data)
        self.bc =self.Param['B']/self.Param['Capacity (veh/h)']**4 
        self.Create_Model()
   
    def Create_Model(self):
        self.model = ConcreteModel()
        self.model.nodes   = Set(initialize=self.NodeSet)
        self.model.arcs    = Set(initialize=self.Links)
        self.model.origins = Set(initialize=self.Origin_set)
        self.model.destinations = Set(initialize=self.Destination_set)
        self.model.x       = Var(self.model.arcs*self.model.origins, within=NonNegativeReals)
            
        def obj_rule(model):
            first_term  = sum(self.Param['Free Flow Time (min)'][e]*(sum(self.model.x[e,o] for o in self.model.origins)) for e in self.model.arcs)
            second_term = sum(self.Param['Free Flow Time (min)'][e]*self.bc[e]*(sum(self.model.x[e,o] for o in self.model.origins))**5/5 for e in self.model.arcs)
            return first_term+second_term
            
        self.model.obj     = Objective(rule=obj_rule,sense=minimize)
            
        def flow_balance_rule(model,n,o):
            pred = [];succ=[]
            for arc in self.model.arcs:
                if arc[1] == n:
                    pred.append(arc[0])
                if arc[0] == n:
                    succ.append(arc[1])
                if n in self.model.origins:
                    q_rs = self.qrs[o-1,n-1] 
                else:
                    q_rs= 0
            return sum(self.model.x[n,p,o] for p in succ)-sum(self.model.x[s,n,o] for s in pred) == q_rs   
                
        self.model.flowbal = Constraint(self.model.nodes,self.model.origins, rule=flow_balance_rule)
            
        def ocon(model,n,p,s,o):
            if n==o:
                if s==n:
                    return self.model.x[p,s,o]==0
                else: return Constraint.Skip 
            else: return Constraint.Skip    
            #model.flowbal.display()
        self.model.ocon= Constraint(self.model.nodes,self.model.arcs,self.model.origins, rule=ocon)
        
    def solve(self):
        
        opt = SolverFactory("IPOPT")
        results = opt.solve(self.model)
        test_flow = 0
        for o in self.model.origins:
            test_flow = self.model.x[1,117,o].value + test_flow
            
        return test_flow    


            
###TEST CODE###        
a=OriginBasedTAP(data_path,od_pair_data)
b=a.solve()