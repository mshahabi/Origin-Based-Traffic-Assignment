
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
### GETTING LINK and TRAVELTIME DATA
Param, Links, LinkNum, NodeSet, From, To = getParam(data_path)
od_pair_demand, od_pair_list, Origin_set, Destination_set, qrs = getDemand(od_pair_data)
### DEFINING the CAP RATIO in BPR FUNCT
bc =Param['B']/Param['Capacity (veh/h)']**4 
###DEFING the FLOW ORIGINATED FROM EVERY NODE:
    
 

model = ConcreteModel()

model.nodes   = Set(initialize=NodeSet)
model.arcs    = Set(initialize=Links)
model.origins = Set(initialize=Origin_set)
model.destinations = Set(initialize=Destination_set)
model.x       = Var(model.arcs*model.origins, within=NonNegativeReals)

def obj_rule(model):
    first_term  = sum(Param['Free Flow Time (min)'][e]*(sum(model.x[e,o] for o in model.origins)) for e in model.arcs)
    second_term = sum(Param['Free Flow Time (min)'][e]*bc[e]*(sum(model.x[e,o] for o in model.origins))**5/5 for e in model.arcs)
    return first_term+second_term

model.obj     = Objective(rule=obj_rule,sense=minimize)

def flow_balance_rule(model,n,o):
    pred = [];succ=[]
    for arc in model.arcs:
        if arc[1] == n:
            pred.append(arc[0])
        if arc[0] == n:
            succ.append(arc[1])
    return sum(model.x[n,p,o] for p in pred)-sum(model.x[s,n,o] for s in succ) == qrs[o-1,n-1]    
    
model.flowbal = Constraint(model.nodes,model.origins, rule=flow_balance_rule)

def ocon(model,n,p,s,o):
    if n==o:
        if s==n:
            return model.x[p,s,o]==0
        else: return Constraint.Skip 
    else: return Constraint.Skip    
#model.flowbal.display()
model.ocon= Constraint(model.nodes,model.arcs,model.origins, rule=ocon)

instance = model.create_instance()

opt = SolverFactory("baron")
## SOLVER SPECIFIC GAP##
##########################
#ATTN: EACH SOLVER HAS ITS OWN OPTION SYNTAX, IF NOT SET PROPERLY ERROR WOULD BE RETURNED
#opt.options["mipgap"] = 0.05
#store the results 
results = opt.solve(instance)

#print and display the solver output and the results 
print(results)
a=0
for o in model.origins:
    a=instance.x[1,2,o].value +a 
    
    
class Origin_Based_TAP:
    
    def __init__(self, data_path,od_pair_data):
        self.Param, self.Links, self.LinkNum, self.NodeSet, self.From, self.To = getParam(data_path)
        self.od_pair_demand, self.od_pair_list, self.Origin_set, self.Destination_set, self.qrs = getDemand(od_pair_data)
        self.bc =Param['B']/Param['Capacity (veh/h)']**4 
        self.Create_Model()
   
    def Create_Model(self):
        self.model = ConcreteModel()
        self.model.nodes   = Set(initialize=NodeSet)
        self.model.arcs    = Set(initialize=Links)
        self.model.origins = Set(initialize=Origin_set)
        self.model.destinations = Set(initialize=Destination_set)
        self.model.x       = Var(model.arcs*model.origins, within=NonNegativeReals)
            
        def obj_rule(model):
            first_term  = sum(self.Param['Free Flow Time (min)'][e]*(sum(self.model.x[e,o] for o in self.model.origins)) for e in self.model.arcs)
            second_term = sum(self.Param['Free Flow Time (min)'][e]*bc[e]*(sum(self.model.x[e,o] for o in self.model.origins))**5/5 for e in self.model.arcs)
            return first_term+second_term
            
        self.model.obj     = Objective(rule=obj_rule,sense=minimize)
            
        def flow_balance_rule(model,n,o):
            pred = [];succ=[]
            for arc in self.model.arcs:
                if arc[1] == n:
                    pred.append(arc[0])
                if arc[0] == n:
                    succ.append(arc[1])
            return sum(self.model.x[n,p,o] for p in pred)-sum(self.model.x[s,n,o] for s in succ) == qrs[o-1,n-1]    
                
        self.model.flowbal = Constraint(model.nodes,model.origins, rule=flow_balance_rule)
            
        def ocon(model,n,p,s,o):
            if n==o:
                if s==n:
                    return self.model.x[p,s,o]==0
                else: return Constraint.Skip 
            else: return Constraint.Skip    
            #model.flowbal.display()
        self.model.ocon= Constraint(self.model.nodes,self.model.arcs,self.model.origins, rule=ocon)
        
    def solve(self):
        instance = self.model.create_instance()
        opt = SolverFactory("baron")
        results = opt.solve(instance)
        test_flow = 0
        for o in model.origins:
            test_flow = instance.x[1,2,o].value + test_flow
            
        return test_flow    
## SOLVER SPECIFIC GAP##
##########################
#ATTN: EACH SOLVER HAS ITS OWN OPTION SYNTAX, IF NOT SET PROPERLY ERROR WOULD BE RETURNED
#opt.options["mipgap"] = 0.05
#store the results 

            
        
a=Origin_Based_TAP(data_path,od_pair_data)
a.solve()