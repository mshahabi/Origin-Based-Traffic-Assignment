
import os

from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition

from scipy import sparse

from getParam import getParam
from getDemand import getDemand
import pickle

### SETTING UP the DATA PATH and Directories
curdir = os.getcwd()
data_path = curdir + "\\Data\\SiouxFalls_net.txt"
partition_network_path = curdir + "\\Data\\SFOutfile.csv"
node_pos_path = curdir + "\\Data\\SiouxFalls_Nodes_Pos.txt"
od_pair_data  = curdir + "\\Data\\SiouxFalls_OD.txt"

class Origin_Based_TAP:
    
    def __init__(self, data_path,od_pair_data):
        self.Param, self.Links, self.LinkNum, self.NodeSet, self.From, self.To = getParam(data_path)
        self.od_pair_demand, self.od_pair_list, self.Origin_set, self.Destination_set, self.qrs = getDemand(od_pair_data)
        self.bc =self.Param['B']/self.Param['Capacity (veh/h)']**4 
        
    def create_model(self):
        self.model = ConcreteModel()
        self.model.nodes   = Set(initialize=self.NodeSet)
        self.model.arcs    = Set(initialize=self.Links)
        self.model.origins = Set(initialize=self.Origin_set)
        self.model.destinations = Set(initialize=self.Destination_set)
        self.model.x       = Var(self.model.arcs*self.model.origins, within=NonNegativeReals)
#        print(self.model.x)
    def create_obj(self):
        def obj_rule(model):
            first_term  = sum(self.Param['Free Flow Time (min)'][e]*(sum(self.model.x[e,o] for o in self.model.origins)) for e in self.model.arcs)
            second_term = sum(self.Param['Free Flow Time (min)'][e]*self.bc[e]*(sum(self.model.x[e,o] for o in self.model.origins))**5/5 for e in self.model.arcs)
            return first_term+second_term
        self.model.obj     = Objective(rule=obj_rule,sense=minimize)    
    
    def flow_con(self):
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
        self.model.flowbal.pprint()
     
    def con(self):
        def ocon(model,n,p,s,o):
            if n==o:
                if s==n:
                    return self.model.x[p,s,o]==0
                else: return Constraint.Skip 
            else: return Constraint.Skip    
        self.model.ocon= Constraint(self.model.nodes,self.model.arcs,self.model.origins, rule=ocon)
    
    def solve(self):
        opt = SolverFactory("baron")
        opt.options["nlpgap"] = 0
        results = opt.solve(self.model)
        print(results)
        test_flow = {}  
        link_flow = {}
        tstt = 0
        for e in self.model.arcs:
            link_flow[e] = 0
            for o in self.model.origins:
                test_flow[e,o] = self.model.x[e,o].value 
                print(e,o,test_flow[e,o])
                link_flow[e]   =  link_flow[e] + self.model.x[e,o].value
        for e in self.model.arcs:        
            first_term  = self.Param['Free Flow Time (min)'][e] 
            second_term = self.Param['Free Flow Time (min)'][e]*self.bc[e]*(link_flow[e])**4
            tstt = tstt + first_term + second_term
            print(e,link_flow[e],self.bc[e],first_term+second_term)
        self.test_flow = test_flow         
        self.link_flow = link_flow 
        self.tstt = tstt
        with open('flow.pickle', 'wb') as handle:
            pickle.dump(test_flow, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return test_flow 


#            
####TEST CODE###        
a=Origin_Based_TAP(data_path,od_pair_data)
a.qrs
a.bc
a.Links
a.Origin_set
a.create_model()

a.create_obj()
a.flow_con()
a.con()
b=a.solve()
print(a.model.obj.expr(),a.tstt)
print(a.link_flow)