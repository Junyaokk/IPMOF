from algorithms.solver import SOLVER
from orders.branch import OrderBranch
from algorithms.pha import ProgressiveHedging
from gurobipy import *

class Assortment_CapConstr():
    def __init__(self, orders: OrderBranch, 
                 solver: SOLVER,
                 pha: ProgressiveHedging):
        self.orders = orders
        self.solver = solver
        self.pha = pha
        self.fdc_id = orders.fdc_id
        self.scenario_set = solver.scenario_set
        self.order_type_set = orders.order_type_set
        self.sku_set = orders.sku_set
        self.fdc_capacity = solver.fdc_capacity
        self.type_size_dict = self.pha.type_size_dict
        self.candidate_sku_set = []
        self.candidate_type_set = []
        self.type_corr_sku_set = self.pha.type_composition_by_sku
        self.A_sol_set = None
        self.final_sol_on_sku = None
    

    def solve_knapscak_problem(self):
        # solve dependent knapscak problem
        model = Model('Opt-Assortment')
        model.params.OutputFlag = 0
        Xsku = model.addVars(self.candidate_sku_set, vtype=GRB.BINARY, name="xsku")
        Xtype = model.addVars(self.candidate_type_set, vtype=GRB.BINARY, name="xtype")

        model.addConstrs((Xtype[type_id] <= Xsku[i] for type_id in self.candidate_type_set
                                                for i in self.type_corr_sku_set[type_id]))
        model.addConstr((quicksum(Xsku[i] for i in self.candidate_sku_set) <= self.fdc_capacity))
        
        reduction_cost = quicksum(self.complete_reduction_cost[type_id] * Xtype[type_id]
                                  for type_id in self.candidate_type_set)

        model.setObjective(reduction_cost, GRB.MAXIMIZE)

        model.optimize()

        A_sol_on_sku = [sku for sku in self.candidate_sku_set if Xsku[sku].x > 0]
        A_sol_on_type = [type_id for type_id in self.candidate_type_set if Xtype[type_id].x > 0]

        
        return A_sol_on_sku, A_sol_on_type
    
    
    def revise_sol_by_assortment(self, input_fdc, A_sol_on_sku, A_sol_on_type, sol_by_sku, sol_by_type):
        if sol_by_sku['F'][input_fdc] > 0:
            for sku in list(set(self.candidate_sku_set) - set(A_sol_on_sku)):
                sol_by_sku['X'][(input_fdc, sku)] = 0
                sol_by_sku['S'][(input_fdc, sku)] = 0
            for sku in A_sol_on_sku:
                sol_by_sku['S'][(input_fdc, sku)] = sum([sol_by_type[type_id]['S']
                                       for type_id in self.orders.sku_type_map_dict[sku]
                                       if type_id in A_sol_on_type])
        return sol_by_sku
    
    
    def __call__(self):
        input_fdc = self.fdc_id
        self.pha()
        sol_by_type = self.pha.solved_type_sol
        sol_k_by_sku = self.pha.aggregate_sol_on_sku(input_fdc, sol_by_type)
        # summarize the candidate set from APH's solution
        self.candidate_type_set = [type_id for type_id in self.order_type_set 
                                    if sol_by_type[type_id]['X'] == 1]
        self.candidate_sku_set = [sku for sku in self.sku_set 
                                    if sol_k_by_sku['X'][(input_fdc, sku)] == 1]
        # evaluated the saving cost for each order type
        self.pha.eva_reduction_cost_given_sol(input_fdc, sol=sol_by_type)
        self.complete_reduction_cost = self.pha.completely_reduction_set
        self.A_sol_on_sku, self.A_sol_on_type = self.solve_knapscak_problem()
        self.final_sol_on_sku = self.revise_sol_by_assortment(input_fdc, 
                                                              A_sol_on_sku=self.A_sol_on_sku,
                                                              A_sol_on_type=self.A_sol_on_type,
                                                              sol_by_sku=sol_k_by_sku,
                                                              sol_by_type=sol_by_type)