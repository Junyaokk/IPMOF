from domain.instance import Instance
from gurobipy import *
import numpy as np
import copy

class Repooling_Opt():
    def __init__(self, instance: Instance,
                 input_sol_on_sku=None,
                 input_sol_on_type=None,
                 demand_aggregation_startegy='F',
                 capacity_type='warehouse'):
        self.instance = instance
        self.input_sol_on_sku = input_sol_on_sku
        self.input_sol_on_type = input_sol_on_type
        self.capacity_type = capacity_type
        self.demand_aggregation_startegy = demand_aggregation_startegy
        self.fdc_id = self.instance.branch.fdc_id
        self.rdc_id = self.instance.branch.rdc_id
        self.warehouse_capacity = None
        if self.capacity_type == 'warehouse':
            self.delta_pool = self.instance.capmanager.CapConstr.delta_pool
            self.mu_pool = self.instance.capmanager.CapConstr.mu_pool
            self.warehouse_capacity = self.instance.solver.fdc_capacity
        elif self.capacity_type == 'mixed':
            self.delta_pool = self.instance.capmanager.CapConstr.warehouse_capconstr.delta_pool
            self.mu_pool = self.instance.capmanager.CapConstr.warehouse_capconstr.mu_pool
            self.warehouse_capacity = self.instance.solver.warehouse_capacity

    def find_opt_mu(self):
        closest_negative = min((num for num in self.delta_pool if num < 0), default=None, key=abs)
        if closest_negative is None:
            return 0
        self.opt_mu = self.mu_pool[self.delta_pool.index(closest_negative)]


    def demand_aggregation(self):
        self.aggregate_demand_dict = {}
        # Get assortment SKU set
        self.assortment_sku_set = [sku for sku in self.instance.branch.sku_set
                                if self.input_sol_on_sku['X'][self.fdc_id, sku] == 1]

        # Select type set based on demand aggregation strategy
        if self.demand_aggregation_startegy == 'F':
            self.selected_type_set = [
                type_id for type_id in self.instance.branch.order_type_set
                if (self.input_sol_on_type[type_id]['X'] == 1) or 
                   (set(self.instance.branch.type_element_pool[type_id]).issubset(self.assortment_sku_set))
            ]
        else:  # demand_aggregation_startegy == 'T'
            self.find_opt_mu()
            self.selected_type_set = [
                type_id for type_id in self.instance.branch.order_type_set
                if self.input_sol_on_type[type_id]['X'] == 1
            ]

        # Calculate transshipment coefficients
        self.transshipment_coef_dict = {
            type_id: (self.instance.branch.first_cost_structure[self.instance.branch.rdc_id][self.instance.branch.fdc_id]['coef'] *
                      self.instance.branch.type_weight_dict[type_id] +
                      (self.opt_mu * self.instance.branch.type_size_dict[type_id] if self.demand_aggregation_startegy == 'T' else 0))
            for type_id in self.selected_type_set
        }

        # Calculate second difference dictionary
        self.type_second_diff_dict = {
            type_id: dict(zip(self.instance.branch.type_region_map[type_id]['region'],
                              self.instance.branch.type_region_map[type_id]['diff'][:-1]))
            for type_id in self.selected_type_set
            if type_id in self.instance.branch.type_region_map
        }

        # Aggregate demand calculation
        for sku in self.assortment_sku_set:
            self.aggregate_demand_dict[sku] = {
                region_id: {scen_id: 0 for scen_id in self.instance.branch.scenario_set}
                for region_id in self.instance.branch.region_set
            }
            for region_id in self.instance.branch.region_set:
                for scen_id in self.instance.branch.scenario_set:
                    sum_demand = sum(
                        self.instance.branch.type_element_dict[type_id][sku] * 
                        self.instance.branch.scenario_sales_by_region[region_id][scen_id].get(type_id, 0)
                        for type_id in set(self.instance.branch.sku_type_map_dict[sku]) & set(self.selected_type_set)
                        if (self.demand_aggregation_startegy == 'F' or 
                            (self.type_second_diff_dict[type_id].get(region_id, 0) - self.transshipment_coef_dict[type_id] > 0))
                    )
                    self.aggregate_demand_dict[sku][region_id][scen_id] = sum_demand

        # Calculate selected involved size dictionary
        self.selected_involved_size_dict = {
            sku: np.mean([self.instance.branch.type_size_dict[type_id]
                          for type_id in set(self.instance.branch.sku_type_map_dict[sku]) & set(self.selected_type_set)
                          for _ in range(self.instance.source.type_sales_dict[type_id]) 
                          if type_id in self.instance.source.type_sales_dict])
            for sku in self.assortment_sku_set
        }


    def generalized_LP(self):
        model = Model('Re-pooling')
        model.params.OutputFlag = 0
        S = model.addVars(self.assortment_sku_set, vtype=GRB.CONTINUOUS, lb=0, name="supply")
        Pk = model.addVars(self.instance.branch.scenario_set, self.assortment_sku_set, self.instance.branch.region_set, lb=0, ub=1, name="pk")

        if self.capacity_type != 'assortment':
            model.addConstr(quicksum(S[i] for i in self.assortment_sku_set) <= self.warehouse_capacity)

        model.addConstrs((quicksum(self.aggregate_demand_dict[i][region_id][s_id] * Pk[(s_id, i, region_id)]
                                     for region_id in self.instance.branch.region_set) <= S[i]
                           for s_id in self.instance.branch.scenario_set
                           for i in self.assortment_sku_set))

        # Define objective function
        first_stage_cost = quicksum(S[i] * self.instance.branch.sku_weight_dict[i] * 
                                     self.instance.branch.first_cost_structure[self.rdc_id][self.fdc_id]['coef'] 
                                     for i in self.assortment_sku_set)

        second_stage_cost = quicksum(self.instance.branch.scenario_prob_dict[s_id] * 
                                      quicksum(self.aggregate_demand_dict[i][region_id][s_id] * 
                                               (Pk[(s_id, i, region_id)] * 
                                                (self.instance.branch.second_cost_structure[self.fdc_id][region_id]['intercept'] * (1 / self.selected_involved_size_dict[i]) + 
                                                 self.instance.branch.sku_weight_dict[i] * self.instance.branch.second_cost_structure[self.fdc_id][region_id]['coef']) +
                                                (1 - Pk[(s_id, i, region_id)]) * 
                                                (self.instance.branch.second_cost_structure[self.rdc_id][region_id]['intercept'] * (1 / self.selected_involved_size_dict[i]) + 
                                                 self.instance.branch.sku_weight_dict[i] * self.instance.branch.second_cost_structure[self.rdc_id][region_id]['coef']))
                                               for i in self.assortment_sku_set
                                               for region_id in self.instance.branch.region_set)    
                                      for s_id in self.instance.branch.scenario_set)

        model.setObjective(first_stage_cost + second_stage_cost, GRB.MINIMIZE)
        model.optimize()
        # Obj_cost = model.ObjVal
        S_val = {(self.fdc_id, i): S[i].x for i in self.assortment_sku_set}

        self.basic_repooling_sol_on_sku = copy.deepcopy(self.input_sol_on_sku)
        self.basic_repooling_sol_on_sku['S'].update(S_val)
    

    def mean_sol(self, sol_a, sol_b):
        avg_sol = copy.deepcopy(sol_a)
        for key, val in sol_a['S'].items():
            avg_sol['S'][key] = np.ceil(0.5 * (val + sol_b['S'][key]))

        if self.warehouse_capacity is not None and sum(avg_sol['S'].values()) > self.warehouse_capacity:
            avg_sol = copy.deepcopy(sol_a)
            for key, val in sol_a['S'].items():
                avg_sol['S'][key] = round(0.5 * (val + sol_b['S'][key]))

        # avg_cost, _ = self.instance.solver.solve_Model_given_placement_plan(avg_sol)

        # avg_gap = (100 * (avg_cost - self.opt_cost) / (self.ub_cost - self.opt_cost)
        #            if self.opt_cost is not None else avg_cost)

        return avg_sol


    def __call__(self):
        self.demand_aggregation()
        self.generalized_LP()
    

    def get_pure_sol(self):
        return self.basic_repooling_sol_on_sku


    def get_mix_sol(self):
        return self.mean_sol(sol_a=self.input_sol_on_sku,
                             sol_b=self.basic_repooling_sol_on_sku)