from utils import *
from default_paras import BigM, TIME_LIMIT
from gurobipy import *
from orders.branch import OrderBranch

class SOLVER():
    def __init__(self, orders: OrderBranch):
        self.orders = orders
        self.scenario_set = orders.scenario_set
        self.scenario_sales_by_region = orders.scenario_sales_by_region
        self.scenario_pool_by_region = orders.scenario_pool_by_region
        self.sku_type_map_dict = orders.sku_type_map_dict
        self.type_element_dict = orders.type_element_dict
        self.type_weight_dict = orders.type_weight_dict
        self.scenario_prob_dict = orders.scenario_prob_dict
        self.sku_demand_dict = orders.sku_demand_dict
        self.sku_sales_dict = orders.sku_sales_dict
        self.scenario_num = len(self.scenario_set)
        self.fdc_capacity = BigM
        self.time_limit = TIME_LIMIT
        self.mixed_flag = False
        self.mixed_paras = None
        

    def reset_capacity_info(self, input_capacity_val):
        self.fdc_capacity = input_capacity_val

    
    def reset_mixed_info(self, mixed_paras=None):
        self.mixed_flag = True
        self.mixed_paras = mixed_paras
        self.assortment_capacity = mixed_paras['assortment']
        self.warehouse_capacity = mixed_paras['warehouse']
        
        
    def solve_PO_model_on_specific_scenario(self, input_scenario_set,
                                        assortment_flag=False,
                                        timelimit_flag=False,
                                        input_timilimit=None):
        scenario_set = input_scenario_set
        rdc_id = self.orders.rdc_id
        fdc_id = self.orders.fdc_id
        sku_set = self.orders.sku_set
        region_set = self.orders.region_set
        sku_weight_dict = self.orders.sku_weight_dict
        type_weight_dict = self.type_weight_dict
        sku_type_map_dict = self.sku_type_map_dict
        type_element_dict = self.type_element_dict
        first_cost_structure = self.orders.first_cost_structure
        second_cost_structure = self.orders.second_cost_structure
        scenario_prob_dict = self.scenario_prob_dict
        scenario_pool_by_region = self.scenario_pool_by_region
        scenario_sales_by_region = self.scenario_sales_by_region
        bigm_per_sku = self.orders.bigm_per_sku

        
        model = Model('PO-Opt')
        model.params.OutputFlag = 0
        if timelimit_flag is True:
            model.params.OutputFlag = 1
            model.setParam('LogToConsole', 0)
            model.setParam('LogFile', './gurobi.log')
            model.Params.TimeLimit = self.time_limit
        else:
            if input_timilimit is not None:
                model.params.OutputFlag = 1
                model.setParam('LogToConsole', 0)
                model.setParam('LogFile', './gurobi.log')
                model.Params.TimeLimit = input_timilimit

        M = BigM
        fdc_capacity = self.fdc_capacity
        zk_index = []
        
        for region_id in region_set:
            for scen_id in scenario_set:
                for type_id in scenario_pool_by_region[region_id][scen_id]:
                    zk_index.append((region_id, scen_id, type_id))
        
        # replenishment decisions
        S = model.addVars(sku_set, vtype=GRB.CONTINUOUS, lb=0, name="supply")
        # assortment decisions
        X = model.addVars(sku_set, vtype=GRB.BINARY, name="x")
        # transshipment decisions
        F = model.addVar(vtype=GRB.BINARY, name="f")
        # fulfillment decisions
        Zk = model.addVars(zk_index, vtype=GRB.CONTINUOUS, lb=0, name="zk")

        # add first-stage constraints
        model.addConstrs((X[i] <= S[i] for i in sku_set))
        # model.addConstrs((S[i] <= M * X[i] for i in sku_set))
        # model.addConstr((quicksum(X[i] for i in sku_set) <= M * F))
        model.addConstrs((S[i] <= bigm_per_sku[i] * X[i] for i in sku_set))
        model.addConstr((quicksum(X[i] for i in sku_set) <= len(sku_set) * F))

        if self.mixed_flag:
            model.addConstr((quicksum(X[i] for i in sku_set) <= self.assortment_capacity))
            model.addConstr((quicksum(S[i] for i in sku_set) <= self.warehouse_capacity))
        else:
            if assortment_flag:
                model.addConstr((quicksum(X[i] for i in sku_set) <= fdc_capacity))
            else:
                model.addConstr((quicksum(S[i] for i in sku_set) <= fdc_capacity))

        # add second-stage constraints
        model.addConstrs((quicksum(scenario_sales_by_region[region_id][scen_id][type_id] \
                                   * type_element_dict[type_id][i] * Zk[region_id, scen_id, type_id]
                                   for region_id in region_set
                                   for type_id in (set(scenario_pool_by_region[region_id][scen_id]) & set(sku_type_map_dict[i])))
                            <= S[i]
                        for i in sku_set
                        for scen_id in scenario_set
                        ))
        model.update()
        model.addConstrs((Zk[region_id, scen_id, type_id] <= X[i]
                    for i in sku_set
                    for region_id in region_set
                    for scen_id in scenario_set
                    for type_id in (set(scenario_pool_by_region[region_id][scen_id]) & set(sku_type_map_dict[i]))
                    ))
        
        # add objective
        first_stage_cost = F * first_cost_structure[rdc_id][fdc_id]['intercept'] \
                            + quicksum([S[i] * sku_weight_dict[i] * first_cost_structure[rdc_id][fdc_id]['coef'] 
                            for i in sku_set])                          
        second_stage_cost = quicksum(scenario_prob_dict[scen_id] * 
                                quicksum(scenario_sales_by_region[region_id][scen_id][type_id] *
                                         (Zk[region_id, scen_id, type_id] * (second_cost_structure[fdc_id][region_id]['intercept'] \
                                                                                     + second_cost_structure[fdc_id][region_id]['coef'] * type_weight_dict[type_id])
                                            + (1 - Zk[region_id, scen_id, type_id]) * (second_cost_structure[rdc_id][region_id]['intercept'] \
                                                                                     + second_cost_structure[rdc_id][region_id]['coef'] * type_weight_dict[type_id]))
                                         for region_id in region_set for type_id in scenario_pool_by_region[region_id][scen_id])
                                        for scen_id in scenario_set)

        model.setObjective(first_stage_cost + second_stage_cost, GRB.MINIMIZE)
        model.optimize()

        # get optimized solution
        Obj_cost = model.ObjVal
        # todo: revise
        S_val = {(fdc_id, i): S[i].x for i in sku_set}
        X_val = {(fdc_id, i): X[i].x for i in sku_set}
        F_val = {fdc_id: F.x }
        sol = {'S': S_val, 'X': X_val, 'F': F_val}
        Zk_val = {zk_i: Zk[zk_i].x for zk_i in zk_index}
        
        return sol, Obj_cost, Zk_val
    

    def solve_PO_model(self, assortment_flag=False, timelimit_flag=False, input_timilimit=None):
        # todo: revise
        sol, _, _ = self.solve_PO_model_on_specific_scenario(self.scenario_set, 
                                                            assortment_flag=assortment_flag,
                                                            timelimit_flag=timelimit_flag,
                                                            input_timilimit=input_timilimit)
        Obj_cost, Zk_val = self.solve_PO_model_given_first_plan(placement_plan=sol)
        return sol, Obj_cost, Zk_val


    def solve_PO_model_given_first_plan(self, placement_plan, input_scenario_set=[]):
        scenario_set = input_scenario_set if len(input_scenario_set) else self.scenario_set
        rdc_id = self.orders.rdc_id
        fdc_id = self.orders.fdc_id
        sku_set = self.orders.sku_set
        region_set = self.orders.region_set
        sku_weight_dict = self.orders.sku_weight_dict
        type_weight_dict = self.type_weight_dict
        sku_type_map_dict = self.sku_type_map_dict
        type_element_dict = self.type_element_dict
        first_cost_structure = self.orders.first_cost_structure
        second_cost_structure = self.orders.second_cost_structure
        scenario_prob_dict = {scen_id: 1 / len(scenario_set) for scen_id in scenario_set}
        scenario_pool_by_region = self.scenario_pool_by_region
        scenario_sales_by_region = self.scenario_sales_by_region

        
        model = Model('PO-Opt-Validation')
        model.params.OutputFlag = 0
        fdc_capacity = self.fdc_capacity
        zk_index = []
        
        for region_id in region_set:
            for scen_id in scenario_set:
                for type_id in scenario_pool_by_region[region_id][scen_id]:
                    zk_index.append((region_id, scen_id, type_id))

        # replenishment solution
        bar_S = placement_plan['S']
        # assortment solution
        bar_X = placement_plan['X']
        # transshipment solution
        bar_F = placement_plan['F'][fdc_id]
        # fulfillment decisions
        Zk = model.addVars(zk_index, vtype=GRB.CONTINUOUS, lb=0, name="zk")


        # add second-stage constraints
        model.addConstrs((quicksum(scenario_sales_by_region[region_id][scen_id][type_id] \
                                   * type_element_dict[type_id][i] * Zk[region_id, scen_id, type_id]
                                   for region_id in region_set
                                   for type_id in (set(scenario_pool_by_region[region_id][scen_id]) & set(sku_type_map_dict[i])))
                            <= bar_S[(fdc_id, i)]
                        for i in sku_set
                        for scen_id in scenario_set
                        ))
        model.update()
        model.addConstrs((Zk[region_id, scen_id, type_id] <= bar_X[(fdc_id, i)]
                    for i in sku_set
                    for region_id in region_set
                    for scen_id in scenario_set
                    for type_id in (set(scenario_pool_by_region[region_id][scen_id]) & set(sku_type_map_dict[i]))
                    ))
        
        # add objective
        placement_cost = bar_F * first_cost_structure[rdc_id][fdc_id]['intercept'] \
                            + quicksum([bar_S[(fdc_id, i)] * sku_weight_dict[i] * first_cost_structure[rdc_id][fdc_id]['coef'] 
                            for i in sku_set])                          
        second_stage_cost = quicksum(scenario_prob_dict[scen_id] * 
                                quicksum(scenario_sales_by_region[region_id][scen_id][type_id] *
                                         (Zk[region_id, scen_id, type_id] * (second_cost_structure[fdc_id][region_id]['intercept'] \
                                                                                     + second_cost_structure[fdc_id][region_id]['coef'] * type_weight_dict[type_id])
                                            + (1 - Zk[region_id, scen_id, type_id]) * (second_cost_structure[rdc_id][region_id]['intercept'] \
                                                                                     + second_cost_structure[rdc_id][region_id]['coef'] * type_weight_dict[type_id]))
                                         for region_id in region_set for type_id in scenario_pool_by_region[region_id][scen_id])
                                        for scen_id in scenario_set)
        
        model.setObjective(second_stage_cost, GRB.MINIMIZE)
        model.optimize()

        # get optimized solution
        Obj_cost = model.ObjVal + placement_cost
        Obj_cost = Obj_cost.getValue()
        Zk_val = {zk_i: Zk[zk_i].x for zk_i in zk_index}

        return Obj_cost, Zk_val

    def solve_OO_model(self, input_scenario_set=[],
                       assortment_flag=False,
                       given_first_plan=None):
        scenario_set = input_scenario_set if len(input_scenario_set) else self.scenario_set
        rdc_id = self.orders.rdc_id
        fdc_id = self.orders.fdc_id
        sku_set = self.orders.sku_set
        region_set = self.orders.region_set
        type_weight_dict = self.type_weight_dict
        sku_type_map_dict = self.sku_type_map_dict
        type_element_dict = self.type_element_dict
        first_cost_structure = self.orders.first_cost_structure
        second_cost_structure = self.orders.second_cost_structure
        scenario_prob_dict = self.scenario_prob_dict
        scenario_pool_by_region = self.scenario_pool_by_region
        scenario_sales_by_region = self.scenario_sales_by_region
        bigm_per_type = self.orders.bigm_per_type
        type_set = self.orders.order_type_set
        type_size_dict = self.orders.type_size_dict
        
        model = Model('OO-Opt')
        model.params.OutputFlag = 0
        # todo: set time limit for oo model
        # model.Params.TimeLimit = self.time_limit
        M = BigM
        fdc_capacity = self.fdc_capacity
        zk_index = []
        for region_id in region_set:
            for scen_id in scenario_set:
                # for type_id in scenario_pool_by_region[region_id][scen_id]:
                for type_id in type_set:
                    zk_index.append((region_id, scen_id, type_id))

        if given_first_plan is None:
            Stype = model.addVars(type_set, vtype=GRB.CONTINUOUS, lb=0, name="supply")
            Xtype = model.addVars(type_set, vtype=GRB.BINARY, name="x")
            Xi = model.addVars(sku_set, vtype=GRB.BINARY, name="xi")
            Ftype = model.addVar(vtype=GRB.BINARY, name="f")
            Zk = model.addVars(zk_index, vtype=GRB.CONTINUOUS, lb=0, name="zk")

            # add first-stage constraints
            # model.addConstrs((Xtype[type_id] <= Stype[type_id] for type_id in type_set))
            model.addConstrs((Stype[type_id] <= bigm_per_type[type_id] * Xtype[type_id] for type_id in type_set))
            model.addConstr((quicksum(Xtype[type_id] for type_id in type_set) <= len(type_set) * Ftype))
            
            if self.mixed_flag:
                model.addConstrs((Xtype[type_id] <= Xi[i] for i in sku_set for type_id in sku_type_map_dict[i]))
                model.addConstr((quicksum(Xi[i] for i in sku_set) <= self.assortment_capacity))
                model.addConstr((quicksum(Stype[type_id] * type_size_dict[type_id] for type_id in type_set) <= self.warehouse_capacity))
            else:
                if assortment_flag:
                    model.addConstrs((Xtype[type_id] <= Xi[i] for i in sku_set for type_id in sku_type_map_dict[i]))
                    model.addConstr((quicksum(Xi[i] for i in sku_set) <= fdc_capacity))
                else:
                    model.addConstr((quicksum(Stype[type_id] * type_size_dict[type_id] for type_id in type_set) <= fdc_capacity))

            # add second-stage constraints
            model.addConstrs((quicksum(scenario_sales_by_region[region_id][scen_id].get(type_id, 0) * 
                                    Zk[region_id, scen_id, type_id]
                                    for region_id in region_set
                                    )
                                <= Stype[type_id]
                            for type_id in type_set
                            for scen_id in scenario_set
                            ))
            model.addConstrs((Zk[region_id, scen_id, type_id] <= Xtype[type_id]
                            for region_id in region_set
                            for scen_id in scenario_set
                            for type_id in type_set)
                            )
            model.update()

            # add objective
            first_stage_cost = Ftype * first_cost_structure[rdc_id][fdc_id]['intercept'] \
                                + quicksum([Stype[type_id] * type_weight_dict[type_id] * first_cost_structure[rdc_id][fdc_id]['coef'] 
                                for type_id in type_set])                          
            second_stage_cost = quicksum(scenario_prob_dict[scen_id] * 
                                    quicksum(scenario_sales_by_region[region_id][scen_id][type_id] *
                                            (Zk[region_id, scen_id, type_id] * (second_cost_structure[fdc_id][region_id]['intercept'] \
                                                                                        + second_cost_structure[fdc_id][region_id]['coef'] * type_weight_dict[type_id])
                                                + (1 - Zk[region_id, scen_id, type_id]) * (second_cost_structure[rdc_id][region_id]['intercept'] \
                                                                                        + second_cost_structure[rdc_id][region_id]['coef'] * type_weight_dict[type_id]))
                                            for region_id in region_set for type_id in scenario_pool_by_region[region_id][scen_id])
                                            for scen_id in scenario_set)
        else:
            bar_Stype = {}
            bar_Xtype = {}
            for type_id, plan in given_first_plan.items():
                bar_Stype[type_id] = plan['S']
                bar_Xtype[type_id] = plan['X']
            
            Xi = model.addVars(sku_set, vtype=GRB.BINARY, name="xi")
            Ftype = model.addVar(vtype=GRB.BINARY, name="f")
            Zk = model.addVars(zk_index, vtype=GRB.CONTINUOUS, lb=0, name="zk")
            model.addConstr((quicksum(bar_Xtype[type_id] for type_id in type_set) <= len(type_set) * Ftype))
            
            model.addConstrs((bar_Xtype[type_id] <= Xi[i] for i in sku_set for type_id in sku_type_map_dict[i]))
            # add second-stage constraints
            model.addConstrs((quicksum(scenario_sales_by_region[region_id][scen_id].get(type_id, 0) * 
                                    Zk[region_id, scen_id, type_id]
                                    for region_id in region_set
                                    )
                                <= bar_Stype[type_id]
                            for type_id in type_set
                            for scen_id in scenario_set
                            ))
            model.addConstrs((Zk[region_id, scen_id, type_id] <= bar_Xtype[type_id]
                            for region_id in region_set
                            for scen_id in scenario_set
                            for type_id in type_set)
                            )
            model.update()

            # add objective
            first_stage_cost = Ftype * first_cost_structure[rdc_id][fdc_id]['intercept'] \
                                + quicksum([bar_Stype[type_id] * type_weight_dict[type_id] * first_cost_structure[rdc_id][fdc_id]['coef'] 
                                for type_id in type_set])                          
            second_stage_cost = quicksum(scenario_prob_dict[scen_id] * 
                                    quicksum(scenario_sales_by_region[region_id][scen_id][type_id] *
                                            (Zk[region_id, scen_id, type_id] * (second_cost_structure[fdc_id][region_id]['intercept'] \
                                                                                        + second_cost_structure[fdc_id][region_id]['coef'] * type_weight_dict[type_id])
                                                + (1 - Zk[region_id, scen_id, type_id]) * (second_cost_structure[rdc_id][region_id]['intercept'] \
                                                                                        + second_cost_structure[rdc_id][region_id]['coef'] * type_weight_dict[type_id]))
                                            for region_id in region_set for type_id in scenario_pool_by_region[region_id][scen_id])
                                            for scen_id in scenario_set)

        model.setObjective(first_stage_cost + second_stage_cost, GRB.MINIMIZE)
        model.optimize()

        # get optimized solution
        Obj_cost = model.ObjVal

        sol = {}
        if given_first_plan is None:
            # todo: revise
            S_val = {(fdc_id, i): sum(Stype[type_id].x * type_element_dict[type_id][i]
                                    for type_id in sku_type_map_dict[i]) for i in sku_set}
            X_val = {(fdc_id, i): 1 if S_val[(fdc_id, i)] > 0 else 0 for i in sku_set}
            F_val = {fdc_id: Ftype.x}
            sol = {'S': S_val, 'X': X_val, 'F': F_val}
        
        Zk_val = {zk_i: Zk[zk_i].x for zk_i in zk_index}

        return sol, Obj_cost, Zk_val
                    
    
    def solve_PP_model(self, input_scenario_set=[], 
                               assortment_flag=False):
        scenario_set = input_scenario_set if len(input_scenario_set) else self.scenario_set
        rdc_id = self.orders.rdc_id
        fdc_id = self.orders.fdc_id
        region_set = self.orders.region_set
        sku_set = self.orders.sku_set
        sku_weight_dict = self.orders.sku_weight_dict
        first_cost_structure = self.orders.first_cost_structure
        second_cost_structure = self.orders.second_cost_structure
        scenario_prob_dict = self.scenario_prob_dict
        fdc_capacity = self.fdc_capacity
        scenario_sales_by_region = self.scenario_sales_by_region
        scenario_pool_by_region = self.scenario_pool_by_region
        type_element_dict = self.orders.type_element_dict
        sku_type_map_dict = self.orders.sku_type_map_dict
        avg_involved_size_dict = self.orders.avg_involved_size_dict
        demand_on_sku = {sku: {region_id: {scen_id: 0 for scen_id in scenario_set} for region_id in region_set} for sku in sku_set}
        bigm_per_sku = self.orders.bigm_per_sku if hasattr(self.orders, 'bigm_per_sku') else {i: BigM for i in sku_set}
        # aggregate the demand from orders
        # demand_on_sku = copy.deepcopy(self.sku_demand_dict)
        # for s_id in scenario_set:
        #     for sku in sku_set:
        #         if sku not in list(demand_on_sku[s_id].keys()):
        #             demand_on_sku[s_id][sku] = {region_id: 0 for region_id in region_set}
        #         else:
        #             for region_id in region_set:
        #                 if region_id not in list(demand_on_sku[s_id][sku].keys()):
        #                     demand_on_sku[s_id][sku][region_id] = 0
        #                 else:
        #                     order_set = demand_on_sku[s_id][sku][region_id]
        #                     demand_on_sku[s_id][sku][region_id] = sum(list(order_set.values()))

        for sku in sku_set:
            for region_id in region_set:
                for scen_id in scenario_set:
                    demand = sum([scenario_sales_by_region[region_id][scen_id][type_id] * type_element_dict[type_id][sku]
                                for type_id in (set(sku_type_map_dict[sku]) & set(scenario_pool_by_region[region_id][scen_id]))])
                    demand_on_sku[sku][region_id][scen_id] = demand
        
        model = Model('PP-Opt')
        M = BigM
        fdc_capacity = self.fdc_capacity
        model.params.OutputFlag = 0
        # model.Params.TimeLimit = self.time_limit
        # replenishment decisions
        S = model.addVars(sku_set, vtype=GRB.CONTINUOUS, lb=0, name="supply")
        # assortment decisions
        X = model.addVars(sku_set, vtype=GRB.BINARY, name="x")
        # transshipment decisions
        F = model.addVar(vtype=GRB.BINARY, name="f")
        # fulfillment decisions
        Pk = model.addVars(scenario_set, sku_set, region_set, lb=0, name="pk")

        # add first-stage constraints
        model.addConstrs((X[i] <= S[i] for i in sku_set))
        # model.addConstrs((S[i] <= M * X[i] for i in sku_set))
        model.addConstrs((S[i] <= bigm_per_sku[i] * X[i] for i in sku_set))
        model.addConstr((quicksum(X[i] for i in sku_set) >= F))
        model.addConstr((quicksum(X[i] for i in sku_set) <= len(sku_set) * F))

        if self.mixed_flag:
            model.addConstr((quicksum(X[i] for i in sku_set) <= self.assortment_capacity))
            model.addConstr((quicksum(S[i] for i in sku_set) <= self.warehouse_capacity))
        else:
            if assortment_flag:
                model.addConstr((quicksum(X[i] for i in sku_set) <= fdc_capacity))
            else:
                model.addConstr((quicksum(S[i] for i in sku_set) <= fdc_capacity))

        model.addConstrs((quicksum([demand_on_sku[i][region_id][s_id] * Pk[(s_id, i, region_id)]
                        for region_id in region_set])
                <= S[i]
                for s_id in scenario_set
                for i in sku_set
                ))
        model.addConstrs((Pk[(s_id, i, region_id)] <= X[i]
                    for s_id in scenario_set
                    for i in sku_set
                    for region_id in region_set
                    ))
        model.update()

        # add objective
        first_stage_cost = F * first_cost_structure[rdc_id][fdc_id]['intercept'] \
                    + quicksum([S[i] * sku_weight_dict[i] * first_cost_structure[rdc_id][fdc_id]['coef'] 
                    for i in sku_set])  
        
        second_stage_cost = quicksum([scenario_prob_dict[s_id] * 
                            quicksum([demand_on_sku[i][region_id][s_id] * (Pk[(s_id, i, region_id)] * (second_cost_structure[fdc_id][region_id]['intercept'] * (1 / avg_involved_size_dict[i]) + \
                                                                                                       sku_weight_dict[i] * second_cost_structure[fdc_id][region_id]['coef'])
                                                                    + (1 - Pk[(s_id, i, region_id)]) * (second_cost_structure[rdc_id][region_id]['intercept'] * (1 / avg_involved_size_dict[i]) + \
                                                                                                        sku_weight_dict[i] * second_cost_structure[rdc_id][region_id]['coef']))
                                    for i in list(set(self.orders.sku_set) & set(self.sku_sales_dict.keys()))  
                                    for region_id in region_set])    
                        for s_id in scenario_set])

        model.setObjective(first_stage_cost + second_stage_cost, GRB.MINIMIZE)
        model.optimize()
        # get optimized solution
        Obj_cost = model.ObjVal
        S_val = {(fdc_id, i): S[i].x for i in sku_set}
        X_val = {(fdc_id, i): X[i].x for i in sku_set}
        F_val = {fdc_id: F.x}
        sol = {'S': S_val, 'X': X_val, 'F': F_val}
        
        return sol, Obj_cost