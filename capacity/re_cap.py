from orders.branch import OrderBranch
from algorithms.solver import SOLVER
from algorithms.pha import ProgressiveHedging
from utils import *
from default_paras import *
from gurobipy import *
import copy

class Warehouse_CapConstr():
    def __init__(self, orders: OrderBranch, 
                 solver: SOLVER,
                 pha: ProgressiveHedging,
                 scaler_parm=SCALER_PARM,
                 fix_step_const=FIX_STEP_CONSTANT,
                 terminate_prop_parm=TERMINATE_PROP_PARM,
                 terminate_gap_parm=TERMINATE_GAP_PARM,
                 terminate_replicate_parm=TERMINATE_REPLICATE_PARM,
                 max_update_num=MAX_UPDATE_NUM):
        self.orders = orders
        self.solver = solver
        self.pha = pha
        self.fdc_id = orders.fdc_id
        self.scenario_set = solver.scenario_set
        self.order_type_set = orders.order_type_set
        self.sku_set = orders.sku_set
        self.fdc_capacity = solver.fdc_capacity
        self.type_size_dict = self.pha.type_size_dict
        self.init_scaler_parm = scaler_parm
        self.cur_scaler_parm = scaler_parm
        self.fix_step_const = fix_step_const
        self.terminate_flag = False
        self.terminate_prop_parm = terminate_prop_parm
        self.terminate_gap_parm = terminate_gap_parm
        self.terminate_replicate_parm = terminate_replicate_parm
        self.max_update_num = max_update_num
        self.mu_lb_info = (BigM, 0)
        self.mu_ub_info = ((-1) * self.fdc_capacity, BigM)
        self.init_mu = 0
        self.cur_mu = 0
        self.cur_step_size = 0
        self.final_sol_on_type = None
        self.final_sol_on_sku = None
        self.LB = 0
        self.UB = 0
        self.mu_pool = [0]
        self.delta_pool = []
        self.val_pool = []

    
    def solve_UB_zero(self):
        zero_sol_on_type = {type_id: {'S': 0, 'X': 0} for type_id in self.order_type_set}
        cost = self.pha.complete_rdc_cost
        self.UB = cost
        self.final_sol_on_type = zero_sol_on_type
    

    def compute_spill_capacity(self, sol, type_flag=True):
        if type_flag is True:
            used_capacity = sum(sol[type_id]['S'] * self.type_size_dict[type_id]
                                for type_id in self.order_type_set)
        else:
            used_capacity = sum(sol['S'].values())
        return used_capacity - self.fdc_capacity
    

    def eva_value_given_sol(self, input_fdc, sol_by_type):
        # evalute the objective cost by type
        cost = self.pha.eva_object_cost_by_type(input_fdc, sol_by_type)
        return cost


    def update_multipliers(self, input_fdc, sol_by_type,
                        # repair_sol_on_sku=None
                        ):
        delta = self.compute_spill_capacity(sol_by_type)
        if (delta > 0) & (delta <= self.mu_lb_info[0]):
            self.mu_lb_info = (delta, self.cur_mu)
        elif (delta < 0) & (delta >= self.mu_ub_info[0]):
            self.mu_ub_info = (delta, self.cur_mu)
        aux_cost = self.eva_value_given_sol(input_fdc, sol_by_type)
        cur_cost = aux_cost + self.cur_mu * delta
        self.val_pool.append(cur_cost)
        self.delta_pool.append(delta)
        gamma = self.cur_step_size * (self.UB - cur_cost) / (delta**2) \
                if delta != 0 else 0
        cur_mu = max(0, self.cur_mu + (delta * gamma))
        if (cur_mu >= self.mu_ub_info[1]) | (cur_mu <= self.mu_lb_info[1]):
            cur_mu = (self.mu_ub_info[1] + self.mu_lb_info[1]) / 2
            self.cur_scaler_parm /= 2
        self.cur_mu= cur_mu
        self.mu_pool.append(self.cur_mu)
        self.LB = max(self.LB, cur_cost)
        if delta <= 0:
            self.UB = min(aux_cost, self.UB)
            self.final_sol_on_type = copy.deepcopy(sol_by_type) if aux_cost == self.UB \
                                else self.final_sol_on_type
            if delta == 0:
                self.terminate_flag = True
    

    def update_origin_multipliers(self, input_fdc, sol_by_type):    
        delta = self.compute_spill_capacity(sol_by_type)
        aux_cost = self.eva_value_given_sol(input_fdc, sol_by_type)
        cur_cost = aux_cost + self.cur_mu * delta
        self.val_pool.append(cur_cost)
        self.delta_pool.append(delta)
        gamma = self.cur_step_size * (self.UB - cur_cost) / (delta**2) \
                if delta != 0 else 0
        cur_mu = max(0, self.cur_mu + (delta * gamma))
        self.cur_mu= cur_mu
        self.mu_pool.append(self.cur_mu)
        self.LB = max(self.LB, cur_cost)
        if delta <= 0:
            self.UB = min(aux_cost, self.UB)
            self.final_sol_on_type = copy.deepcopy(sol_by_type) if aux_cost == self.UB \
                                else self.final_sol_on_type
            if delta == 0:
                self.terminate_flag = True

    
    def check_terminate(self):
        # check terminate criterion
        return ((self.UB / self.LB) < (1 + self.terminate_prop_parm)) or \
               ((self.mu_ub_info[1] - self.mu_lb_info[1]) < self.terminate_gap_parm) or \
               ((len(self.delta_pool) - len(set(self.delta_pool))) > self.terminate_replicate_parm)
    

    def update_type_info(self, bar_S):
        pos_type_set_by_scenario = {scen_id: [] for scen_id in self.orders.scenario_set}
        neg_type_set_by_scenario = {scen_id: [] for scen_id in self.orders.scenario_set}
        for type_id in self.orders.order_type_set:
            if type_id not in self.orders.demand_lb_dict.keys():
                    continue
            for scen_id in self.orders.scenario_set:
                if self.orders.demand_lb_dict[type_id][scen_id] >= bar_S[type_id]:
                    pos_type_set_by_scenario[scen_id].append(type_id)
                else:
                    neg_type_set_by_scenario[scen_id].append(type_id)
        return pos_type_set_by_scenario, neg_type_set_by_scenario
        
    
    def __call__(self):
        input_fdc = self.fdc_id
        self.update_counter = 0
        self.solve_UB_zero()
        while (self.update_counter < self.max_update_num) & (self.terminate_flag == False):
            self.update_counter += 1
            self.cur_step_size = self.cur_scaler_parm * (1 + self.fix_step_const) / (self.update_counter + self.fix_step_const)
            self.pha(mu=self.cur_mu, mu_flag=True)
            sol_by_type = copy.deepcopy(self.pha.solved_type_sol)
            self.update_multipliers(input_fdc, sol_by_type)
            self.terminate_flag = self.check_terminate()
        self.final_sol_on_sku = self.pha.aggregate_sol_on_sku(input_fdc, self.final_sol_on_type)

    
    def get_algorithm_paras(self):
        paras_info = {'max_update_num': self.max_update_num,
                 'fix_step_const': self.fix_step_const,
                 'scaler_parm': self.init_scaler_parm,
                 'terminate_prop_parm': self.terminate_prop_parm,
                 'terminate_gap_parm': self.terminate_gap_parm,
                 'terminate_replicate_parm': self.terminate_replicate_parm
                 }
        return paras_info