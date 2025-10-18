from algorithms.solver import SOLVER
from orders.branch import OrderBranch
from default_paras import *
from utils import *
import numpy as np
import copy
import concurrent.futures
import os
import math

class ProgressiveHedging():
    def __init__(self, orders: OrderBranch, 
                 solver: SOLVER,
                 init_penalty_parm_x=INIT_PENALTY_PARM_X,
                 init_penalty_parm_s=INIT_PENALTY_PARM_S, 
                 penalty_multiplier=PENALTY_MULTIPLIER,
                 max_iteration_num=MAX_ITERATION_NUM, 
                 max_terminate_s=MAX_TERMINATE_S,
                 max_tolerance_x=MAX_TOLERANCE_X,
                 max_tolerance_s=MAX_TOLERANCE_S
                 ):
        self.orders = orders
        self.solver = solver
        self.rdc_id = orders.rdc_id
        self.fdc_id = orders.fdc_id
        self.sku_set = orders.sku_set
        self.scenario_set = solver.scenario_set
        self.scenario_prob = solver.scenario_prob_dict
        
        # 预计算并缓存概率数组
        self.scenario_probs = np.array([solver.scenario_prob_dict[s] for s in self.scenario_set])
        
        self.order_type_set = orders.order_type_set
        self.sku_type_map_dict = orders.sku_type_map_dict
        self.reverse_type_map_dict = {v: k for k, v in orders.type_map_dict.items()}
        self.type_size_dict = {k: len(v) for k, v in self.reverse_type_map_dict.items()}
        
        # 预计算并缓存type相关数据
        self.type_weights = np.array([orders.type_weight_dict[t] for t in self.order_type_set])
        self.type_sizes = np.array([self.type_size_dict[t] for t in self.order_type_set])
        
        self.type_sales_dict = orders.type_sales_dict
        self.type_weight_dict = orders.type_weight_dict
        self.UB_cost = self.orders.UB_cost
        self.complete_rdc_cost = self.orders.complete_rdc_cost
        self.demand_diff_pairs = self.orders.demand_diff_pairs
        self.worst_second_cost = self.orders.worst_second_cost
        
        self.init_penalty_parm_x = init_penalty_parm_x
        self.init_penalty_parm_s = init_penalty_parm_s 
        self.penalty_multiplier = penalty_multiplier
        self.max_iteration_num = max_iteration_num
        self.max_terminate_s = max_terminate_s
        self.max_tolerance_s = max_tolerance_s
        self.max_tolerance_x = max_tolerance_x
        
        self.penalty_flag = False
        self.penalty_parm_x_dict = None
        self.penalty_parm_s_dict = None
        self.final_sol = {}
        self.conv_record_set = []
        
        # 使用numpy数组优化初始化
        self.init_zero_structures()
        self.get_type_composition_by_sku()
    
    def init_zero_structures(self):
        """使用numpy数组初始化数据结构"""
        n_scenarios = len(self.scenario_set)
        n_types = len(self.order_type_set)
        
        # 初始化dual和solution结构
        self.zero_dual_x = np.zeros((n_scenarios, n_types))
        self.zero_dual_s = np.zeros((n_scenarios, n_types))
        self.zero_sol_x = np.zeros(n_types)
        self.zero_sol_s = np.zeros(n_types)
        
        # 转换为字典格式
        self.zero_dual = {s_id: {'X': {}, 'S': {}} for s_id in self.scenario_set}
        self.zero_sol = {'X': {}, 'S': {}}
        
        for i, s_id in enumerate(self.scenario_set):
            for j, t_id in enumerate(self.order_type_set):
                self.zero_dual[s_id]['X'][t_id] = 0
                self.zero_dual[s_id]['S'][t_id] = 0
                if i == 0:
                    self.zero_sol['X'][t_id] = 0
                    self.zero_sol['S'][t_id] = 0
        
        # 初始化reduction相关结构
        self.second_reduction_dict = {s_id: {type_id: 0 for type_id in self.order_type_set} 
                                    for s_id in self.scenario_set}
        self.completely_reduction_set = {type_id: 0 for type_id in self.order_type_set}


    def update_dual_value_per_iteration(self, sol, sol_by_scenario, last_dual):
        """使用numpy向量化计算更新dual values"""
        n_scenarios = len(self.scenario_set)
        n_types = len(self.unsolved_type_set)
        
        # 构建numpy数组
        penalty_x = np.array([self.penalty_parm_x_dict[t] for t in self.unsolved_type_set])
        penalty_s = np.array([self.penalty_parm_s_dict[t] for t in self.unsolved_type_set])
        
        sol_x = np.zeros((n_scenarios, n_types))
        sol_s = np.zeros((n_scenarios, n_types))
        target_x = np.array([sol['X'][t] for t in self.unsolved_type_set])
        target_s = np.array([sol['S'][t] for t in self.unsolved_type_set])
        
        # 填充数组
        for i, s_id in enumerate(self.scenario_set):
            for j, t_id in enumerate(self.unsolved_type_set):
                sol_x[i,j] = sol_by_scenario[s_id]['X'][t_id]
                sol_s[i,j] = sol_by_scenario[s_id]['S'][t_id]
        
        # 向量化计算
        add_dual_x = penalty_x[None,:] * (sol_x - target_x[None,:])
        add_dual_s = penalty_s[None,:] * (sol_s - target_s[None,:])
        
        # 转回字典格式
        dual = {s_id: {'X': {}, 'S': {}} for s_id in self.scenario_set}
        for i, s_id in enumerate(self.scenario_set):
            for j, t_id in enumerate(self.unsolved_type_set):
                if len(last_dual) > 0:
                    dual[s_id]['X'][t_id] = last_dual[s_id]['X'][t_id] + add_dual_x[i,j]
                    dual[s_id]['S'][t_id] = last_dual[s_id]['S'][t_id] + add_dual_s[i,j]
                else:
                    dual[s_id]['X'][t_id] = add_dual_x[i,j]
                    dual[s_id]['S'][t_id] = add_dual_s[i,j]
        
        return dual
    

    def eva_second_cost_by_type(self, type_id, scenario_id, sol):
        """优化第二阶段成本计算"""
        sol_s = sol[type_id]['S']
        worst_second_cost = self.worst_second_cost[type_id][scenario_id]
        
        if sol_s == 0:
            return worst_second_cost
            
        demand_diff_pairs = self.demand_diff_pairs[type_id][scenario_id]
        corr_demand_set = np.array(demand_diff_pairs['demand'])
        corr_diff_set = np.array(demand_diff_pairs['diff'])
        
        # 使用numpy的cumsum优化累加计算
        cum_demands = np.cumsum(corr_demand_set)
        idx = np.searchsorted(cum_demands, sol_s)
        
        if idx == len(corr_demand_set):
            reduction_cost = np.sum(corr_demand_set * corr_diff_set)
        else:
            if idx == 0:
                remain = sol_s
            else:
                remain = sol_s - cum_demands[idx-1]
            reduction_cost = np.sum(corr_demand_set[:idx] * corr_diff_set[:idx]) + \
                           remain * corr_diff_set[idx]
        
        self.second_reduction_dict[scenario_id][type_id] = reduction_cost
        return worst_second_cost - reduction_cost

    

    def get_type_composition_by_sku(self):
        # summary the related SKUs for each order type
        type_composition_by_sku = {type_id: [] for type_id in self.order_type_set}
        for k, v in self.orders.type_map_dict.items():
            corr_sku_set = [t[0] for t in k]
            type_composition_by_sku[v] = corr_sku_set
        self.type_composition_by_sku = type_composition_by_sku
    
    
    def eva_reduction_cost_given_sol(self, input_fdc, sol=None):
        rdc_id = self.orders.rdc_id
        first_coef = self.orders.first_cost_structure[rdc_id][input_fdc]['coef']

        sol = sol if sol is not None else self.solved_type_sol
        for type_id in self.order_type_set:
            for s_id in self.scenario_set:
                self.eva_second_cost_by_type(type_id, scenario_id=s_id, sol=sol)

            reduction_cost = sum([self.scenario_prob[s_id] * self.second_reduction_dict[s_id][type_id]
                                for s_id in self.scenario_set]) - (first_coef * sol[type_id]['S'] * self.type_weight_dict[type_id])
            self.completely_reduction_set[type_id] = reduction_cost
    

    def eva_object_cost_by_type(self, input_fdc, sol):
        rdc_id = self.orders.rdc_id
        first_coef = self.orders.first_cost_structure[rdc_id][input_fdc]['coef']
        first_intercept = self.orders.first_cost_structure[rdc_id][input_fdc]['intercept']

        sol_f = 1 if any(sol[type_id]['X'] > 0 for type_id in self.order_type_set) else 0

        first_cost = sol_f * first_intercept + sum(sol[type_id]['S'] * self.orders.type_weight_dict[type_id]
                                                    for type_id in self.order_type_set) * first_coef

        second_cost = sum(self.scenario_prob[s_id] * self.eva_second_cost_by_type(type_id, s_id, sol)
                        for s_id in self.scenario_set
                        for type_id in self.order_type_set)

        return first_cost + second_cost

        
    def get_optimal_sol_of_lagrange_per_unit(self, input_fdc, type_id, scenario_id, sol, dual, max_demand):
        rdc_id = self.orders.rdc_id
        first_coef = self.orders.first_cost_structure[rdc_id][input_fdc]['coef']
        weight = self.orders.type_weight_dict[type_id]
        type_len = self.type_size_dict[type_id]

        penalty_parm_s = self.penalty_parm_s_dict[type_id]
        penalty_parm_x = self.penalty_parm_x_dict[type_id]

        corr_demand_set = self.demand_diff_pairs[type_id][scenario_id]['demand']
        corr_diff_set = self.demand_diff_pairs[type_id][scenario_id]['diff']

        remain_demand = max_demand - sum(corr_demand_set)
        corr_diff_set.append(0)
        corr_demand_set.append(remain_demand)

        worst_second_cost = self.worst_second_cost[type_id][scenario_id]

        dual_s = dual[scenario_id]['S'][type_id]
        dual_x = dual[scenario_id]['X'][type_id]
        sol_s = sol['S'][type_id]
        sol_x = sol['X'][type_id]

        def _func_derivate_of_lagrange(t, s):
            return first_coef * weight + dual_s + self.penalty_flag * penalty_parm_s * (s - sol_s) - corr_diff_set[t+1] + self.mu_flag * self.cur_mu * type_len

        # use binary search to find the optimal solution
        def _adaptive_binary_search():
            lb = 0
            ub = len(corr_demand_set) - 1
            s_bs = 0

            while ub > lb:
                tb = int(np.floor((lb + ub) / 2))
                sb = sum(corr_demand_set[:tb+1])
                lhdv = _func_derivate_of_lagrange(tb-1, sb)
                rhdv = _func_derivate_of_lagrange(tb, sb)

                if (lhdv <= 0) and (rhdv >= 0):
                    return sb, tb-1
                else:
                    if (lhdv <= 0) and (rhdv <= 0):
                        temp = _func_derivate_of_lagrange(tb, sb + corr_demand_set[tb+1]) if self.penalty_flag else 0

                        if temp > 0:
                            s_bs = sol_s + ((corr_diff_set[tb+1] - first_coef * weight - dual_s -
                                    self.mu_flag * self.cur_mu * type_len) / penalty_parm_s)
                            return s_bs, tb
                        else:
                            lb = tb + 1
                    else:
                        ub = tb

        if _func_derivate_of_lagrange(0, 0) >= 0:
            subopt_s = 0
            subopt_t = 0
        elif _func_derivate_of_lagrange(len(corr_demand_set)-2, sum(corr_demand_set)) <= 0:
            subopt_s = sum(corr_demand_set)
            subopt_t = len(corr_demand_set) - 2
        else:
            subopt_s, subopt_t = _adaptive_binary_search()

        subopt_t = int(subopt_t)
        subopt_s = int(subopt_s)

        if not self.penalty_flag:
            return subopt_s

        if subopt_s == 0:
            opt_s = 0
        else:
            value_of_lagrange_on_zero = (-1) * dual_x * sol_x + (penalty_parm_x / 2) * (sol_x)**2 - dual_s * sol_s + (penalty_parm_s / 2) * (sol_s)**2 + worst_second_cost
            value_of_lagrange_on_subopt = first_coef * subopt_s * weight + dual_x * (1 - sol_x) + (penalty_parm_x / 2) * (1 - sol_x)**2 + dual_s * (subopt_s - sol_s) + (penalty_parm_s / 2) * (subopt_s - sol_s)**2 + worst_second_cost - (subopt_s * corr_diff_set[subopt_t+1] + sum([corr_demand_set[j] * (corr_diff_set[j] - corr_diff_set[subopt_t+1]) for j in range(subopt_t+1)])) + self.mu_flag * self.cur_mu * type_len * subopt_s

            if value_of_lagrange_on_zero <= value_of_lagrange_on_subopt:
                opt_s = 0
            else:
                opt_s = subopt_s

        return opt_s


    # def process_type_group(self, type_group, input_fdc, sol, dual):
    #     group_results = []
    #     for type_id in type_group:
    #         max_demand = max(self.type_sales_dict[type_id][s_id]['total'] for s_id in self.scenario_set)
    #         s_updates_by_scenario = {}
    #         x_updates_by_scenario = {}
    #         for s_id in self.scenario_set:
    #             opt_s = self.get_optimal_sol_of_lagrange_per_unit(input_fdc, type_id, s_id, sol, dual, max_demand)
    #             opt_x = 1 if opt_s > 0 else 0
    #             s_updates_by_scenario[s_id] = opt_s
    #             x_updates_by_scenario[s_id] = opt_x

    #         agg_s = round(sum(s_updates_by_scenario[s_id] * self.scenario_prob[s_id] for s_id in self.scenario_set))
    #         agg_x = round(sum(x_updates_by_scenario[s_id] * self.scenario_prob[s_id] for s_id in self.scenario_set))
    #         result = type_id, s_updates_by_scenario, x_updates_by_scenario, agg_s, agg_x
    #         group_results.append(result)
    #     return group_results


    # def update_expected_sol_per_iteration(self, input_fdc, sol, dual):
    #     import concurrent.futures
    #     import os
    #     max_workers = os.cpu_count() // 2
    #     type_list = list(self.unsolved_type_set)
    #     chunk_size = len(type_list) // max_workers
    #     type_groups = [type_list[i:i + chunk_size] for i in range(0, len(type_list), chunk_size)]
        
    #     update_sol_by_scenario = {s_id: {'S': {}, 'X': {}} for s_id in self.scenario_set}
    #     update_sol = {'S': {}, 'X': {}}
        
    #     with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
    #         futures = {executor.submit(self.process_type_group, group, input_fdc, sol, dual): i 
    #                 for i, group in enumerate(type_groups)}
            
    #         for future in concurrent.futures.as_completed(futures):
    #             group_results = future.result()
    #             for type_id, s_updates, x_updates, agg_s, agg_x in group_results:
    #                 update_sol['S'][type_id] = agg_s
    #                 update_sol['X'][type_id] = agg_x
    #                 for s_id in self.scenario_set:
    #                     update_sol_by_scenario[s_id]['S'][type_id] = s_updates[s_id]
    #                     update_sol_by_scenario[s_id]['X'][type_id] = x_updates[s_id]

    #     update_dual = self.update_dual_value_per_iteration(
    #         sol=update_sol, sol_by_scenario=update_sol_by_scenario, last_dual=dual)
    #     return update_sol_by_scenario, update_sol, update_dual


    def update_expected_sol_per_iteration(self, input_fdc, sol, dual):
        # update the expected solution after each iteration
        update_sol_by_scenario = {s_id: {'S': {}, 'X': {}} for s_id in self.scenario_set}
        update_sol = {'S': {}, 'X': {}}
        for type_id in self.unsolved_type_set:
            max_demand = max(self.type_sales_dict[type_id][s_id]['total'] for s_id in self.scenario_set)
            for s_id in self.scenario_set:
                opt_s = self.get_optimal_sol_of_lagrange_per_unit(input_fdc, type_id, s_id, sol, dual, max_demand)
                opt_x = 1 if opt_s > 0 else 0
                update_sol_by_scenario[s_id]['S'][type_id] = opt_s
                update_sol_by_scenario[s_id]['X'][type_id] = opt_x
            update_sol['S'][type_id] = round(sum(update_sol_by_scenario[s_id]['S'][type_id] * self.scenario_prob[s_id] for s_id in self.scenario_set))
            # update_sol['S'][type_id] = sum(update_sol_by_scenario[s_id]['S'][type_id] * self.scenario_prob[s_id] for s_id in self.scenario_set)
            update_sol['X'][type_id] = round(sum(update_sol_by_scenario[s_id]['X'][type_id] * self.scenario_prob[s_id] for s_id in self.scenario_set))
        update_dual = self.update_dual_value_per_iteration(sol=update_sol, sol_by_scenario=update_sol_by_scenario, last_dual=dual)
        return update_sol_by_scenario, update_sol, update_dual


    def update_unsolved_order_type(self, sol_by_scenario, sol, last_sol):
        # check terminate criterion and update solved order type
        scenario_set = list(sol_by_scenario.keys())
        solved_type_sol = {}
        unsolved_type_set = copy.deepcopy(self.unsolved_type_set)

        for type_id in unsolved_type_set:
            if sol['S'][type_id] <= self.max_terminate_s:
                self.unsolved_type_set.remove(type_id)
                sol_s = int(round(sol['S'][type_id]))
                sol_x = 0 if sol_s == 0 else 1
                solved_type_sol[type_id] = {'S': sol_s, 'X': sol_x}
                continue
            tol_s = np.sqrt(sum(abs((sol_by_scenario[s_id]['S'][type_id] - last_sol['S'][type_id])) * self.scenario_prob[s_id] for s_id in scenario_set))
            tol_x = np.sqrt(sum(abs((sol_by_scenario[s_id]['X'][type_id] - last_sol['X'][type_id])) * self.scenario_prob[s_id] for s_id in scenario_set))
            if (tol_s <= self.max_tolerance_s) and (tol_x <= self.max_tolerance_x):
                self.unsolved_type_set.remove(type_id)
                sol_s = int(round(sol['S'][type_id]))
                sol_x = 0 if sol_s == 0 else 1
                solved_type_sol[type_id] = {'S': sol_s, 'X': sol_x}

        return solved_type_sol

            
    def __call__(self, mu=0, mu_flag=False):
        input_fdc = self.fdc_id
        self.iteration_counter = 0
        self.solved_type_sol = {}
        self.sol_by_step_on_type = []
        self.sol_by_step_on_sku = []
        self.unsolved_type_set = list(self.order_type_set)  # Convert to list for easier manipulation

        self.penalty_parm_x_dict = {type_id: self.init_penalty_parm_x for type_id in self.order_type_set}
        self.penalty_parm_s_dict = {type_id: self.init_penalty_parm_s for type_id in self.order_type_set}

        # if need, incorporate mu to iteration
        self.mu_flag = mu_flag
        self.cur_mu = mu
        # initialization
        self.penalty_flag = False

        init_sol_by_scenario, init_sol, init_dual = self.update_expected_sol_per_iteration(input_fdc,
                                                                                          sol=self.zero_sol,
                                                                                          dual=self.zero_dual)
        self.penalty_flag = True
        cur_solved_type_sol = self.update_unsolved_order_type(init_sol_by_scenario, init_sol, init_sol)
        self.solved_type_sol.update(cur_solved_type_sol)
        cur_sol = init_sol
        cur_dual = init_dual
        cur_sol_by_scenario = init_sol_by_scenario
        last_sol = copy.deepcopy(cur_sol)
        self.pha_initial_sol = init_sol

        while self.iteration_counter <= self.max_iteration_num and self.unsolved_type_set:
            self.iteration_counter += 1
            last_sol = copy.deepcopy(cur_sol)
            # update the penalty parameters
            for type_id in self.unsolved_type_set:
                self.penalty_parm_x_dict[type_id] /= self.penalty_multiplier
                self.penalty_parm_s_dict[type_id] /= self.penalty_multiplier

            cur_sol_by_scenario, cur_sol, cur_dual = self.update_expected_sol_per_iteration(input_fdc, sol=cur_sol, dual=cur_dual)
            cur_solved_type_sol = self.update_unsolved_order_type(cur_sol_by_scenario, cur_sol, last_sol)
            self.solved_type_sol.update(cur_solved_type_sol)

        
        if self.unsolved_type_set:
            for type_id in self.unsolved_type_set:
                self.solved_type_sol[type_id] = {'S': cur_sol['S'][type_id], 'X': 0 if cur_sol['S'][type_id] == 0 else 1}

        self.algorithm_sol_on_sku = self.aggregate_sol_on_sku(input_fdc, sol=self.solved_type_sol)


    def aggregate_sol_on_sku(self, input_fdc, sol, type_flag=True):
        # aggregate the solution from order type level to SKU level
        algorithm_sol_on_sku = {'S': {(input_fdc, sku): 0 for sku in self.sku_set},
                      'X': {(input_fdc, sku): 0 for sku in self.sku_set}}

        for k, v in self.orders.type_map_dict.items():
            if v in self.order_type_set:
                S = sol[v]['S'] if type_flag else sol['S'][v]
                corr_sku_set = [t[0] for t in k]

                for sku in corr_sku_set:
                    algorithm_sol_on_sku['S'][(input_fdc, sku)] += S

        algorithm_sol_on_sku['X'] = {g: 1 if algorithm_sol_on_sku['S'][g] > 0 else 0 for g in algorithm_sol_on_sku['S']}
        final_sol = {'S': algorithm_sol_on_sku['S'],
                     'X': algorithm_sol_on_sku['X'],
                     'F': {input_fdc: 1 if sum(algorithm_sol_on_sku['X'].values()) > 0 else 0}}
        return final_sol


    def get_pha_paras(self):
        paras = {'init_penalty_parm_x': self.init_penalty_parm_x,
                 'init_penalty_parm_s': self.init_penalty_parm_s,
                 'penalty_multiplier': self.penalty_multiplier,
                 'max_iteration': self.max_iteration_num,
                 'max_terminate_s': self.max_terminate_s,
                 'max_tolerance_x': self.max_tolerance_x,
                 'max_tolerance_s': self.max_tolerance_s}
        return paras