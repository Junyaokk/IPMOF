import datetime
import os
import pandas as pd
from typing import Optional
from domain.instance import InstanceBase, Instance
from algorithms.repooling import Repooling_Opt

_pha_paras = ['init_penalty_parm_x',
                 'init_penalty_parm_s',
                 'penalty_multiplier',
                 'max_iteration',
                 'max_terminate_s',
                 'max_tolerance_x',
                 'max_tolerance_s']
class Task():
    def __init__(self, task_id,
                 instance_base: InstanceBase, 
                 capacity_constr=False, 
                 capacity_type=None,
                 need_solver=True,
                 need_opt=False,
                 need_benchmark=False,
                 input_fdc_capacity=None,
                 mixed_paras=None,
                 # todo: add
                 scenario_num=None,
                 algorithm_paras: Optional[dict] = None):
        self.task_id = task_id
        self.instance_base = instance_base
        self.capacity_constr = capacity_constr
        self.capacity_type = capacity_type
        self.input_fdc_capacity = input_fdc_capacity
        self.need_solver = need_solver
        self.need_opt = need_opt
        self.need_benchmark = need_benchmark
        self.scenario_num = scenario_num
        self.mixed_paras=mixed_paras
        self.pha_paras = None
        if algorithm_paras is None:
            self.algorithm_paras = {}
            self.pha_paras = {}
        else:
            self.pha_paras = {k: v for k, v in algorithm_paras.items() 
                              if k in _pha_paras}
            self.algorithm_paras = {k: algorithm_paras[k] 
                                   for k in algorithm_paras if k not in self.pha_paras}
        
        self.instance = Instance(self.instance_base)
        if (self.capacity_constr) & (self.capacity_type != 'mixed'):
            if input_fdc_capacity <= 1:
                if self.capacity_type == 'warehouse':
                    raise ValueError('Warehouse capacity does not support value lower than 1')
                else:
                    self.input_fdc_capacity = int(input_fdc_capacity * len(self.instance.source.sku_set))
            else:
                self.input_fdc_capacity = input_fdc_capacity

        self.task_info = None
        self.task_time = None
        self.benchmark_time = None
        self.opt_time = None
        self.algorithm_sol_on_sku = None
        self.opt_sol = None
        self.benchmark_sol = None
        self.algorithm_cost = None
        self.benchmark_cost = None
        self.opt_cost = None
        self.ub_cost = None
        self.opt_gap = None
        self.benchmark_gap = None
        self.ub_opt_gap = None
        
        self.instance.generation(scenario_num=scenario_num)
    

    def reset_algorithm_paras(self, algorithm_paras):
        self.pha_paras = {k: v for k, v in algorithm_paras.items() 
                            if k in _pha_paras}
        self.algorithm_paras = {k: algorithm_paras[k] 
                                for k in algorithm_paras if k not in self.pha_paras}

    
    def run(self,
            data_dir=None):
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        self.task_dir = data_dir + self.task_id + '/'
        self.instance.preprocess(self.capacity_constr, self.pha_paras)


        time0 = datetime.datetime.now()
        if self.instance.capmanager is None:
            self.instance.pha()
            self.algorithm_sol_on_sku = self.instance.pha.algorithm_sol_on_sku
            self.algorithm_sol_on_type = self.instance.pha.solved_type_sol
        else:
            self.instance.capmanager.invoke_CapConstr(self.input_fdc_capacity, 
                                                      cap_type=self.capacity_type, 
                                                      algorithm_paras=self.algorithm_paras,
                                                      mixed_paras=self.mixed_paras)
            self.algorithm_sol_on_sku = self.instance.capmanager.opt_sol_on_sku
            self.algorithm_sol_on_type = self.instance.capmanager.opt_sol_on_type
        
        self.task_time = (datetime.datetime.now() - time0).total_seconds()
        self.ub_cost = self.instance.pha.complete_rdc_cost
        
        if self.need_benchmark:
            self.get_benchmark_sol()
        if self.need_solver:
            self.compute_algorithm_cost()
        if self.need_opt:
            self.compute_opt_cost()

        self.gen_task_info()
        self.write_to_csv(data_dir=self.task_dir)

    
    def run_repooling(self, repooling_methods_type='F'):
        if repooling_methods_type == 'F':
            repooling_optimizer = Repooling_Opt(self.instance, 
                                                self.algorithm_sol_on_sku,
                                                self.algorithm_sol_on_type,
                                                'F',
                                                self.capacity_type)
        else:
            repooling_optimizer = Repooling_Opt(self.instance, 
                                                self.algorithm_sol_on_sku,
                                                self.algorithm_sol_on_type,
                                                'T',
                                                self.capacity_type)
        
        time0 = datetime.datetime.now()
        repooling_optimizer()
        repooling_time = (datetime.datetime.now() - time0).total_seconds()
        
        pure_sol = repooling_optimizer.get_pure_sol()
        mix_sol = repooling_optimizer.get_mix_sol()
        pure_cost, _ = self.instance.solver.solve_Model_given_placement_plan(pure_sol)
        mix_cost, _ = self.instance.solver.solve_Model_given_placement_plan(mix_sol)
        
        if self.opt_sol is not None:
            pure_gap = 100 * (pure_cost - self.opt_cost) / (self.ub_cost - self.opt_cost)
            mix_gap = 100 * (mix_cost - self.opt_cost) / (self.ub_cost - self.opt_cost)
        
        repooling_info = {
            'repooling_methods_type': repooling_methods_type,
            'repooling_time': repooling_time,
            'pure_cost': pure_cost,
            'mix_cost': mix_cost,
            'pure_gap (%)': pure_gap,
            'mix_gap (%)': mix_gap
        }
        
        repooling_info_df = pd.DataFrame.from_dict(repooling_info, orient='index').reset_index()
        repooling_info_df.columns = ['repooling_info', 'value']
        repooling_file_path = self.task_dir + str(self.task_id) + '_repooling_info.csv'
        header = not os.path.exists(repooling_file_path) or os.path.getsize(repooling_file_path) == 0
        repooling_info_df.to_csv(repooling_file_path, mode='a', header=header, index=False)
        

    def get_benchmark_sol(self):
        assortment_flag = True if self.capacity_type == 'assortment' else False
        time0 = datetime.datetime.now()
        self.benchmark_sol, _ = self.instance.solver.solve_benchmark_method(assortment_flag=assortment_flag)
        self.benchmark_time = (datetime.datetime.now() - time0).total_seconds()

    
    def compute_algorithm_cost(self):
        time0 = datetime.datetime.now()
        algorithm_cost, algorithm_second_stage_sol = self.instance.solver.solve_Model_given_placement_plan(self.algorithm_sol_on_sku)
        self.algorithm_validate_time = (datetime.datetime.now() - time0).total_seconds()
        self.algorithm_cost = algorithm_cost
        self.algorithm_second_stage_sol = algorithm_second_stage_sol
        
        if self.benchmark_sol is not None:
            benchmark_cost, benchmark_second_stage_sol = self.instance.solver.solve_Model_given_placement_plan(self.benchmark_sol)
            self.benchmark_second_stage_sol = benchmark_second_stage_sol
            self.benchmark_cost = benchmark_cost
            self.benchmark_algorithm_gap = ((self.benchmark_cost / self.algorithm_cost) - 1) * 100


    def compute_opt_cost(self):
        time0 = datetime.datetime.now()
        assortment_flag = True if self.capacity_type == 'assortment' else False
        self.opt_sol, self.opt_cost, self.opt_second_stage_sol = self.instance.solver.solve_Model_on_sku_completely(assortment_flag=assortment_flag)
        self.opt_time = (datetime.datetime.now() - time0).total_seconds()
        
        if self.algorithm_cost is not None:
            # todo: revise
            # self.opt_gap = ((self.algorithm_cost / self.opt_cost) - 1) * 100
            self.opt_gap = 100 * (self.algorithm_cost - self.opt_cost) / (self.ub_cost - self.opt_cost)
        
        if self.benchmark_cost is not None:
            # todo: revise
            # self.benchmark_gap = ((self.benchmark_cost / self.opt_cost) - 1) * 100
            self.benchmark_gap = 100 * (self.benchmark_cost - self.opt_cost) / (self.ub_cost - self.opt_cost)

    
    def gen_task_info(self):
        task_info = {
            'task_id': self.task_id,
            'instance_id': self.instance.retailer_id,
            'num.order': self.instance.branch.size,
            'scenario_num': self.instance.branch.scenario_num,
            'scenario_unit': self.instance.branch.scenario_unit,
            'capacity_constr': self.capacity_constr,
            'capacity_type': self.capacity_type,
            'capacity_info': self.input_fdc_capacity,
            'task_time': self.task_time,
            'opt_time': self.opt_time,
            'benchmark_time': self.benchmark_time,
            'algorithm_cost': self.algorithm_cost,
            'benchmark_cost': self.benchmark_cost,
            'ub_cost': self.ub_cost,
            'opt_cost': self.opt_cost,
            'opt_gap (algorithm) %': self.opt_gap,
            'opt_gap (benchmark) %': self.benchmark_gap
        }

        task_info.update(self.instance.pha.get_pha_paras())
        
        if self.capacity_constr:
            if self.capacity_type != 'assortment':
                task_info.update(self.instance.capmanager.CapConstr.get_algorithm_paras())
        self.task_info = task_info
    
    
    def write_to_csv(self, data_dir):
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        
        input_fdc = self.instance.source.fdc_id
        if self.opt_sol is None:
            sol_list = [(sku, self.algorithm_sol_on_sku['X'][(input_fdc, sku)], self.algorithm_sol_on_sku['S'][(input_fdc, sku)])
                    for sku in self.instance.source.sku_set]
            sol_df = pd.DataFrame(sol_list, columns=['sku_id', 'X', 'S'])
        else:
            sol_list = [(sku, self.algorithm_sol_on_sku['X'][(input_fdc, sku)], self.opt_sol['X'][(input_fdc, sku)],
                        self.algorithm_sol_on_sku['S'][(input_fdc, sku)], self.opt_sol['S'][(input_fdc, sku)])
                        for sku in self.instance.source.sku_set]
            sol_df = pd.DataFrame(sol_list, columns=['sku_id', 'X', 'opt_X', 'S', 'opt_S'])

        task_info_df = pd.DataFrame.from_dict(self.task_info, orient='index').reset_index()
        task_info_df.columns = ['task_info', 'value']

        sol_df.to_csv(data_dir + str(self.task_id) + '_sol.csv', index=False)
        task_info_df.to_csv(data_dir + str(self.task_id) + '_task_info.csv', index=False)