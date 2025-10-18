import datetime
import os
import pandas as pd
from typing import Optional
from domain.instance import Instance
from algorithms.repooling import Repooling_Opt
from tqdm import tqdm

_pha_paras = ['init_penalty_parm_x',
                 'init_penalty_parm_s',
                 'penalty_multiplier',
                 'max_iteration',
                 'max_terminate_s',
                 'max_tolerance_x',
                 'max_tolerance_s']
class Task():
    def __init__(self, task_id,
                 instance: Instance, 
                 capacity_constr=False, 
                 capacity_type=None,
                 need_solver=True,
                 need_opt_po=False,
                 need_opt_pp=False,
                 need_opt_oo=False,
                 need_limited_po=False,
                 input_fdc_capacity=None,
                 mixed_paras=None,
                 scenario_num=None,
                 algorithm_paras: Optional[dict] = None):
        self.task_id = task_id
        self.instance = instance
        self.capacity_constr = capacity_constr
        self.capacity_type = capacity_type
        self.input_fdc_capacity = input_fdc_capacity
        self.need_solver = need_solver
        self.need_opt_po = need_opt_po
        self.need_opt_pp = need_opt_pp
        self.need_opt_oo = need_opt_oo
        self.need_limited_po = need_limited_po
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

        self.opt_po_time = None
        self.opt_po_sol = None
        self.opt_po_cost = None

        self.limited_po_time = None
        self.limited_po_sol = None
        self.limited_po_cost = None
        self.ub_opt_limited_po_gap = None
        self.opt_limite_po_gap = None

        self.algorithm_sol_on_sku = None
        self.algorithm_cost = None
        self.algo_gap = None
        self.ub_algo_gap = None

        self.opt_pp_time = None
        self.opt_pp_sol = None
        self.opt_pp_cost = None
        self.opt_pp_gap = None
        self.ub_opt_pp_gap = None

        self.ub_cost = None
        self.ub_opt_gap = None

        self.opt_oo_time = None
        self.opt_oo_cost = None
        self.opt_oo_gap = None
        self.ub_opt_oo_gap = None
        self.opt_oo_sol = None

        self.repooling_info = {}
        
        # self.instance.generation(scenario_num=scenario_num)
        
    

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

        print(f'Run task {self.task_id} ...')

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

        methods_to_run = []
        if self.need_opt_pp:
            methods_to_run.append(('opt_pp', self.get_opt_pp_sol))
        if self.need_opt_oo:
            methods_to_run.append(('opt_oo', self.get_opt_oo_sol))
        if self.need_solver:
            methods_to_run.append(('algorithm', self.compute_algorithm_cost))
        if self.need_limited_po:
            methods_to_run.append(('limited_po', lambda: self.compute_opt_cost(limited_flag=True)))
        if self.need_opt_po:
            methods_to_run.append(('opt_po', lambda: self.compute_opt_cost()))

        for method_name, method_func in tqdm(methods_to_run, 
                                           desc="Computing solutions", 
                                           unit="method"):
            print(f"\nProcessing {method_name}...")
            method_func()

        self.gen_task_info()
        self.write_to_file()

    
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
        pure_cost, _ = self.instance.solver.solve_PO_model_given_first_plan(pure_sol)
        pure_gap = None

        if self.opt_po_sol is not None:
            pure_gap = 100 * (pure_cost - self.opt_po_cost) / self.opt_po_cost
        
        repooling_info = {
            'repooling_methods_type': repooling_methods_type,
            'repooling_time': repooling_time,
            'pure_cost': pure_cost,
            'pure_gap (%)': pure_gap,
        }
        
        self.repooling_info.update({f'{repooling_methods_type}_time': repooling_time,
                                    f'{repooling_methods_type}_cost': pure_cost,
                                    f'{repooling_methods_type}_gap (%)': pure_gap,
                                    })
        
        repooling_info_df = pd.DataFrame.from_dict(repooling_info, orient='index').reset_index()
        repooling_info_df.columns = ['repooling_info', 'value']
        repooling_file_path = self.task_dir + str(self.task_id) + '_repooling_info.csv'
        header = not os.path.exists(repooling_file_path) or os.path.getsize(repooling_file_path) == 0
        repooling_info_df.to_csv(repooling_file_path, mode='a', header=header, index=False)
        

    def get_opt_pp_sol(self):
        assortment_flag = True if self.capacity_type == 'assortment' else False
        time0 = datetime.datetime.now()
        self.opt_pp_sol, _ = self.instance.solver.solve_PP_model(assortment_flag=assortment_flag)
        self.opt_pp_time = (datetime.datetime.now() - time0).total_seconds()


    def get_opt_oo_sol(self):
        time0 = datetime.datetime.now()
        assortment_flag = True if self.capacity_type == 'assortment' else False
        self.opt_oo_sol, self.exact_cost_of_oo, _ = self.instance.solver.solve_OO_model(assortment_flag=assortment_flag)
        self.opt_oo_time = (datetime.datetime.now() - time0).total_seconds()
        self.approx_cost_of_oo = self.instance.solver.solve_OO_model(assortment_flag=assortment_flag,
                                                                     given_first_plan=self.algorithm_sol_on_type)[1]
        self.approx_gap_of_algorithm = 100 * (self.approx_cost_of_oo - self.exact_cost_of_oo) / self.exact_cost_of_oo

    
    def compute_algorithm_cost(self):
        time0 = datetime.datetime.now()
        algorithm_cost, algorithm_second_stage_sol = self.instance.solver.solve_PO_model_given_first_plan(self.algorithm_sol_on_sku)
        self.algorithm_validate_time = (datetime.datetime.now() - time0).total_seconds()
        self.algorithm_cost = algorithm_cost
        self.algorithm_second_stage_sol = algorithm_second_stage_sol
        
        if self.opt_pp_sol is not None:
            opt_pp_cost, opt_pp_second_stage_sol = self.instance.solver.solve_PO_model_given_first_plan(self.opt_pp_sol)
            self.opt_pp_second_stage_sol = opt_pp_second_stage_sol
            self.opt_pp_cost = opt_pp_cost
            self.opt_pp_algorithm_gap = ((self.opt_pp_cost / self.algorithm_cost) - 1) * 100

    def compute_opt_cost(self, limited_flag=False, opt_limited_time=None):
        time0 = datetime.datetime.now()
        assortment_flag = True if self.capacity_type == 'assortment' else False
        if limited_flag:
            self.limited_po_sol, self.limited_po_cost, _ = self.instance.solver.solve_PO_model(assortment_flag=assortment_flag,
                                                                                             timelimit_flag=True)

        else:
            self.opt_po_sol, self.opt_po_cost, self.opt_po_second_stage_sol = self.instance.solver.solve_PO_model(assortment_flag=assortment_flag,
                                                                                                                  input_timilimit=opt_limited_time)

            if self.algorithm_cost is not None:
                self.ub_algo_gap = 100 * (self.algorithm_cost - self.opt_po_cost) / (self.ub_cost - self.opt_po_cost)
                self.algo_gap = 100 * (self.algorithm_cost - self.opt_po_cost) / self.opt_po_cost

            if self.opt_pp_cost is not None:
                self.ub_opt_pp_gap = 100 * (self.opt_pp_cost - self.opt_po_cost) / (self.ub_cost - self.opt_po_cost)
                self.opt_pp_gap = 100 * (self.opt_pp_cost - self.opt_po_cost) / self.opt_po_cost

            if self.limited_po_cost is not None:
                self.ub_opt_limited_po_gap = 100 * (self.limited_po_cost - self.opt_po_cost) / (self.ub_cost - self.opt_po_cost)
                self.opt_limite_po_gap = 100 * (self.limited_po_cost - self.opt_po_cost) / self.opt_po_cost
            
        if limited_flag:
            self.limited_po_time = (datetime.datetime.now() - time0).total_seconds()
        else:
            self.opt_po_time = (datetime.datetime.now() - time0).total_seconds()
        
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
            'opt_po_time': self.opt_po_time,
            'opt_oo_time': self.opt_oo_time,
            'opt_pp_time': self.opt_pp_time,
            'approx_cost_on_oo (algorithm)': self.approx_cost_of_oo,
            'exact_cost_on_oo': self.exact_cost_of_oo,
            'approx_gap (oo) %': self.approx_gap_of_algorithm,
            'algorithm_cost': self.algorithm_cost,
            'opt_po_cost': self.opt_po_cost,
            'opt_pp_cost': self.opt_pp_cost,
            'limited_po_cost': self.limited_po_cost,
            'ub_cost': self.ub_cost,
            'opt_gap (algorithm) %': self.algo_gap,
            'opt_gap (opt_pp) %': self.opt_pp_gap,
        }

        task_info.update(self.instance.pha.get_pha_paras())
        
        if self.capacity_constr:
            if self.capacity_type != 'assortment':
                task_info.update(self.instance.capmanager.CapConstr.get_algorithm_paras())
        self.task_info = task_info
    
    def write_to_file(self):
        import pickle
        task_info_df = pd.DataFrame.from_dict(self.task_info, orient='index').reset_index()
        task_info_df.columns = ['task_info', 'value']
        task_info_df.to_csv(self.task_dir + str(self.task_id) + '_task_info.csv', index=False)

        
        complete_sol_dict = {'algorithm': self.algorithm_sol_on_sku,
                    'opt_po': self.opt_po_sol,
                    'opt_oo': self.opt_oo_sol,
                    'opt_pp': self.opt_pp_sol,
                    'limited_po': self.limited_po_sol}
        pickle_path = self.task_dir + str(self.task_id) + '_solution.pkl'
        with open(pickle_path, 'wb') as f:
            pickle.dump(complete_sol_dict, f)
        