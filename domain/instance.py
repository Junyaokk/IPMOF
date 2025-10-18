from orders.source import OrderSource
from orders.branch import OrderBranch
from algorithms.solver import SOLVER
from algorithms.pha import ProgressiveHedging
from capacity.cap_manager import CapManager
import numpy as np


class Instance:
    def __init__(self, retailer_id: str, 
                 data_dir: str, 
                 random_seed: int, 
                 scenario_unit: int,
                 scenario_num: int):
        self.retailer_id = retailer_id
        self.scenario_unit = scenario_unit
        self.source = OrderSource(data_dir, scenario_unit)
        self.source.get_orders_by_region_pool(generator_type='bootstrap')
        self.random_seed = random_seed
        self.random_state = np.random.default_rng(seed=self.random_seed)
        # self.random_state = np.random.RandomState(seed=self.random_seed)
        self.scenario_num = scenario_num
        self.generation()

# class Instance:
#     def __init__(self, instance_base: InstanceBase, random_seed: int):
#         self.retailer_id = instance_base.retailer_id
#         self.scenario_unit = instance_base.scenario_unit
#         self.source = instance_base.source
#         self.random_seed = random_seed
#         self.random_state = np.random.default_rng(seed=self.random_seed)


    def generation(self):
        # generate problem instance
        self.branch = OrderBranch(source=self.source)
        self.branch.generate_random_scenario(scenario_num=self.scenario_num,
                                             random_state=self.random_state)


    def preprocess(self, capacity_constr, pha_paras={}):
        # preprocess solver and algorithms
        self.solver = SOLVER(orders=self.branch)
        pha_paras.update({'orders': self.branch, 'solver': self.solver})
        self.pha = ProgressiveHedging(**pha_paras)
        if capacity_constr:
            self.capmanager = CapManager(orders=self.branch, solver=self.solver, pha=self.pha)
        else:
            self.capmanager = None