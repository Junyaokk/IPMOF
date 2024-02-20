from orders.source import OrderSource
from orders.branch import OrderBranch
from algorithms.solver import SOLVER
from algorithms.pha import ProgressiveHedging
from capacity.cap_manager import CapManager



class InstanceBase:
    def __init__(self, retailer_id, data_dir, scenario_unit=7):
        self.retailer_id = retailer_id
        self.scenario_unit = scenario_unit
        self.source = OrderSource(data_dir, scenario_unit)
        self.source.get_orders_by_region_pool(generator_type='bootstrap')

class Instance:
    def __init__(self, instance_base: InstanceBase):
        self.retailer_id = instance_base.retailer_id
        self.scenario_unit = instance_base.scenario_unit
        self.source = instance_base.source


    def generation(self, scenario_num=100):
        # generate problem instance
        self.branch = OrderBranch(source=self.source)
        self.branch.generate_random_scenario(scenario_num=scenario_num)


    def preprocess(self, capacity_constr, pha_paras={}):
        # preprocess solver and algorithms
        self.solver = SOLVER(orders=self.branch)
        pha_paras.update({'orders': self.branch, 'solver': self.solver})
        self.pha = ProgressiveHedging(**pha_paras)
        if capacity_constr:
            self.capmanager = CapManager(orders=self.branch, solver=self.solver, pha=self.pha)
        else:
            self.capmanager = None
    


