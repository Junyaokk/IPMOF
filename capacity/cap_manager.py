import copy
from capacity.re_cap import Warehouse_CapConstr
from capacity.assmt_cap import Assortment_CapConstr
from orders.branch import OrderBranch
from algorithms.solver import SOLVER
from algorithms.pha import ProgressiveHedging


class CapManager:
    def __init__(self, 
                 orders: OrderBranch, 
                 solver: SOLVER, 
                 pha: ProgressiveHedging):
        self.orders = orders
        self.solver = solver
        self.pha = pha
        self.opt_sol = {}

    def invoke_CapConstr(self, input_capacity_val, 
                         cap_type='warehouse', 
                         algorithm_paras={}, 
                         mixed_paras={}):
        self.solver.reset_capacity_info(input_capacity_val)
        kwargs = {**algorithm_paras, 'orders': self.orders, 'solver': self.solver, 'pha': self.pha}
        cap_classes = {'warehouse': Warehouse_CapConstr, 'assortment': Assortment_CapConstr,
                       'mixed': mixed_CapConstr}
        
        if cap_type in cap_classes:
            self.CapConstr = cap_classes[cap_type](**kwargs)
            if cap_type == 'mixed':
                self.CapConstr(mixed_paras=mixed_paras)
            else:
                self.CapConstr()
            self.opt_sol_on_sku = self.CapConstr.final_sol_on_sku
            self.opt_sol_on_type = self.CapConstr.final_sol_on_type
        else:
            raise ValueError("Invalid cap_type. Supported values are 'warehouse' and 'assortment'.")
        

class mixed_CapConstr():
    def __init__(self, orders: OrderBranch, 
                 solver: SOLVER,
                 pha: ProgressiveHedging):
        self.orders = orders
        self.solver = solver
        self.pha = pha
        self.fdc_id = orders.fdc_id
        self.assortment_capconstr = Assortment_CapConstr(self.orders, 
                                                         self.solver, 
                                                         self.pha)
        
    
    def __call__(self, mixed_paras):
        self.solver.reset_mixed_info(mixed_paras)
        self.assortment_capconstr.fdc_capacity = self.solver.assortment_capacity
        self.assortment_capconstr()
        self.select_type_set = self.assortment_capconstr.A_sol_on_type
        self.select_sku_set = self.assortment_capconstr.A_sol_on_sku
        
        self.update_orders = copy.deepcopy(self.orders)
        self.update_orders.sku_set = self.select_sku_set
        self.update_orders.order_type_set = self.select_type_set
        self.update_solver = SOLVER(orders=self.update_orders)
        self.update_solver.reset_capacity_info(input_capacity_val=self.solver.warehouse_capacity)
        self.update_pha = ProgressiveHedging(orders=self.update_orders, solver=self.update_solver)
        self.warehouse_capconstr = Warehouse_CapConstr(self.update_orders,
                                                       self.update_solver,
                                                       self.update_pha)
        self.warehouse_capconstr.fdc_capacity = self.solver.warehouse_capacity
        self.warehouse_capconstr()

        self.final_sol_on_type = {type_id: {'S': 0, 'X': 0} for type_id in self.orders.order_type_set}
        self.final_sol_on_type.update(self.warehouse_capconstr.final_sol_on_type)
        
        self.final_sol_on_sku = self.pha.aggregate_sol_on_sku(self.fdc_id, 
                                                              self.final_sol_on_type)
    

    def get_algorithm_paras(self):
        return self.warehouse_capconstr.get_algorithm_paras()      