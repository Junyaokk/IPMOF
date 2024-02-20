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
                         both_paras={}):
        self.solver.reset_capacity_info(input_capacity_val)
        kwargs = {**algorithm_paras, 'orders': self.orders, 'solver': self.solver, 'pha': self.pha}
        cap_classes = {'warehouse': Warehouse_CapConstr, 'assortment': Assortment_CapConstr,
                       'both': Both_CapConstr}
        
        if cap_type in cap_classes:
            self.CapConstr = cap_classes[cap_type](**kwargs)
            if cap_type == 'both':
                self.CapConstr(both_paras=both_paras)
            else:
                self.CapConstr()
            self.opt_sol = self.CapConstr.final_sol_on_sku
        else:
            raise ValueError("Invalid cap_type. Supported values are 'warehouse' and 'assortment'.")
        

class Both_CapConstr():
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
        
    
    def __call__(self, both_paras):
        self.solver.reset_both_info(both_paras)
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

        self.opt_type_sol = {type_id: {'S': 0, 'X': 0} for type_id in self.orders.order_type_set}
        self.opt_type_sol.update(self.warehouse_capconstr.final_opt_sol)
        
        self.final_sol_on_sku = self.pha.aggregate_sol_on_sku(self.fdc_id, 
                                                              self.opt_type_sol)
    

    def get_algorithm_paras(self):
        return self.warehouse_capconstr.get_algorithm_paras()
        
        

        
    
        
        
        