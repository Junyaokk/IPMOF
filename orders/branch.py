import os
import copy
from utils import *
from orders.source import OrderSource

class OrderBranch():
    def __init__(self, source: OrderSource):
        # initialize attributes from the source
        self.rdc_id = source.rdc_id
        self.fdc_id = source.fdc_id
        self.region_set = source.region_set
        self.scenario_unit = source.scenario_unit
        self.sku_set = source.sku_set
        self.sku_weight_dict = source.sku_weight_dict
        self.first_cost_structure = source.first_cost_structure
        self.second_cost_structure = source.second_cost_structure
        self.type_collector = source.type_collector
        self.type_element_dict = source.type_element_dict
        self.type_element_pool = source.type_element_pool
        self.type_weight_dict = source.type_weight_dict
        self.type_map_dict = source.type_map_dict
        self.type_size_dict = source.type_size_dict
        self.order_type_set = source.order_type_set
        self.sku_type_map_dict = source.sku_type_map_dict
        self.sku_sales_dict = source.sku_sales_dict
        self.generator_pool = source.generator_pool
        self.order_unit_pool = source.order_unit_pool
        self.multi_item_type_set = source.multi_item_type_set
        self.single_item_type_set = source.single_item_type_set
        self.avg_involved_size_dict = source.avg_involved_size_dict
        self.sum_max_first_cost = None
        self.sum_min_second_cost = None
        self.complete_fdc_cost = None
        self.complete_rdc_cost = None
        self.UB_cost = None
        self.size = None
        self.scenario_num = None
        self.demand_diff_pairs = {}
        self.worst_second_cost = {}


    def generate_scenario_by_region(self, scenario_num=100):
        scenario_sales_by_region = {}
        scenario_pool_by_region = {}
        order_info_per_type = {}
        order_no_counter = 0
        _func_order_no = lambda n: 'ORD'+str(n).zfill(5)
        
        for region_id in self.region_set:
            generator = self.generator_pool[region_id]
            _scenarios, _type_pool = generator(scenario_num=scenario_num)
            scenario_sales_by_region[region_id] = _scenarios
            scenario_pool_by_region[region_id] = _type_pool
            
            _temp_dict = {}
            for s_id, v in _scenarios.items():
                _temp_dict[s_id] = {}
                for type_id, sales in v.items():
                    start_counter = order_no_counter + 1
                    end_counter = start_counter + sales - 1
                    _temp_dict[s_id][type_id] = [_func_order_no(i) for i in range(start_counter, end_counter+1)]
                    order_no_counter = end_counter
            order_info_per_type[region_id] = _temp_dict

        return scenario_sales_by_region, scenario_pool_by_region, order_info_per_type
    

    def aggregate_sales_by_order_type(self, scenario_sales_by_region):
        # aggregate sales by order type
        type_sales_dict = {type_id: {} for type_id in self.order_type_set}
        for type_id in self.order_type_set:
            for scenario_id in self.scenario_set:
                type_sales_dict[type_id][scenario_id] = {}
                _temp_dict = {'total': 0, 'regions': {}}
                total_sales = 0
                for region_id in self.region_set:
                    if scenario_id not in scenario_sales_by_region[region_id]:
                        continue
                    sales = scenario_sales_by_region[region_id][scenario_id].get(type_id, 0)
                    total_sales += sales
                    _temp_dict['regions'][region_id] = sales
                _temp_dict['total'] = total_sales
                type_sales_dict[type_id][scenario_id] = _temp_dict

        return type_sales_dict
    

    def construct_cost_function(self):
        # declare the functions that compute the corresponding cost
        self._func_cost_diff = lambda rdc, fc, region, u: round(max(0, 
            self.second_cost_structure[rdc][region]['intercept'] - self.second_cost_structure[fc][region]['intercept'] +
            (self.second_cost_structure[rdc][region]['coef'] - self.second_cost_structure[fc][region]['coef']) * u), 2)

        self._func_cost_rdc = lambda rdc, region, u: round(self.second_cost_structure[rdc][region]['intercept'] +
            self.second_cost_structure[rdc][region]['coef'] * u, 2)

        self._func_cost_fdc = lambda fdc, region, u: round(self.second_cost_structure[fdc][region]['intercept'] +
            self.second_cost_structure[fdc][region]['coef'] * u, 2)


    def transform_data_to_order(self, order_info_per_type):
        # transform the generated scenario to order data
        scenario_dict = {}
        order_info_per_sku_dict = {}
        sku_demand_dict = {}
        reverse_type_map_dict = {v: k for k, v in self.type_map_dict.items()}
        for s_id in self.scenario_set:
            _temp_dict = {}
            sku_demand_dict[s_id] = {i: {region_id: {} for region_id in self.region_set} for i in self.sku_set}

            for region_id in self.region_set:
                if s_id not in order_info_per_type[region_id]:
                    continue
                for type_id, order_index in order_info_per_type[region_id][s_id].items():
                    sku_set = {t[0]: t[1] for t in reverse_type_map_dict[type_id]}
                    type_weight_dict = self.type_weight_dict[type_id]

                    for sku, num in sku_set.items():
                        sku_dict = {order_no: num for order_no in order_index}
                        sku_demand_dict[s_id][sku][region_id].update(sku_dict)

        remove_empty_dict_keys(sku_demand_dict)

        order_set_per_sku_dict = {}

        return scenario_dict, order_info_per_sku_dict, order_set_per_sku_dict, sku_demand_dict


    def generate_random_scenario(self, scenario_num=100):
        self.scenario_num = scenario_num
        scenario_prob = 1 / scenario_num
        self.scenario_set = ['S'+str(index+1).zfill(3) for index in range(scenario_num)]
        self.scenario_prob_dict = {s_id: scenario_prob for s_id in self.scenario_set}
        scenario_sales_by_region, scenario_pool_by_region, order_info_per_type = self.generate_scenario_by_region(scenario_num=scenario_num)
        type_sales_dict = self.aggregate_sales_by_order_type(scenario_sales_by_region)
        scenario_dict, order_info_per_sku_dict, \
        order_set_per_sku_dict, sku_demand_dict = self.transform_data_to_order(order_info_per_type)
        demand_diff_pairs, worst_second_cost = self.get_demand_diff_pairs(type_sales_dict)

        self.scenario_sales_by_region = scenario_sales_by_region
        self.scenario_pool_by_region = scenario_pool_by_region
        self.order_info_per_type = order_info_per_type
        self.type_sales_dict = type_sales_dict
        self.scenario_dict = scenario_dict
        self.order_info_per_sku_dict = order_info_per_sku_dict
        self.sku_demand_dict = sku_demand_dict
        self.order_set_per_sku_dict = order_set_per_sku_dict
        self.demand_diff_pairs = demand_diff_pairs
        self.worst_second_cost = worst_second_cost
        self.compute_UB_cost()
        self.count_orders_size()
        self.prepare_sorted_demand_info()


    def count_orders_size(self):
        sum_order = 0
        for region_id in self.region_set:
            num = count_unique_orders(self.order_info_per_type, region_id)
            sum_order += num
        self.size = sum_order


    def get_demand_diff_pairs(self, type_sales_dict):
        # calculate type-region marginal saving cost if the order is satisfied by fdc
        rdc_id = self.rdc_id
        fdc_id = self.fdc_id
        demand_diff_pairs = {}
        worst_second_cost = {}
        sum_min_second_cost = 0
        self.construct_cost_function()

        for type_id in self.order_type_set:
            demand_diff_pairs[type_id] = {}
            worst_second_cost[type_id] = {}
            weight = self.type_weight_dict[type_id]

            for s_id in self.scenario_set:
                region_diff_dict = {}
                region_demand_dict = {}
                total_region_demand_dict = {}
                
                for region_id in self.region_set:
                    diff = self._func_cost_diff(rdc_id, fdc_id, region_id, weight)
                    demand = type_sales_dict[type_id][s_id]['regions'].get(region_id, 0)
                    
                    if demand:
                        total_region_demand_dict[region_id] = demand
                    
                    if demand and diff:
                        region_diff_dict[region_id] = diff
                        region_demand_dict[region_id] = demand
                
                sorted_diff_dict = dict(sorted(region_diff_dict.items(), key=lambda x: x[1], reverse=True))
                sorted_info_dict = {k: (v, region_demand_dict[k]) for k, v in sorted_diff_dict.items()}
                sorted_region_index = list(sorted_diff_dict.keys())
                corr_demand_set = [0] + [v[1] for v in sorted_info_dict.values()]
                corr_diff_set = [v[0] for v in sorted_info_dict.values()]
                
                if len(sorted_region_index) > 0:
                    corr_diff_set.insert(0, corr_diff_set[0])
                else:
                    corr_diff_set.insert(0, 0)
                
                worst_cost = sum([self._func_cost_rdc(rdc_id, region_id, weight) * total_region_demand_dict[region_id] 
                                  for region_id in total_region_demand_dict.keys()])
                
                sum_min_second_cost += self.scenario_prob_dict[s_id] * sum([self._func_cost_fdc(fdc_id, region_id, weight) * total_region_demand_dict[region_id] 
                                                                           for region_id in total_region_demand_dict.keys()])

                demand_diff_pairs[type_id][s_id] = {'demand': corr_demand_set, 'diff': corr_diff_set}
                worst_second_cost[type_id][s_id] = worst_cost
        
        self.sum_min_second_cost = sum_min_second_cost
        
        return demand_diff_pairs, worst_second_cost


    def compute_UB_cost(self):
        # compute prior upper bound cost
        fdc_id = self.fdc_id
        
        if self.sum_min_second_cost == 0:
            self.get_demand_diff_pairs(self.type_sales_dict)
        
        result = {}
        for s_id, s_dict in self.sku_demand_dict.items():
            result[s_id] = {}
            
            for i, i_dict in s_dict.items():
                result[s_id][i] = {}
                sum_num = 0
                
                for region_dict in i_dict.values():
                    for num in region_dict.values():
                        sum_num += num
                
                result[s_id][i] = sum_num

        supply_dict = {sku: max([result[s_id].get(sku, 0) for s_id in self.scenario_set]) for sku in self.sku_set}
        sum_max_first_cost = sum([self.sku_weight_dict[sku] * supply_dict[sku] for sku in self.sku_set]) * \
                             self.first_cost_structure[self.rdc_id][fdc_id]['coef'] + \
                             self.first_cost_structure[self.rdc_id][fdc_id]['intercept']
        self.sum_max_first_cost = sum_max_first_cost
        self.complete_fdc_cost = self.sum_max_first_cost + self.sum_min_second_cost

        self.complete_rdc_cost_dict = {s_id: 0 for s_id in self.scenario_set}
        self.complete_rdc_cost = 0
        for type_id in self.order_type_set:
            for s_id in self.scenario_set:
                self.complete_rdc_cost_dict[s_id] += self.worst_second_cost[type_id][s_id]
                self.complete_rdc_cost += self.scenario_prob_dict[s_id] * self.worst_second_cost[type_id][s_id]
        self.UB_cost = self.complete_rdc_cost
    

    def prepare_sorted_demand_info(self):
        self.union_scenarios_pool_by_region = {}
        for key1, value1 in self.scenario_pool_by_region.items():
            cur_pool = []
            for value2 in value1.values():
                cur_pool.extend(value2)
            self.union_scenarios_pool_by_region[key1] = cur_pool

        self.type_region_map = {}
        self.type_scen_demand_info = {}

        for type_id in self.order_type_set:
            weight = self.type_weight_dict[type_id]
            corr_region_set = []
            corr_diff_set = []
            for region_id in self.region_set:
                if type_id in self.union_scenarios_pool_by_region[region_id]:
                    diff = self._func_cost_diff(self.rdc_id, self.fdc_id, region_id, weight)
                    corr_region_set.append(region_id)
                    corr_diff_set.append(diff)
            if len(corr_region_set) == 0:
                continue
            self.type_region_map[type_id] = {}
            self.type_scen_demand_info[type_id] = {}
            sorted_diff_set = sorted(corr_diff_set, reverse=True)
            sorted_region_set = [x for _, x in sorted(zip(corr_diff_set, corr_region_set), 
                                                    key=lambda pair: pair[0], reverse=True)]
            sorted_diff_set.append(0)
            max_demand = 0
            self.type_region_map[type_id]['region'] = sorted_region_set
            self.type_region_map[type_id]['diff'] = sorted_diff_set
            for scen_id in self.scenario_set:
                corr_demand_set = []
                for region_id in sorted_region_set:
                    demand = self.scenario_sales_by_region[region_id][scen_id].get(type_id, 0)
                    corr_demand_set.append(demand)
                # cumulative interval
                cum_corr_demand_set = [0]
                for i in range(len(corr_demand_set)):
                    cum_corr_demand_set.append(cum_corr_demand_set[-1] + corr_demand_set[i])
                self.type_scen_demand_info[type_id][scen_id] = copy.deepcopy(cum_corr_demand_set)
                if cum_corr_demand_set[-1] > max_demand:
                    max_demand = copy.deepcopy(cum_corr_demand_set[-1])
            for scen_id in self.scenario_set:
                self.type_scen_demand_info[type_id][scen_id].append(max_demand)
        
        self.demand_lb_dict = {}
        for type_id in self.order_type_set:
            if type_id not in self.type_region_map.keys():
                continue
            else:
                self.demand_lb_dict[type_id] = {}
            region_index_set = [j+1 for j, val in enumerate(self.type_region_map[type_id]['diff']) if val != 0]
            for scen_id in self.scenario_set:
                corr_demand = sum([self.type_scen_demand_info[type_id][scen_id][region_index]
                                for region_index in region_index_set])
                self.demand_lb_dict[type_id][scen_id] = corr_demand
        
