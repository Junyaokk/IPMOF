import copy
import random
from utils import *
from orders.generator import OrderUnit, OrderGenerator
from default_paras import MULTT_ITEM_PROP, SINGLE_ITEM_PROP

class OrderSource:
    def __init__(self, data_dir, scenario_unit):
        self.data_dir = data_dir
        self.scenario_unit = scenario_unit
        self.order_unit_pool = {}
        self.generator_pool = {}
        self.preprocess_data()


    def preprocess_data(self):
        # load data
        sku_info_df, region_info_df, ord_info_df, trans_info_df = load_data(data_dir=self.data_dir)
        self.sku_set, self.sku_weight_dict = get_sku_info(sku_info_df)
        self.rdc_id, self.fdc_id, self.region_set = get_supply_info(region_info_df)
        self.data = ord_info_df
        self.data['ord_date'] = pd.to_datetime(self.data['ord_date'])
        self.first_cost_structure, self.second_cost_structure = get_trans_info(trans_info_df, self.rdc_id, self.fdc_id)
        
        # preprocess order types from order data
        # todo: determine corresponding select set
        self.type_collector, self.type_element_dict, self.type_weight_dict = self.order_type_separation()
        self.type_element_pool = {type_id: list(info.keys()) for type_id, info in self.type_element_dict.items()}
        self.order_type_set = list(self.type_map_dict.values())
        
        # todo: add
        self.multi_item_type_set = [type_id for type_id in self.order_type_set
                       if (len(self.type_element_pool[type_id]) != 1) 
                       or (sum(self.type_element_dict[type_id].values()) != 1)]
        self.single_item_type_set = [type_id for type_id in self.order_type_set
                            if (len(self.type_element_pool[type_id]) == 1) 
                            & (sum(self.type_element_dict[type_id].values()) == 1)]
        
        self.selected_type_set = random.sample(self.multi_item_type_set, 
                                    int(len(self.multi_item_type_set) * MULTT_ITEM_PROP)) \
                                + random.sample(self.single_item_type_set, 
                                    int(len(self.single_item_type_set) * SINGLE_ITEM_PROP))
        self.data = self.data[self.data['ord_type'].isin(self.selected_type_set)]
        self.type_element_dict = {key: val for key, val in self.type_element_dict.items() if key in self.selected_type_set}
        self.type_weight_dict = {key: val for key, val in self.type_weight_dict.items() if key in self.selected_type_set}
        self.type_element_pool = {key: val for key, val in self.type_element_pool.items() if key in self.selected_type_set}
        self.order_type_set = list(set(self.order_type_set) & set(self.selected_type_set))
        
        self.sku_type_map_dict = self.summary_sku_by_order_type()
        self.sku_sales_dict = self.data.groupby('sku_id')['sku_num'].sum().to_dict()

        self.reverse_type_map_dict = {v: k for k, v in self.type_map_dict.items()}
        self.type_size_dict = {k: sum(pool.values()) for k, pool in self.type_element_dict.items()}
        unique_data = self.data.drop_duplicates(subset=['ord_no'], inplace=False)
        self.type_sales_dict = unique_data['ord_type'].value_counts().to_dict()

        self.avg_involved_size_dict = {}
        for sku, type_set in self.sku_type_map_dict.items():
            self.avg_involved_size_dict[sku] = np.mean([self.type_size_dict.get(type_id, 0)
                                                for type_id in type_set
                                                for _ in range(self.type_sales_dict.get(type_id, 0)) 
                                                if type_id in self.type_sales_dict.keys()])


    def order_type_separation(self):
        # separate orders into unique order types and calculate type weights
        def get_order_type_id(sku_info):
            sorted_info = tuple(sorted(sku_info, key=lambda x: x[0]))
            return self.type_map_dict[sorted_info]

        # get unique order type compositions
        self.data['sku_info'] = list(zip(self.data['sku_id'], self.data['sku_num']))
        order_type_info = self.data.groupby('ord_no')['sku_info'].apply(tuple)
        type_collector = copy.deepcopy(set(order_type_info))
        unique_order_type = sorted(set([tuple(sorted(t, key=lambda x: x[0])) for t in type_collector]))
        # construct order type - id mapping
        self.type_map_dict = {t: f'TYPE{index+1:03d}' for index, t in enumerate(unique_order_type)}

        order_type_info = self.data.groupby('ord_no')['sku_info'].apply(lambda x: tuple(x)).apply(get_order_type_id)
        order_type_dict = dict(zip(order_type_info.index, order_type_info.values))
        self.data['ord_type'] = self.data['ord_no'].apply(lambda x: order_type_dict[x])
        self.data = self.data.drop(columns='sku_info')

        type_element_dict = {}
        for k, v in self.type_map_dict.items():
            inner_dict = {}
            for key in k:
                inner_dict[key[0]] = key[1]
            type_element_dict[v] = inner_dict

        # compute the weight for order types
        type_weight_dict = {}
        for type_info in unique_order_type:
            sum_weight = 0
            type_id = self.type_map_dict[type_info]
            for t in type_info:
                sum_weight += self.sku_weight_dict[t[0]] * t[1]
            # todo: revise
            # type_weight_dict[type_id] = np.ceil(sum_weight)
            type_weight_dict[type_id] = sum_weight

        return unique_order_type, type_element_dict, type_weight_dict


    def summary_sku_by_order_type(self):
        # summarize SKU by corresponding order type
        sku_type_map_dict = {sku: [] for sku in self.sku_set}
        # for type_info in self.type_collector:
        #     type_id = self.type_map_dict[type_info]
        #     for t in type_info:
        #         sku_type_map_dict[t[0]].append(type_id)
        for type_id, type_pool in self.type_element_pool.items():
            for sku in type_pool:
                sku_type_map_dict[sku].append(type_id)
        return sku_type_map_dict


    def get_orders_by_region_pool(self, generator_type='bootstrap'):
        # generate OrderUnit and OrderGenerator pools on regions
        generator_pool = {}

        for region_id in self.region_set:
            if region_id in self.order_unit_pool:
                order_unit_pool = self.order_unit_pool[region_id]
            else:
                order_unit_pool = OrderUnit(order_data=self.data, region_id=region_id, scenario_unit=self.scenario_unit)
                self.order_unit_pool[region_id] = order_unit_pool
            generator = OrderGenerator(order_unit_pool, generator_type=generator_type)
            generator_pool[region_id] = generator
        
        self.generator_pool = generator_pool