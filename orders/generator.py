import random
from collections import Counter
from utils import *


class OrderUnit:
    def __init__(self, order_data, region_id, scenario_unit=1):
        self.data = order_data[order_data['region_id'] == region_id]
        self.sku_set = self.data['sku_id'].unique().tolist()
        self.scenario_unit = scenario_unit
        self.region_id = region_id
        self.get_origin_scenario(self.data)

    def get_origin_scenario(self, data):
        data['ord_date'] = pd.to_datetime(data['ord_date'])
        scenario_unit_timedelta = pd.to_timedelta(self.scenario_unit, unit='D')
        data['scen_id'] = ((data['ord_date'] - data['ord_date'].min()) // scenario_unit_timedelta) + 1
        data['scen_id'] = 'OS' + data['scen_id'].astype(str).str.zfill(3)
        data.drop_duplicates(subset='ord_no', inplace=True)
        scenario_order_type_dict = data.groupby('scen_id')['ord_type'].apply(list).to_dict()
        scenario_sales_dict = {k: dict(Counter(v)) for k, v in scenario_order_type_dict.items()}
        self.origin_scenario_set = list(scenario_sales_dict.keys())
        self.origin_scenario_pool = scenario_order_type_dict


class OrderGenerator:
    def __init__(self, orders: OrderUnit, generator_type):
        self.orders = orders
        self.type = generator_type
        self.generate_size_sample_pool()

    def __call__(self, scenario_num):
        if self.type == 'bootstrap':
            _scenarios = self.generate_scenario_by_bootstrap(scenario_num=scenario_num)
        else:
            raise AttributeError('Wrong generation type')
        return _scenarios


    def generate_size_sample_pool(self):
        # generate a sample pool based on order size data
        sale_by_day = self.orders.data.groupby('ord_date')['ord_no'].nunique().to_dict()
        
        if len(sale_by_day) != 0:
            sale_by_day = fill_missing_dates(sale_by_day)
            sale_set_by_day = list(sale_by_day.values())   
            if (len(sale_set_by_day) - self.orders.scenario_unit) < 20:
                sample_pool_by_unit = generate_new_data_by_sample(data=sale_set_by_day,
                                                                  sample_size=self.orders.scenario_unit,
                                                                  num_samples=100)
            else:
                sample_pool_by_unit = generate_new_data_by_sequence(data=sale_set_by_day,
                                                                    num_elements=self.orders.scenario_unit)
            
            self.size_sample_pool = sample_pool_by_unit
        else:
            self.size_sample_pool = [0]


    def generate_scenario_by_bootstrap(self, scenario_num):
        # sample order from transaction data
        def bootstrap_step(size):
            sample_set = []
            for _ in range(size):
                s_id = random.choice(self.orders.origin_scenario_set)
                order = random.choice(self.orders.origin_scenario_pool[s_id])
                sample_set.append(order)
            return sample_set
        
        # sample sales from historial data
        orders_num_by_scenario = [random.choice(self.size_sample_pool) for _ in range(scenario_num)]
        _scenarios = {}
        _type_pool = {}
        
        for index, orders_num in enumerate(orders_num_by_scenario):
            sample_set = bootstrap_step(size=orders_num)
            counter = dict(Counter(sample_set))
            order_type_sales = {k: int(v) for k, v in counter.items()}
            _scenarios['S'+str(index+1).zfill(3)] = order_type_sales
            _type_pool['S'+str(index+1).zfill(3)] = list(order_type_sales.keys())

        return _scenarios, _type_pool